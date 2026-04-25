#!/usr/bin/env python

"""
Cluster job scheduler for ForceBalance / Tinker / PySCF workloads.

Submits jobs to Ren Lab CPU or GPU nodes via SSH, polling resource
availability and retrying until all jobs are placed.

Supported job types:
  GPU  - Tinker9 dynamics, OpenMM, GROMACS, TeraChem, PySCF (CUDA)
  CPU  - Tinker analyze, CPU dynamics

Usage examples:
  submitTinker.py -t GPU -x run_md.sh
  submitTinker.py -t CPU -c "dynamic water.xyz" -n 4
  submitTinker.py -t GPU -x calc.py -p /scratch/project1

Chengwen Liu — Feb 2022 (refactored Mar 2026)
"""

import os
import sys
import time
import logging
import argparse
import subprocess
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NODE_LIST_FILE = "/home/liuchw/bin/TinkerGPU2022/nodes.dat"
PYSCF_NODE_FILE = "/home/liuchw/bin/pyscf_node"

CPU_USAGE_LIMIT = 0.6          # only use 60 % of a node's cores
SSH_TIMEOUT = 10.0             # seconds
CPU_SETTLE_TIME = 15.0         # wait after CPU submission round
GPU_SETTLE_TIME = 30.0         # wait after GPU submission round
RETRY_INTERVAL = 5.0           # wait between retry rounds

# Processes on a GPU that indicate a card is occupied
GPU_OCCUPANT_KEYWORDS = (
    "tinker9", "dynamic", "dynamic9", "bar9",  # Tinker9
    "gmx",                                      # GROMACS
    "terachem",                                  # TeraChem
    "python",                                    # OpenMM / PySCF / etc.
)

# ---------------------------------------------------------------------------
# Logging (replaces bare ANSI prints)
# ---------------------------------------------------------------------------

GREEN = "\033[92m"
YELLOW = "\033[93m"
RED = "\033[91m"
ENDC = "\033[0m"


class ColorFormatter(logging.Formatter):
    LEVEL_COLORS = {
        logging.INFO: GREEN,
        logging.WARNING: YELLOW,
        logging.ERROR: RED,
    }

    def format(self, record):
        color = self.LEVEL_COLORS.get(record.levelno, "")
        msg = super().format(record)
        return f"{color}{msg}{ENDC}" if color else msg


log = logging.getLogger("submitTinker")
log.setLevel(logging.DEBUG)
_handler = logging.StreamHandler()
_handler.setFormatter(ColorFormatter("%(levelname)s: %(message)s"))
log.addHandler(_handler)

# ---------------------------------------------------------------------------
# SSH helper
# ---------------------------------------------------------------------------

def ssh_output(node, remote_cmd, timeout=SSH_TIMEOUT):
    """Run *remote_cmd* on *node* via SSH and return stdout lines.

    Returns an empty list (never raises) so callers can treat a
    failed / unreachable node the same as an idle one.
    """
    cmd = f'ssh -o StrictHostKeyChecking=no {node} "{remote_cmd}"'
    try:
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=timeout,
        )
        if result.returncode != 0:
            log.debug("SSH to %s returned code %d: %s", node, result.returncode, result.stderr.strip())
            return []
        return result.stdout.splitlines()
    except subprocess.TimeoutExpired:
        log.warning("SSH to %s timed out after %.0fs", node, timeout)
        return []
    except OSError as exc:
        log.warning("SSH to %s failed: %s", node, exc)
        return []

# ---------------------------------------------------------------------------
# Resource checks
# ---------------------------------------------------------------------------

def check_cpu_avail(node, nproc_required):
    """Return True if *node* has enough idle CPU cores for the job."""

    # Get total core count
    lines = ssh_output(node, "nproc")
    if not lines or not lines[0].strip().isdigit():
        log.debug("Cannot determine core count on %s — skipping", node)
        return False
    total_cores = int(lines[0].strip())
    usable_cores = int(total_cores * CPU_USAGE_LIMIT)

    # Parse `top` to estimate occupied cores
    lines = ssh_output(node, "top -n1 -b")
    occupied = 0.0
    for line in lines:
        # Match lines with R (running) or S (sleeping) state
        for state in (" R ", " S "):
            if state not in line:
                continue
            token = line.split(state)[-1].split()[0] if line.split(state)[-1].split() else "0"
            if token.replace(".", "", 1).isdigit():
                occupied += float(token) / 100.0

    occupied = round(occupied)
    available = usable_cores - occupied
    log.debug("%s: %d usable cores, ~%d occupied, %d available (need %d)",
              node, usable_cores, occupied, available, nproc_required)
    return available >= nproc_required


def get_available_gpus(node):
    """Return a list of available GPU card indices (as strings) on *node*."""

    detail_lines = ssh_output(node, "nvidia-smi -a 2>/dev/null")
    summary_lines = ssh_output(node, "nvidia-smi 2>/dev/null")

    if not detail_lines:
        log.debug("nvidia-smi unreachable on %s — skipping", node)
        return []

    # Count total GPUs
    num_gpus = 0
    for line in detail_lines:
        if "Attached GPU" in line:
            try:
                num_gpus = int(line.split()[-1])
            except ValueError:
                pass

    all_cards = [str(i) for i in range(num_gpus)]

    # Find occupied cards by scanning for known GPU-using processes.
    # nvidia-smi's per-process table puts the GPU index at parts[1]; newer
    # driver versions occasionally emit "N/A" / "MiB" there instead, so
    # require the token to be a plain integer before trusting it.
    occupied = []
    for line in summary_lines:
        if any(kw in line for kw in GPU_OCCUPANT_KEYWORDS):
            parts = line.split()
            if len(parts) >= 2 and parts[1].isdigit():
                occupied.append(parts[1])

    # Remove one slot per occupant (handles multi-GPU nodes correctly)
    available = list(all_cards)
    for card_id in occupied:
        if card_id in available:
            available.remove(card_id)

    log.debug("%s: %d GPUs total, %d occupied, %d available",
              node, num_gpus, len(occupied), len(available))
    return available

# ---------------------------------------------------------------------------
# Node list
# ---------------------------------------------------------------------------

def read_node_list(is_pyscf_job=False):
    """Read GPU and CPU node lists from config files.

    Returns (gpu_nodes, cpu_nodes).
    """
    gpu_nodes, cpu_nodes = [], []

    if not Path(NODE_LIST_FILE).is_file():
        log.error("Node list file not found: %s", NODE_LIST_FILE)
        sys.exit(1)

    with open(NODE_LIST_FILE) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            if "GPU" in line:
                gpu_nodes.append(parts[1])
            if "CPU" in line:
                cpu_nodes.append(parts[1])

    # PySCF jobs use a separate GPU node list
    if is_pyscf_job:
        if not Path(PYSCF_NODE_FILE).is_file():
            log.error("PySCF node file not found: %s", PYSCF_NODE_FILE)
            sys.exit(1)
        gpu_nodes = []
        with open(PYSCF_NODE_FILE) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                gpu_nodes.append(line)

    log.info("GPU nodes: %s", gpu_nodes or "(none)")
    log.info("CPU nodes: %s", cpu_nodes or "(none)")
    return gpu_nodes, cpu_nodes

# ---------------------------------------------------------------------------
# Job submission
# ---------------------------------------------------------------------------

def _ssh_submit(node, remote_cmd):
    """Fire-and-forget a command on *node* via SSH."""
    full_cmd = f"ssh -o StrictHostKeyChecking=no {node} '{remote_cmd}' &"
    log.info("  --> %s : %s", node, remote_cmd)
    subprocess.run(full_cmd, shell=True)


def submit_round(pending_cmds, job_type, gpu_nodes, cpu_nodes, nproc_required):
    """Try to place as many *pending_cmds* as possible in one pass.

    Returns the list of commands that could NOT be placed this round.
    """
    remaining = list(pending_cmds)

    if job_type == "CPU":
        for node in cpu_nodes:
            if not remaining:
                break
            if check_cpu_avail(node, nproc_required):
                _ssh_submit(node, remaining.pop(0))
        if len(remaining) < len(pending_cmds):
            time.sleep(CPU_SETTLE_TIME)

    else:  # GPU
        for node in gpu_nodes:
            if not remaining:
                break
            for card_id in get_available_gpus(node):
                if not remaining:
                    break
                cmd = remaining[0]
                # PySCF / Python jobs manage CUDA internally. Check any
                # whitespace-separated token for a .py suffix so a working
                # directory containing ".py" cannot trigger a false match.
                is_py_job = any(tok.endswith(".py") for tok in cmd.split())
                if is_py_job:
                    remote_cmd = cmd
                else:
                    remote_cmd = (
                        f"export CUDA_DEVICE_ORDER=PCI_BUS_ID; "
                        f'export CUDA_VISIBLE_DEVICES="{card_id}"; '
                        f"{cmd}"
                    )
                _ssh_submit(node, remote_cmd)
                remaining.pop(0)
        if len(remaining) < len(pending_cmds):
            time.sleep(GPU_SETTLE_TIME)

    return remaining

# ---------------------------------------------------------------------------
# Command building helpers
# ---------------------------------------------------------------------------

def expand_paths(paths, job_count):
    """Ensure *paths* has exactly *job_count* entries.

    - No paths given  → use cwd for every job.
    - One path given  → replicate it for every job.
    - Otherwise       → must match job_count exactly.
    """
    if len(paths) == 0:
        return [os.getcwd()] * job_count
    if len(paths) == 1:
        return paths * job_count
    if len(paths) != job_count:
        sys.exit(f"Error: got {len(paths)} paths but {job_count} jobs.")
    return list(paths)


def build_commands(jobshs, jobcmds, paths):
    """Convert script names and raw commands into final shell strings."""
    commands = []

    # Raw commands: just prepend `cd`
    for cmd, workdir in zip(jobcmds, paths[: len(jobcmds)]):
        commands.append(f"cd {workdir}; {cmd}")

    # Script files
    for script, workdir in zip(jobshs, paths[len(jobcmds):]):
        if script.endswith(".py"):
            stem = script.removesuffix(".py")
            commands.append(
                f"cd {workdir}; conda activate PySCF; "
                f"python {script} > {stem}.log 2>err"
            )
        elif script.endswith(".sh"):
            commands.append(f"cd {workdir}; bash {script}")
        else:
            sys.exit(f"Error: unsupported script type '{script}' (must be .sh or .py)")

    return commands

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Submit Tinker / PySCF jobs to Ren Lab cluster nodes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-x", dest="jobshs", nargs="+", default=[],
                        help=".sh or .py scripts to run")
    parser.add_argument("-c", dest="jobcmds", nargs="+", default=[],
                        help="Raw shell commands to run")
    parser.add_argument("-p", dest="paths", nargs="+", default=[],
                        help="Working directories (one per job, or one for all)")
    parser.add_argument("-t", dest="type", required=True, type=str.upper,
                        choices=["CPU", "GPU"], help="Job type")
    parser.add_argument("-n", dest="nproc", type=int, default=2,
                        help="CPUs required per job (default: 2)")
    parser.add_argument("-nodes", dest="nodes", nargs="+", default=[],
                        help="Override node list")
    args = parser.parse_args()

    if not args.jobshs and not args.jobcmds:
        parser.error("Provide jobs via -x (scripts) or -c (commands).")

    job_count = len(args.jobcmds) + len(args.jobshs)
    paths = expand_paths(args.paths, job_count)
    commands = build_commands(args.jobshs, args.jobcmds, paths)

    log.info("=== Submitting %d %s job(s) to Ren Lab Clusters ===", len(commands), args.type)

    # Resolve node lists
    if args.nodes:
        gpu_nodes = cpu_nodes = args.nodes
    else:
        is_pyscf = any(s.endswith(".py") for s in args.jobshs)
        gpu_nodes, cpu_nodes = read_node_list(is_pyscf_job=is_pyscf)

    # Submit loop: keep retrying until every job lands
    pending = commands
    round_num = 0
    while pending:
        round_num += 1
        if round_num > 1:
            log.info("Retry round %d — %d job(s) still pending …", round_num, len(pending))
            time.sleep(RETRY_INTERVAL)
        pending = submit_round(pending, args.type, gpu_nodes, cpu_nodes, args.nproc)

    log.info("All jobs submitted.")


if __name__ == "__main__":
    main()
