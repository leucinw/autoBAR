#!/usr/bin/env python3
"""
Used to do parameter optimization targeting any combination of expt. HFE,
neat liquid density (at one or more temperatures), and dimer interaction
energies. Each target is optional and enabled by its settings.yaml keys:
  - HFE                     : expt_hfe
  - neat liquid density     : expt_density / expt_densities
  - dimer interaction energy: dimer_data
  - dimer binding (relaxed) : dimeropt_start
At least one target must be enabled. Necessary settings are in the main
settings.yaml file.

Example usage: python parmOPT.py

- Chengwen Liu
- 2025/11

"""

import glob as globmod
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import ruamel.yaml as yaml
from scipy.optimize import least_squares

# Module-level config dict, replacing scattered global declarations.
# Populated once in main() and read by model_func / write_prm / jacobian_fd.
_config = {}

log = logging.getLogger(__name__)


def _log_step_table(step, params, hfe=None, densities=None, dimer_eints=None, dimeropt_bind=None):
    """Log a per-step summary table comparing current vs. target properties.

    Disabled targets are passed as None and skipped.
    """
    col = (24, 12, 12, 12, 16)   # property, target, current, diff, wt_norm_residual
    hdr = (f"{'Property':<{col[0]}}"
           f"{'Target':>{col[1]}}"
           f"{'Current':>{col[2]}}"
           f"{'Diff':>{col[3]}}"
           f"{'WtNormRes':>{col[4]}}")
    sep = "-" * sum(col)

    rows = []
    wt_res_list = []
    if hfe is not None:
        expt_hfe = _config["expt_hfe"]
        hfe_diff = hfe - expt_hfe
        hfe_wt_res = _config["hfe_weight"] * hfe_diff / _config["hfe_denom"]
        wt_res_list.append(hfe_wt_res)
        rows.append(
            f"{'HFE (kcal/mol)':<{col[0]}}"
            f"{expt_hfe:>{col[1]}.4f}"
            f"{hfe:>{col[2]}.4f}"
            f"{hfe_diff:>{col[3]+1}.4f}"
            f"{hfe_wt_res:>{col[4]}.4f}"
        )
    if densities is not None:
        density_denom = _config["density_denom"]
        for T, rho, rho_tgt, d_weight in zip(_config["temperatures"], densities,
                                             _config["expt_densities"],
                                             _config["density_weights"]):
            label = f"Density@{T:.1f}K (kg/m³)"
            diff = rho - rho_tgt
            wt_res = d_weight * diff / density_denom
            wt_res_list.append(wt_res)
            rows.append(
                f"{label:<{col[0]}}"
                f"{rho_tgt:>{col[1]}.3f}"
                f"{rho:>{col[2]}.3f}"
                f"{diff:>{col[3]+1}.3f}"
                f"{wt_res:>{col[4]}.4f}"
            )
    if dimer_eints is not None and _config.get("dimer_enabled"):
        dw = _config["dimer_weight"]
        dd = _config["dimer_denom"]
        for p, e in zip(_config["dimer_points"], dimer_eints):
            diff = e - p["qm"]
            wt_res = dw * p["weight"] * diff / dd
            wt_res_list.append(wt_res)
            label = f"Dimer {p['tag']} (w={p['weight']:g})"
            rows.append(
                f"{label:<{col[0]}}"
                f"{p['qm']:>{col[1]}.3f}"
                f"{e:>{col[2]}.3f}"
                f"{diff:>{col[3]+1}.3f}"
                f"{wt_res:>{col[4]}.4f}"
            )
    if dimeropt_bind is not None and _config.get("dimeropt_enabled"):
        tgt = _config["dimeropt_target"]
        diff = dimeropt_bind - tgt
        wt_res = _config["dimeropt_weight"] * diff / _config["dimeropt_denom"]
        wt_res_list.append(wt_res)
        rows.append(
            f"{'Dimer bind (opt geom)':<{col[0]}}"
            f"{tgt:>{col[1]}.3f}"
            f"{dimeropt_bind:>{col[2]}.3f}"
            f"{diff:>{col[3]+1}.3f}"
            f"{wt_res:>{col[4]}.4f}"
        )
    sum_wt_res = sum(wt_res_list)
    sum_row = (
        f"{'Sum':<{col[0]}}"
        f"{'':>{col[1]}}"
        f"{'':>{col[2]}}"
        f"{'':>{col[3]+1}}"
        f"{sum_wt_res:>{col[4]}.4f}"
    )

    opt_entries = _config["opt_entries"]
    param_parts = []
    for e in opt_entries:
        full = _reconstruct_entry_params(params, e)
        vals = [f"{v:.6g}" + ("" if is_free else "(fixed)")
                for v, is_free in zip(full, e["free_mask"])]
        param_parts.append(f"{e['term_idx']}=[{', '.join(vals)}]")
    log.info(f"--- Step {step} | {' | '.join(param_parts)} ---")
    log.info(hdr)
    log.info(sep)
    for row in rows:
        log.info(row)
    log.info(sep)
    log.info(sum_row)
    log.info(sep)


def _setup_logging(log_file="parmOPT.log"):
    """Write progress to both the console and *log_file* (overwritten each run)."""
    fmt = logging.Formatter(
        "%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh = logging.FileHandler(log_file, mode='w')
    fh.setFormatter(fmt)
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    log.setLevel(logging.DEBUG)
    log.addHandler(fh)
    log.addHandler(ch)

# Smallest positive lower bound we allow for any parameter — prevents the
# optimizer from driving vdw-style quantities (sigma, epsilon) through zero
# into unphysical territory.
_MIN_LOWER_BOUND = 1e-4

# kB in kcal/(mol·K)
_KB = 0.001987204

# Keywords present in HFE key files that must be absent from a neat-liquid key.
# Tinker keywords are case-insensitive; match at the start of a line.
_HFE_ONLY_RE = re.compile(
    r'^\s*(vdw-annihilate|vdw-lambda|ele-lambda|ligand)\b',
    re.IGNORECASE,
)

# Combined factor for density: ρ (kg/m³) = M (g/mol) / (0.0006022140857 × V (Å³))
_DENSITY_FACTOR = 0.0006022140857


def _remove_matching(pattern):
    """Remove every file or directory that matches *pattern* (glob syntax)."""
    for path in globmod.glob(pattern):
        try:
            if os.path.isdir(path):
                shutil.rmtree(path)
            else:
                os.remove(path)
        except OSError:
            pass


def _run_autobar(background=False):
    """Invoke ``autoBAR.py auto`` under the current interpreter.

    When *background* is True the Popen object is returned immediately so the
    caller can submit other jobs (e.g. neat-liquid MD) in parallel.  The
    caller must then call ``proc.wait()`` and check the return code.

    When *background* is False (default) the call blocks until autoBAR
    finishes and raises on failure or missing result.txt.
    """
    autobar_path = _config["autobar_path"]
    proc = subprocess.Popen([sys.executable, autobar_path, 'auto'])
    if background:
        return proc
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(
            f"autoBAR.py exited with code {rc}; aborting optimization"
        )
    if not os.path.isfile('result.txt'):
        raise RuntimeError(
            "autoBAR.py returned 0 but result.txt was not produced"
        )


def _restore_param_file():
    """Reset the in-tree parameter file to its pristine pre-opt snapshot."""
    shutil.copy(_config["param_file_snapshot"], _config["param_file"])


def _update_shared_key(current_prm):
    """Rewrite the shared neat-liquid key file with *current_prm* as PARAMETERS.

    All per-temperature .sh files reference the same key base name, so only
    one key file needs updating per optimizer step.
    """
    liquid_key_file = _config["liquid_key_file"]
    with open(liquid_key_file) as f:
        contents = f.read()
    prm_abs = str(Path(current_prm).resolve())
    contents = re.sub(
        r'(?im)^PARAMETERS\s+.*$',
        f'PARAMETERS   {prm_abs}',
        contents,
    )
    if not re.search(r'(?im)^archive\s*$', contents):
        contents = contents.rstrip('\n') + '\narchive\n'
    with open(liquid_key_file, 'w') as f:
        f.write(contents)


def _count_arc_frames(arc_path):
    """Return the number of trajectory frames in a Tinker .arc file, or 0 on error."""
    if not os.path.isfile(arc_path):
        return 0
    try:
        with open(arc_path, 'rb') as f:
            first_line = f.readline().decode(errors='replace').split()
            if not first_line:
                return 0
            n_atoms = int(first_line[0])
            second_line = f.readline().decode(errors='replace').split()
            if not second_line:
                return 1
            # stride = n_atoms+1 if first data line starts with "1" (no box),
            # else n_atoms+2 (with box line)
            stride = (n_atoms + 1) if second_line[0] == "1" else (n_atoms + 2)
        result = subprocess.run(
            ['wc', '-l', arc_path], capture_output=True, text=True
        )
        total_lines = int(result.stdout.split()[0])
        with open(arc_path, 'rb') as f:
            f.seek(-1, os.SEEK_END)
            if f.read(1) != b'\n':
                total_lines += 1
        return total_lines // stride
    except (OSError, ValueError, IndexError):
        return 0


def _submit_neat_liquid_to_cluster(is_initial=False):
    """Submit per-temperature neat-liquid MD .sh files to the GPU cluster.

    When *is_initial* is True, checks the existing .arc for each temperature:
      - complete (>= n_total frames): skip submission entirely.
      - partial with a .dyn checkpoint: rewrite the .sh with the remaining
        step count and append-redirect (>>) so Tinker continues seamlessly.
      - empty / no arc: submit the normal .sh for a fresh run.

    When *is_initial* is False (optimizer trial steps), the existing .arc and
    .dyn are deleted first to ensure a clean run with the new parameters.
    """
    liquid_dir   = _config["liquid_dir"]
    liquid_base  = _config["liquid_base"]
    sh_names     = _config["neat_liquid_sh_names"]
    temperatures = _config["temperatures"]
    submitexe    = _config["submitexe"]
    nodes        = _config.get("nodes", [])
    n_equil      = _config["n_equil"]
    n_production = _config["n_production"]
    md_dt        = _config["md_dt"]
    md_t_out     = _config["md_t_out"]
    steps_per_frame = round(md_t_out * 1000.0 / md_dt)
    n_total         = n_equil + n_production
    total_steps     = n_total * steps_per_frame

    to_submit = []
    for sh_name, T in zip(sh_names, temperatures):
        tag      = f"_{T:.0f}K"
        arc_path = os.path.join(liquid_dir, f"{liquid_base}{tag}.arc")
        dyn_path = os.path.join(liquid_dir, f"{liquid_base}{tag}.dyn")
        log_name = f"{liquid_base}{tag}.log"

        if not is_initial:
            # Non-initial step: clear stale arc+dyn so Tinker starts fresh
            Path(arc_path).unlink(missing_ok=True)
            Path(dyn_path).unlink(missing_ok=True)
            to_submit.append(sh_name)
            continue

        n_done = _count_arc_frames(arc_path)
        if n_done >= n_total:
            log.info("[Skip]   T=%.1fK: arc already complete (%d/%d frames)",
                     T, n_done, n_total)
            continue

        if n_done > 0 and os.path.isfile(dyn_path):
            remaining        = n_total - n_done
            remaining_steps  = remaining * steps_per_frame
            log.info("[Resume] T=%.1fK: %d/%d frames done, %d steps remaining",
                     T, n_done, n_total, remaining_steps)
            sh_path = os.path.join(liquid_dir, sh_name)
            with open(sh_path) as fh:
                sh_content = fh.read()
            # Replace total step count in the DYNAMIC9 command line
            sh_content = re.sub(
                rf'(\$DYNAMIC9\s+\S+\s+-k\s+\S+\s+){total_steps}(\s+)',
                rf'\g<1>{remaining_steps}\2',
                sh_content,
            )
            # Append log so earlier density frames are preserved
            sh_content = sh_content.replace(f' > {log_name}', f' >> {log_name}', 1)
            resume_sh_name = sh_name.replace('.sh', '-resume.sh')
            resume_sh_path = os.path.join(liquid_dir, resume_sh_name)
            with open(resume_sh_path, 'w') as fh:
                fh.write(sh_content)
            to_submit.append(resume_sh_name)
        else:
            log.info("[New]    T=%.1fK: submitting fresh MD run", T)
            to_submit.append(sh_name)

    if not to_submit:
        log.info("All neat liquid MD trajectories already complete; skipping submission.")
        return

    cmd_args = ['python', submitexe, '-x'] + to_submit + ['-t', 'GPU']
    if nodes:
        cmd_args.extend(['-nodes'] + nodes)
    cmd_args.extend(['-p', liquid_dir])
    log.info("Submitting neat liquid MD to cluster: %s", " ".join(cmd_args))
    rc = subprocess.run(cmd_args, cwd=liquid_dir).returncode
    if rc != 0:
        log.warning("submitTinker for neat liquid MD exited with code %d", rc)


def _count_log_frames(log_path):
    """Return the number of completed MD frames recorded in a Tinker log file."""
    if not os.path.isfile(log_path):
        return 0
    count = 0
    try:
        with open(log_path) as f:
            for line in f:
                if 'Frame Number' in line:
                    count += 1
    except OSError:
        pass
    return count


def _wait_for_neat_liquid_mds():
    """Block until every per-temperature neat-liquid MD arc has enough frames."""
    liquid_dir    = _config["liquid_dir"]
    liquid_base   = _config["liquid_base"]
    temperatures  = _config["temperatures"]
    n_equil       = _config["n_equil"]
    n_production  = _config["n_production"]
    check_interval = _config.get("check_interval", 60)
    n_total = n_equil + n_production

    while True:
        pending = []
        for T in temperatures:
            arc_path = os.path.join(liquid_dir, f"{liquid_base}_{T:.0f}K.arc")
            n_done = _count_arc_frames(arc_path)
            if n_done < n_total:
                pending.append(f"T={T:.1f}K ({n_done}/{n_total} frames)")
        if not pending:
            log.info("All neat liquid MD jobs completed.")
            break
        log.info("Waiting for neat liquid MD: %s", ", ".join(pending))
        time.sleep(check_interval)


# ---------------------------------------------------------------------------
# Liquid density helpers
# ---------------------------------------------------------------------------

def _parse_system_mass(xyz_file, prm_file):
    """Return total system mass (g/mol) by summing atomic masses from the prm file."""
    # Parse atom type indices from xyz (column index 5, 0-based):
    # line format: idx  symbol  x  y  z  atomtype  [bonded_atoms...]
    atom_types = []
    with open(xyz_file) as f:
        lines = f.readlines()

    # Skip header (line 0) and optional box-dimension line (line 1 if 6 floats)
    start = 1
    parts1 = lines[1].split()
    if len(parts1) >= 6:
        try:
            [float(p) for p in parts1[:6]]
            start = 2
        except ValueError:
            pass

    for line in lines[start:]:
        parts = line.split()
        if len(parts) >= 6:
            atom_types.append(int(parts[5]))

    # Parse masses from prm: AMOEBA+ atom record format:
    #   atom  type  class  symbol  "description"  atomicnum  mass  valence
    type_mass = {}
    atom_re = re.compile(
        r'^\s*atom\s+(\d+)\s+\d+\s+\S+\s+"[^"]+"\s+\d+\s+([0-9.]+)\s+\d+'
    )
    with open(prm_file) as f:
        for line in f:
            m = atom_re.match(line)
            if m:
                type_mass[int(m.group(1))] = float(m.group(2))

    total_mass = 0.0
    for atype in atom_types:
        if atype not in type_mass:
            raise KeyError(
                f"Atom type {atype} from {xyz_file} not found in {prm_file}"
            )
        total_mass += type_mass[atype]
    return total_mass


def _write_liquid_sh(liquid_dir, liquid_base, liquid_key,
                     n_equil, n_production, md_dt, md_t_out,
                     md_int_type, temperatures, md_pressure, tinkerenv):
    """Write one .sh per temperature for neat-liquid NPT MD.

    Each temperature gets its own coordinate symlink
    ({liquid_base}_{T}K.xyz → {liquid_base}.xyz) and run script, but all
    share the same key file ({liquid_key}.key) so that parmOPT only needs to
    update one key file per optimizer step.  Scripts are submitted to separate
    GPU cards so temperatures run in parallel.

    Returns a list of .sh basenames written (one per temperature).
    """
    steps_per_frame = round(md_t_out * 1000.0 / md_dt)
    total_steps = (n_equil + n_production) * steps_per_frame

    sh_names = []
    for T in temperatures:
        tag = f"_{T:.0f}K"
        xyz_name = f"{liquid_base}{tag}.xyz"
        key_name = f"{liquid_base}{tag}.key"
        sh_name  = f"{liquid_base}{tag}.sh"
        log_name = f"{liquid_base}{tag}.log"
        xyz_path = str(Path(liquid_dir) / xyz_name)
        key_path = str(Path(liquid_dir) / key_name)
        sh_path  = str(Path(liquid_dir) / sh_name)

        # Per-temperature symlinks with matching base names so Tinker writes
        # {liquid_base}{tag}.arc for each temperature independently.
        # The key symlink points to the shared key so _update_shared_key only
        # needs to rewrite one file per optimizer step.
        if not os.path.islink(xyz_path) and not os.path.isfile(xyz_path):
            os.symlink(f"{liquid_base}.xyz", xyz_path)
        if not os.path.islink(key_path):
            os.symlink(f"{liquid_key}.key", key_path)

        lines = [
            "#!/bin/bash",
            f"source {tinkerenv}",
            "# Neat-liquid NPT MD — auto-generated by parmOPT.py",
            "# Edit settings.yaml to change these parameters.",
            "#",
            f"# timestep     : {md_dt} fs",
            f"# output freq  : {md_t_out} ps  ({steps_per_frame} steps/frame)",
            f"# integrator   : {md_int_type}",
            f"# temperature  : {T} K",
            f"# pressure     : {md_pressure} atm",
            f"# n_equil      : {n_equil} frames  ({n_equil * steps_per_frame * md_dt / 1e6:.4g} ns)",
            f"# n_production : {n_production} frames  ({n_production * steps_per_frame * md_dt / 1e6:.4g} ns)",
            f"# total steps  : {total_steps}",
            "",
            f"$DYNAMIC9 {xyz_name} -k {key_name}"
            f" {total_steps} {md_dt} {md_t_out}"
            f" {md_int_type} {T} {md_pressure}"
            f" > {log_name}",
            "",
        ]
        with open(sh_path, 'w') as f:
            f.write("\n".join(lines))
        sh_names.append(sh_name)
    return sh_names


def _run_liquid_md(prm_file, temperature=None):
    """Run neat-liquid NPT MD at *temperature* with prm_file parameters.

    Returns (arc_path, log_path). When multiple temperatures are configured,
    the arc and dyn files are stored with temperature-tagged names so that
    parallel analyze jobs and dyn-file reuse work correctly across temperatures.

    Tinker names its output after the coordinate base name (``liquid_base``);
    for multi-temperature runs we rename those files immediately after each run
    so each temperature keeps its own checkpoint.
    """
    dynamic_cmd = os.environ.get('DYNAMIC9')
    if not dynamic_cmd:
        raise RuntimeError("$DYNAMIC9 environment variable is not set")

    liquid_dir = _config["liquid_dir"]
    liquid_base = _config["liquid_base"]
    n_equil = _config["n_equil"]
    n_production = _config["n_production"]
    md_dt = _config["md_dt"]
    md_t_out = _config["md_t_out"]
    md_int_type = _config["md_int_type"]
    md_pressure = _config["md_pressure"]
    liquid_key_file = _config["liquid_key_file"]
    temperatures = _config["temperatures"]

    if temperature is None:
        temperature = temperatures[0]

    # Tinker writes arc/dyn using liquid_base as the stem; for multi-temp runs
    # we save per-temperature copies so each temperature keeps its own checkpoint.
    multi_temp = len(temperatures) > 1
    arc_tinker = os.path.join(liquid_dir, f"{liquid_base}.arc")
    dyn_tinker = os.path.join(liquid_dir, f"{liquid_base}.dyn")

    if multi_temp:
        tag = f"_{temperature:.0f}K"
        arc_path = os.path.join(liquid_dir, f"{liquid_base}{tag}.arc")
        dyn_tagged = os.path.join(liquid_dir, f"{liquid_base}{tag}.dyn")
        log_path = os.path.join(liquid_dir, f"{liquid_base}{tag}.log")
        # Restore per-temperature .dyn so Tinker can continue from it, or
        # remove it to force a fresh start.
        if _config.get("use_dyn", False) and os.path.isfile(dyn_tagged):
            shutil.copy(dyn_tagged, dyn_tinker)
        else:
            Path(dyn_tinker).unlink(missing_ok=True)
    else:
        arc_path = arc_tinker
        dyn_tagged = dyn_tinker
        log_path = os.path.join(liquid_dir, f"{liquid_base}.log")
        if not _config.get("use_dyn", False):
            Path(dyn_tinker).unlink(missing_ok=True)

    Path(arc_tinker).unlink(missing_ok=True)   # clear stale arc before fresh run

    # Build a temp key file with PARAMETERS pointing to prm_file
    with open(liquid_key_file) as f:
        contents = f.read()
    prm_abs = str(Path(prm_file).resolve())
    contents = re.sub(
        r'(?im)^PARAMETERS\s+.*$',
        f'PARAMETERS   {prm_abs}',
        contents
    )
    # Place the key file in liquid_dir so $DYNAMIC can find it there
    with tempfile.NamedTemporaryFile(
        mode='w', suffix='.key', dir=liquid_dir, delete=False
    ) as tmp:
        tmp.write(contents)
        tmp_key_path = tmp.name
    tmp_key_base = Path(tmp_key_path).stem   # Tinker wants name without .key

    steps_per_frame = round(md_t_out * 1000.0 / md_dt)
    n_steps = (n_equil + n_production) * steps_per_frame

    try:
        with open(log_path, 'w') as log_f:
            rc = subprocess.run(
                [dynamic_cmd, liquid_base, '-k', tmp_key_base,
                 str(n_steps), str(md_dt), str(md_t_out),
                 md_int_type, str(temperature), str(md_pressure)],
                cwd=liquid_dir,
                stdout=log_f,
                stderr=subprocess.STDOUT,
            ).returncode
    finally:
        try:
            os.unlink(tmp_key_path)
        except OSError:
            pass

    if rc != 0:
        raise RuntimeError(f"$DYNAMIC exited with code {rc}")
    if not os.path.isfile(arc_tinker):
        raise RuntimeError(f"{liquid_base}.arc not produced by $DYNAMIC")

    if multi_temp:
        shutil.move(arc_tinker, arc_path)
        if os.path.isfile(dyn_tinker):
            shutil.copy(dyn_tinker, dyn_tagged)   # save per-temperature checkpoint

    return arc_path, log_path


def _parse_liquid_trajectory(log_file, n_equil):
    """Parse liquid MD log → rho_frames (kg/m³) after skipping n_equil frames."""
    total_mass = _config["total_mass"]
    rho_list = []
    current_lattice_a = None
    seen_potential = False

    with open(log_file) as f:
        for line in f:
            if 'Current Potential' in line:
                seen_potential = True
            elif 'Lattice Lengths' in line:
                current_lattice_a = float(line.split()[2])
            elif 'Frame Number' in line:
                if seen_potential and current_lattice_a is not None:
                    V = current_lattice_a ** 3      # Å³ (cubic box)
                    rho_list.append(total_mass / (_DENSITY_FACTOR * V))
                seen_potential = False
                current_lattice_a = None

    rho_arr = np.array(rho_list[n_equil:])
    if len(rho_arr) == 0:
        raise RuntimeError("No trajectory frames remaining after equilibration cut")
    n_production = _config["n_production"]
    if len(rho_arr) != n_production:
        raise RuntimeError(
            f"{log_file}: expected {n_production} production density frames "
            f"(n_equil={n_equil}), got {len(rho_arr)} "
            f"(total parsed={len(rho_list)}). "
            f"Check that 'Frame Number' appears after 'Current Potential' and "
            f"'Lattice Lengths' in the Tinker9 log."
        )
    return rho_arr


def _find_hfe_liquid_key(settings, autobar_path):
    """Return the path to the HFE liquid key template.

    Preference order:
    1. ``hfe_liquid_key`` in settings.yaml (explicit override)
    2. ``<autobar_repo>/dat/liquid.key`` (repo default)
    """
    explicit = settings.get("hfe_liquid_key")
    if explicit:
        if not os.path.isfile(explicit):
            raise FileNotFoundError(f"hfe_liquid_key not found: {explicit}")
        return explicit
    dat_key = str(Path(autobar_path).parent / 'dat' / 'liquid.key')
    if not os.path.isfile(dat_key):
        raise FileNotFoundError(
            f"Default HFE key not found: {dat_key}. "
            "Provide hfe_liquid_key in settings.yaml or create the neat-liquid key manually."
        )
    return dat_key


def _derive_liquid_key(src_key, dest_key):
    """Write *dest_key* from *src_key* with HFE-only keywords removed.

    Strips vdw-annihilate, vdw-lambda, ele-lambda, and ligand lines so the
    result is suitable for an unperturbed neat-liquid NPT simulation.
    """
    with open(src_key) as f:
        lines = f.readlines()
    kept = [ln for ln in lines if not _HFE_ONLY_RE.match(ln)]
    # Collapse any run of blank lines left by stripping to a single blank line
    cleaned = []
    prev_blank = False
    for ln in kept:
        is_blank = ln.strip() == ''
        if is_blank and prev_blank:
            continue
        cleaned.append(ln)
        prev_blank = is_blank
    # Ensure file ends with exactly one newline
    while cleaned and cleaned[-1].strip() == '':
        cleaned.pop()
    # archive is required for Tinker9 to write the .arc trajectory file
    if not any(ln.strip().lower() == 'archive' for ln in cleaned):
        cleaned.append('archive\n')
    else:
        cleaned.append('\n')
    os.makedirs(str(Path(dest_key).parent), exist_ok=True)
    with open(dest_key, 'w') as f:
        f.writelines(cleaned)


def _trim_arc_to_production(full_arc, prod_arc, n_equil):
    """Write a production-only arc by dropping the first n_equil frames.

    For NPT trajectories the per-frame stride is n_atoms+2 (header + box + coords);
    for NVT it is n_atoms+1.  Frames are identified by counting lines.
    """
    if n_equil == 0:
        if full_arc != prod_arc:
            shutil.copy(full_arc, prod_arc)
        return
    with open(full_arc, 'rb') as f:
        first = f.readline().split()
        if not first:
            raise RuntimeError(f"Empty arc file: {full_arc}")
        n_atoms = int(first[0])
        second = f.readline().split()
        stride = (n_atoms + 1) if second and second[0] == b'1' else (n_atoms + 2)
    skip = n_equil * stride
    with open(full_arc, 'rb') as f_in, open(prod_arc, 'wb') as f_out:
        for i, line in enumerate(f_in):
            if i >= skip:
                f_out.write(line)


def _write_analyze_sh(prm_file, pidx, temperature):
    """Write a per-(prm, temperature) analyze .sh and its named key file.

    Returns (sh_name, log_path).  Both files land in liquid_dir so that
    submitTinker.py can submit the script with -p liquid_dir.
    """
    liquid_dir       = _config["liquid_dir"]
    liquid_base      = _config["liquid_base"]
    liquid_key_file  = _config["liquid_key_file"]
    tinkerenv        = _config["tinkerenv"]

    tag      = f"_{temperature:.0f}K"
    label    = f"{liquid_base}{tag}-prm{pidx:02d}"
    prod_arc = f"{liquid_base}{tag}-prod.arc"
    key_name = f"{label}.key"
    sh_name  = f"{label}-analyze.sh"
    log_name = f"{label}-analyze.log"

    # Named key file with PARAMETERS pointing to this perturbation's prm
    with open(liquid_key_file) as fh:
        key_contents = fh.read()
    prm_abs = str(Path(prm_file).resolve())
    key_contents = re.sub(
        r'(?im)^PARAMETERS\s+.*$',
        f'PARAMETERS   {prm_abs}',
        key_contents,
    )
    with open(os.path.join(liquid_dir, key_name), 'w') as fh:
        fh.write(key_contents)

    Path(os.path.join(liquid_dir, log_name)).unlink(missing_ok=True)

    lines = [
        "#!/bin/bash",
        f"source {tinkerenv}",
        f"$ANALYZE9 {prod_arc} -k {key_name} E > {log_name}",
        "",
    ]
    with open(os.path.join(liquid_dir, sh_name), 'w') as fh:
        fh.write("\n".join(lines))

    return sh_name, os.path.join(liquid_dir, log_name)


def _submit_analyze_to_cluster(sh_names):
    """Submit analyze .sh files to the GPU cluster via submitTinker.py."""
    liquid_dir = _config["liquid_dir"]
    submitexe  = _config["submitexe"]
    nodes      = _config.get("nodes", [])
    cmd_args = ['python', submitexe, '-x'] + sh_names + ['-t', 'GPU']
    if nodes:
        cmd_args.extend(['-nodes'] + nodes)
    cmd_args.extend(['-p', liquid_dir])
    log.info("Submitting analyze jobs to cluster: %s", " ".join(cmd_args))
    rc = subprocess.run(cmd_args, cwd=liquid_dir).returncode
    if rc != 0:
        log.warning("submitTinker for analyze exited with code %d", rc)


def _wait_for_analyze_logs(log_paths):
    """Block until every analyze log has n_production energy entries."""
    n_production   = _config["n_production"]
    check_interval = _config.get("check_interval", 60)
    while True:
        pending = []
        for p in log_paths:
            if not os.path.isfile(p):
                pending.append(os.path.basename(p))
                continue
            count = sum(1 for line in open(p) if 'Total Potential Energy' in line)
            if count < n_production:
                pending.append(f"{os.path.basename(p)} ({count}/{n_production})")
        if not pending:
            log.info("All analyze jobs completed.")
            break
        log.info("Waiting for analyze (%d pending): %s",
                 len(pending), ", ".join(pending))
        time.sleep(check_interval)


def _parse_analyze_log(log_path):
    """Return a numpy array of per-frame total potential energies from an ANALYZE9 log."""
    energies = []
    with open(log_path) as fh:
        for line in fh:
            if 'Total Potential Energy' in line:
                m = re.search(r'Total Potential Energy\s*:\s*([-\d.]+)', line)
                if m:
                    energies.append(float(m.group(1)))
    if not energies:
        raise RuntimeError(f"No energy entries found in analyze log: {log_path}")
    return np.array(energies)


def _reconstruct_entry_params(free_params, entry):
    """Merge optimizer's free params with fixed values for one opt entry.

    Returns the full parameter list (free + fixed) in original order,
    suitable for writing to the parameter file.
    """
    full = list(entry["all_params"])
    fi = entry["free_start"]
    for k, is_free in enumerate(entry["free_mask"]):
        if is_free:
            full[k] = free_params[fi]
            fi += 1
    return full


def _density_jacobian_col(rho_frames, E_plus, E_minus, beta, diff_step):
    """Compute d<ρ>/dλ_j via Equation 4 of Wang et al. (2013).

        d<ρ>/dλ = -β ( <ρ · dE/dλ> - <ρ> · <dE/dλ> )

    dE/dλ is evaluated per-frame via central finite difference on the
    trajectory from the most recent liquid MD run at the current params.
    """
    dEdlambda = (E_plus - E_minus) / (2.0 * diff_step)
    rho_mean = rho_frames.mean()
    dEdl_mean = dEdlambda.mean()
    rho_dEdl_mean = (rho_frames * dEdlambda).mean()
    return -beta * (rho_dEdl_mean - rho_mean * dEdl_mean)


# ---------------------------------------------------------------------------
# Dimer interaction-energy target (optional 3rd reference)
#   E_int = E_dimer - E_mon1 - E_mon2  via Tinker AMOEBA `analyze E`,
#   compared per-point against QM reference energies (kcal/mol).
#   Cheap (gas-phase small clusters) -> evaluated directly, and its Jacobian
#   rows are obtained by central finite difference on the same analyze calls.
# ---------------------------------------------------------------------------
def _read_txyz_atoms(path):
    """Parse a Tinker .xyz -> [[idx,sym,x,y,z,type,[bonded 1-based ...]], ...]."""
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for ln in lines[1:1 + n]:
        s = ln.split()
        atoms.append([int(s[0]), s[1], float(s[2]), float(s[3]), float(s[4]),
                      int(s[5]), [int(b) for b in s[6:]]])
    return atoms


def _write_txyz_atoms(path, atoms, title="monomer"):
    with open(path, 'w') as f:
        f.write(f"{len(atoms)}  {title}\n")
        for i, a in enumerate(atoms):
            _, sym, x, y, z, typ, bonds = a
            conn = ' '.join(str(b) for b in bonds)
            f.write(f"{i+1:>3} {sym:<3} {x:12.6f} {y:12.6f} {z:12.6f} {typ:>5}  {conn}\n")


def _split_dimer_monomers(atoms, n1):
    """Split a dimer (first n1 atoms = frag1, rest = frag2); renumber bonds."""
    def build(sel):
        old2new = {old0 + 1: k + 1 for k, old0 in enumerate(sel)}   # 1-based map
        mon = []
        for old0 in sel:
            _, sym, x, y, z, typ, bonds = atoms[old0]
            nb = [old2new[b] for b in bonds if b in old2new]
            mon.append([old2new[old0 + 1], sym, x, y, z, typ, nb])
        return mon
    return build(list(range(0, n1))), build(list(range(n1, len(atoms))))


def _dimer_write_key(prm_file):
    key = os.path.join(_config["dimer_workdir"], 'dimer.key')
    with open(key, 'w') as f:
        f.write(f"parameters {os.path.abspath(prm_file)}\npolar-eps 0.00001\n")
    return key


def _analyze_dimer_energy(txyz, key):
    exe = _config["dimer_analyze"]
    cwd = _config.get("dimer_workdir") or _config.get("dimeropt_workdir")
    out = subprocess.run([exe, txyz, '-k', key, 'E'],
                         capture_output=True, text=True, cwd=cwd).stdout
    for line in out.splitlines():
        if 'Total Potential Energy' in line:
            for t in line.replace(':', ' ').split():
                try:
                    return float(t)
                except ValueError:
                    continue
    raise RuntimeError(f"analyze produced no energy for {txyz}\n{out[-400:]}")


def _dimer_eval(prm_file):
    """AMOEBA E_int = E_dimer - E_mon1 - E_mon2 (kcal/mol) per point for prm_file."""
    key = _dimer_write_key(prm_file)
    eints = []
    for p in _config["dimer_points"]:
        Ed = _analyze_dimer_energy(p["dimer"], key)
        E1 = _analyze_dimer_energy(p["mon1"], key)
        E2 = _analyze_dimer_energy(p["mon2"], key)
        eints.append(Ed - E1 - E2)
    return np.array(eints)


def _dimer_residuals(eints):
    dw, dd = _config["dimer_weight"], _config["dimer_denom"]
    return np.array([dw * p["weight"] * (e - p["qm"]) / dd
                     for p, e in zip(_config["dimer_points"], eints)])


def _dimer_setup(settings):
    """Parse the optional dimer target; pre-build monomer .xyz (geometry-only)."""
    raw = settings.get("dimer_data")
    if not raw:
        _config["dimer_enabled"] = False
        return
    if isinstance(raw, str):
        raw = [raw]
    workdir = os.path.abspath(settings.get("dimer_dir", "dimer_fit"))
    os.makedirs(workdir, exist_ok=True)
    n1 = int(settings.get("dimer_frag1_natoms", 2))
    exe = settings.get("dimer_analyze") or os.environ.get("ANALYZE8")
    if not exe:
        sys.exit("[Error] dimer target: $ANALYZE8 unset and no 'dimer_analyze' given")
    points = []
    for entry in raw:
        toks = entry.split()
        dpath = os.path.abspath(toks[0])
        qm = float(toks[1])
        w = float(toks[2]) if len(toks) > 2 else 1.0
        atoms = _read_txyz_atoms(dpath)
        mon1, mon2 = _split_dimer_monomers(atoms, n1)
        tag = os.path.splitext(os.path.basename(dpath))[0]
        d_local = os.path.join(workdir, f"{tag}_dimer.xyz")
        m1_local = os.path.join(workdir, f"{tag}_mon1.xyz")
        m2_local = os.path.join(workdir, f"{tag}_mon2.xyz")
        _write_txyz_atoms(d_local, atoms, f"{tag} dimer")
        _write_txyz_atoms(m1_local, mon1, f"{tag} mon1")
        _write_txyz_atoms(m2_local, mon2, f"{tag} mon2")
        points.append({"tag": tag, "dimer": d_local, "mon1": m1_local,
                       "mon2": m2_local, "qm": qm, "weight": w})
    qms = np.array([p["qm"] for p in points])
    denom_default = float(np.sqrt(np.mean(qms ** 2))) if len(qms) else 1.0
    _config.update({
        "dimer_enabled": True,
        "dimer_points": points,
        "dimer_weight": float(settings.get("dimer_weight", 1.0)),
        "dimer_denom": float(settings.get("dimer_denom", denom_default)),
        "dimer_workdir": workdir,
        "dimer_analyze": exe,
    })


# ---------------------------------------------------------------------------
# Dimer-OPTIMIZATION binding-energy target (optional).
#   Let AMOEBA relax the dimer with the trial params (its own geometry), then
#   E_bind = E_dimer(opt) - E_mon1 - E_mon2, compared to a QM binding energy
#   (De). Monomer intramolecular energies are independent of the fitted vdW
#   (1-2 vdW excluded; water is a different type) so they are cached once.
# ---------------------------------------------------------------------------
def _tinker_minimize(exe_min, start_xyz, key, workdir, grad):
    for stale in (start_xyz + '_2', start_xyz + '_3'):
        if os.path.exists(stale):
            os.remove(stale)
    subprocess.run([exe_min, start_xyz, '-k', key, str(grad)],
                   capture_output=True, text=True, cwd=workdir)
    opt = start_xyz + '_2'
    if not os.path.isfile(opt):
        raise RuntimeError(f"minimize produced no {opt}")
    return opt


def _dimeropt_write_key(prm_file):
    key = os.path.join(_config["dimeropt_workdir"], 'dopt.key')
    with open(key, 'w') as f:
        f.write(f"parameters {os.path.abspath(prm_file)}\npolar-eps 0.00001\n")
    return key


def _dimeropt_eval(prm_file):
    """AMOEBA-relax the dimer with prm_file and return E_bind (kcal/mol)."""
    wd = _config["dimeropt_workdir"]
    exe_min = _config["dimeropt_minimize"]
    grad = _config["dimeropt_grad"]
    key = _dimeropt_write_key(prm_file)
    d_opt = _tinker_minimize(exe_min, _config["dimeropt_start"], key, wd, grad)
    Ed = _analyze_dimer_energy(d_opt, key)   # reuse the analyze helper
    return Ed - _config["dimeropt_emon"]


def _dimeropt_residual(ebind):
    return (_config["dimeropt_weight"] * (ebind - _config["dimeropt_target"])
            / _config["dimeropt_denom"])


def _dimeropt_setup(settings):
    """Parse the optional dimer-optimization binding-energy target."""
    start = settings.get("dimeropt_start")
    if not start:
        _config["dimeropt_enabled"] = False
        return
    wd = os.path.abspath(settings.get("dimeropt_dir", "dimeropt_fit"))
    os.makedirs(wd, exist_ok=True)
    exe_min = settings.get("dimeropt_minimize") or os.environ.get("MINIMIZE8")
    exe_ana = settings.get("dimeropt_analyze") or os.environ.get("ANALYZE8")
    if not (exe_min and exe_ana):
        sys.exit("[Error] dimeropt: need $MINIMIZE8/$ANALYZE8 or explicit paths")
    _config["dimer_analyze"] = _config.get("dimer_analyze", exe_ana)  # analyze helper
    n1 = int(settings.get("dimeropt_frag1_natoms", 2))
    grad = float(settings.get("dimeropt_grad", 0.01))
    atoms = _read_txyz_atoms(os.path.abspath(start))
    mon1, mon2 = _split_dimer_monomers(atoms, n1)
    d_local = os.path.join(wd, "dopt_dimer.xyz")
    m1 = os.path.join(wd, "dopt_mon1.xyz")
    m2 = os.path.join(wd, "dopt_mon2.xyz")
    _write_txyz_atoms(d_local, atoms, "dimeropt start")
    _write_txyz_atoms(m1, mon1, "mon1 start")
    _write_txyz_atoms(m2, mon2, "mon2 start")
    _config.update({
        "dimeropt_enabled": True,
        "dimeropt_start": d_local,
        "dimeropt_workdir": wd,
        "dimeropt_minimize": exe_min,
        "dimeropt_grad": grad,
        "dimeropt_target": float(settings["dimeropt_target"]),
        "dimeropt_weight": float(settings.get("dimeropt_weight", 1.0)),
        "dimeropt_denom": float(settings.get("dimeropt_denom",
                                             np.sqrt(abs(float(settings["dimeropt_target"]))))),
    })
    # Cache the (vdW-independent) relaxed monomer energies once.
    key = _dimeropt_write_key(_config["param_file_snapshot"])
    e1 = _analyze_dimer_energy(_tinker_minimize(exe_min, m1, key, wd, grad), key)
    e2 = _analyze_dimer_energy(_tinker_minimize(exe_min, m2, key, wd, grad), key)
    _config["dimeropt_emon"] = e1 + e2


# ---------------------------------------------------------------------------
# Core optimizer functions
# ---------------------------------------------------------------------------

def model_func(params):
    """Compute residuals for every enabled target at the current params.

    Residual vector layout (disabled targets contribute no entries):
    [hfe, density_T0, density_T1, ..., dimer_points..., dimeropt]
    """
    param_file = _config["param_file"]
    initial_params = _config["initial_params"]
    hfe_enabled = _config["hfe_enabled"]
    density_enabled = _config["density_enabled"]

    if hfe_enabled:
        Path('result.txt').unlink(missing_ok=True)

    _config["step"] = _config.get("step", 0) + 1
    is_initial = np.array_equal(params, initial_params)

    if is_initial:
        write_prm(params, param_file)
        current_prm = param_file
    else:
        for i in range(1, 2 * len(params) + 2):
            _remove_matching(f'*/{param_file}_{i:02d}')
            Path(f'{param_file}_{i:02d}').unlink(missing_ok=True)
            _remove_matching(f'*/*e{100 + i * 10}*')
            _remove_matching(f'*/FEP_{i:02d}')
        write_prm(params, param_file + "_01")
        current_prm = param_file + "_01"

    # --- Submit HFE (autoBAR) and neat liquid MD in parallel ---
    # Update the shared key file with the current parameter file, then
    # launch autoBAR 'auto' in the background so we can immediately submit
    # the neat-liquid MD jobs to the GPU cluster without waiting for HFE.
    if density_enabled:
        _update_shared_key(current_prm)
    autobar_proc = _run_autobar(background=True) if hfe_enabled else None
    if density_enabled:
        _submit_neat_liquid_to_cluster(is_initial=is_initial)

    # --- Wait for HFE and parse result ---
    hfe = None
    if hfe_enabled:
        rc = autobar_proc.wait()
        if rc != 0:
            raise RuntimeError(f"autoBAR.py exited with code {rc}; aborting optimization")
        if not os.path.isfile('result.txt'):
            raise RuntimeError("autoBAR.py returned 0 but result.txt was not produced")

        with open('result.txt') as f:
            lines = f.readlines()

        fe0 = fe1 = None
        for line in lines:
            if 'SUM OF ' in line:
                fe0 = float(line.split()[-2])
            if 'FEP_001' in line:
                fe1 = float(line.split()[-1])

        if is_initial:
            if fe0 is None:
                raise RuntimeError("result.txt missing 'SUM OF ' line for reference point")
            hfe = fe0
        else:
            if fe1 is None:
                raise RuntimeError("result.txt missing 'FEP_001' line for trial point")
            hfe = fe1

    # --- Wait for neat liquid MD and parse density ---
    densities = None
    if density_enabled:
        _wait_for_neat_liquid_mds()
        liquid_dir = _config["liquid_dir"]
        liquid_base = _config["liquid_base"]
        n_equil = _config["n_equil"]
        rho_frames_list = []
        densities = []
        for T in _config["temperatures"]:
            log_path = os.path.join(liquid_dir, f"{liquid_base}_{T:.0f}K.log")
            rho_frames = _parse_liquid_trajectory(log_path, n_equil)
            rho_frames_list.append(rho_frames)
            densities.append(rho_frames.mean())

        _config["rho_frames_list"] = rho_frames_list

    # --- Dimer interaction energies at the current params (cheap, gas-phase) ---
    dimer_eints = _dimer_eval(current_prm) if _config.get("dimer_enabled") else None
    # --- Dimer-optimization binding energy (AMOEBA relaxes its own geometry) ---
    dimeropt_bind = _dimeropt_eval(current_prm) if _config.get("dimeropt_enabled") else None

    _log_step_table(_config["step"], params, hfe, densities, dimer_eints, dimeropt_bind)

    residuals = []
    if hfe is not None:
        residuals.append(
            _config["hfe_weight"] * (hfe - _config["expt_hfe"]) / _config["hfe_denom"]
        )
    if densities is not None:
        density_denom = _config["density_denom"]
        for rho, rho_tgt, d_weight in zip(densities, _config["expt_densities"],
                                          _config["density_weights"]):
            residuals.append(d_weight * (rho - rho_tgt) / density_denom)
    if dimer_eints is not None:
        residuals.extend(_dimer_residuals(dimer_eints).tolist())
    if dimeropt_bind is not None:
        residuals.append(_dimeropt_residual(dimeropt_bind))
    residuals = np.array(residuals)

    current_cost = float(np.dot(residuals, residuals))
    best_cost = _config.get("best_cost", np.inf)
    if current_cost < best_cost:
        _config["best_cost"] = current_cost
        log.info(f'Cost improved ({current_cost:.6f} < {best_cost:.6f})')
    else:
        log.info(f'Cost did not improve ({current_cost:.6f} >= {best_cost:.6f})')

    return residuals


def write_prm(params, fname):
    """Write one opt line per parameter group into *fname*.

    The in-tree ``param_file`` is rebuilt from the pristine snapshot before
    every write, so optimizer iterations never accumulate duplicate opt
    lines. Sidecar files (``param_file_NN``) are also derived from the
    snapshot, not from the (possibly modified) in-tree file.
    """
    snapshot = _config["param_file_snapshot"]
    opt_entries = _config["opt_entries"]

    if not os.path.isfile(snapshot):
        raise FileNotFoundError(f"Parameter snapshot not found: {snapshot}")

    shutil.copy(snapshot, fname)
    with open(fname, 'a') as f:
        for entry in opt_entries:
            full = _reconstruct_entry_params(params, entry)
            line = (entry["term_idx"].replace('-', '   ')
                    + '  ' + '  '.join(str(p) for p in full) + '\n')
            f.write(line)


def write_final_prm(params):
    """Write optimized parameters to ``<param_file>.final``.

    The output file is a standalone parameter file containing the snapshot
    plus the optimized opt lines — suitable for use in subsequent simulations.
    """
    param_file = _config["param_file"]
    fname = param_file + ".final"
    write_prm(params, fname)
    log.info("Optimized parameters written to %s", fname)


def jacobian_fd(params):
    """Compute the Jacobian, one row per enabled residual (same order as model_func).

    HFE row       = HFE sensitivity (FD via autoBAR FEP).
    Density rows  = density sensitivity at each temperature (Eq. 4 of
                    Wang et al. 2013, applied to per-temperature trajectories).
    Dimer rows    = central FD on the cheap gas-phase analyze calls.

    All $ANALYZE calls for the density rows are spawned in parallel across
    every (param perturbation, temperature) pair and run concurrently with
    the autoBAR call for the HFE row.
    """
    param_file = _config["param_file"]
    initial_params = _config["initial_params"]
    diff_step = _config["diff_step"]
    hfe_enabled = _config["hfe_enabled"]
    density_enabled = _config["density_enabled"]
    temperatures = _config["temperatures"]
    liquid_dir = _config["liquid_dir"]

    n_temps = len(temperatures) if density_enabled else 0
    n_hfe = 1 if hfe_enabled else 0
    n_params = len(params)
    params = np.atleast_1d(params)

    J = np.zeros((n_hfe + n_temps, n_params))
    step = diff_step * np.ones(n_params)

    if hfe_enabled:
        Path('result.txt').unlink(missing_ok=True)

    is_initial = np.array_equal(params, initial_params)
    perturb_idx = 1 if is_initial else 2

    created_indices = []
    param_perturb_map = {}   # j → (plus_idx, minus_idx)

    # Perturbed prm sidecars are consumed by autoBAR (HFE row) and by the
    # analyze jobs (density rows); neither is needed for dimer-only fits.
    if hfe_enabled or density_enabled:
        for j in range(n_params):
            params_plus = params.copy()
            params_minus = params.copy()
            dp = np.zeros_like(params)
            dp[j] = step[j]

            # plus finite difference
            lambda_str = f"{100 + perturb_idx * 10}"
            param_file_p = param_file + f'_{perturb_idx:02d}'
            _remove_matching(f'*/{param_file_p}')
            Path(param_file_p).unlink(missing_ok=True)
            _remove_matching(f'*/FEP_{perturb_idx:02d}')
            _remove_matching(f'*/*e{lambda_str}*')
            params_plus += dp
            write_prm(params_plus, param_file_p)
            created_indices.append(perturb_idx)
            plus_idx = perturb_idx
            perturb_idx += 1

            # minus finite difference
            lambda_str = f"{100 + perturb_idx * 10}"
            param_file_p = param_file + f'_{perturb_idx:02d}'
            _remove_matching(f'*/{param_file_p}')
            Path(param_file_p).unlink(missing_ok=True)
            _remove_matching(f'*/FEP_{perturb_idx:02d}')
            _remove_matching(f'*/*e{lambda_str}*')
            params_minus -= dp
            write_prm(params_minus, param_file_p)
            created_indices.append(perturb_idx)
            minus_idx = perturb_idx
            perturb_idx += 1

            param_perturb_map[j] = (plus_idx, minus_idx)

    # --- Build production arcs, write analyze .sh files, submit to cluster ---
    analyze_log_map = {}   # (pidx, temp_i) → log_path
    if density_enabled:
        liquid_base = _config["liquid_base"]
        n_equil     = _config["n_equil"]
        for T in temperatures:
            tag      = f"_{T:.0f}K"
            full_arc = os.path.join(liquid_dir, f"{liquid_base}{tag}.arc")
            prod_arc = os.path.join(liquid_dir, f"{liquid_base}{tag}-prod.arc")
            log.info("Trimming production arc for T=%.1fK (%d equil frames dropped)", T, n_equil)
            _trim_arc_to_production(full_arc, prod_arc, n_equil)
            Path(full_arc).unlink(missing_ok=True)

        sh_names = []
        for j in range(n_params):
            for pidx in param_perturb_map[j]:
                prm_k = param_file + f'_{pidx:02d}'
                for temp_i, T in enumerate(temperatures):
                    if (pidx, temp_i) not in analyze_log_map:
                        sh_name, log_path = _write_analyze_sh(prm_k, pidx, T)
                        analyze_log_map[(pidx, temp_i)] = log_path
                        sh_names.append(sh_name)
        _submit_analyze_to_cluster(sh_names)

    if hfe_enabled:
        _run_autobar()

    if density_enabled:
        _wait_for_analyze_logs(list(analyze_log_map.values()))
        E_by_pidx_temp = {
            key: _parse_analyze_log(log_path)
            for key, log_path in analyze_log_map.items()
        }

    # --- HFE Jacobian (row 0) from result.txt ---
    if hfe_enabled:
        with open('result.txt') as f:
            lines = f.readlines()

        feps = []
        for line in lines:
            if 'FEP_' in line:
                feps.append(float(line.split()[-1]))

        # Non-initial calls reuse the trial point (prm_01 → FEP_01) produced by
        # the preceding model_func call, so result.txt has one extra FEP row.
        expected = 2 * n_params if is_initial else 2 * n_params + 1
        if len(feps) != expected:
            raise RuntimeError(
                f"result.txt has {len(feps)} FEP rows, expected {expected} "
                f"(is_initial={is_initial}, created perturb_idx={created_indices})"
            )

        hfe_weight = _config["hfe_weight"]
        hfe_denom = _config["hfe_denom"]
        for j in range(n_params):
            if is_initial:
                r_plus = feps[j * 2]
                r_minus = feps[j * 2 + 1]
            else:
                r_plus = feps[j * 2 + 1]
                r_minus = feps[j * 2 + 2]
            J[0, j] = hfe_weight * (r_plus - r_minus) / (2 * diff_step) / hfe_denom

    # --- Density Jacobian rows via Eq. 4 ---
    if density_enabled:
        density_weights = _config["density_weights"]
        density_denom = _config["density_denom"]
        beta_list = _config["beta_list"]
        rho_frames_list = _config["rho_frames_list"]
        for temp_i, (d_weight, beta, rho_frames) in enumerate(
            zip(density_weights, beta_list, rho_frames_list)
        ):
            for j in range(n_params):
                plus_idx, minus_idx = param_perturb_map[j]
                J[n_hfe + temp_i, j] = d_weight * _density_jacobian_col(
                    rho_frames,
                    E_by_pidx_temp[(plus_idx, temp_i)],
                    E_by_pidx_temp[(minus_idx, temp_i)],
                    beta,
                    diff_step,
                ) / density_denom

    # --- Dimer Jacobian rows (central FD on the cheap analyze calls) ---
    if _config.get("dimer_enabled"):
        Nd = len(_config["dimer_points"])
        dw = _config["dimer_weight"]
        dd = _config["dimer_denom"]
        weights = np.array([p["weight"] for p in _config["dimer_points"]])
        Jd = np.zeros((Nd, n_params))
        prm_p = param_file + ".dimerfd_p"
        prm_m = param_file + ".dimerfd_m"
        for j in range(n_params):
            pp = params.copy(); pp[j] += diff_step
            pm = params.copy(); pm[j] -= diff_step
            write_prm(pp, prm_p)
            write_prm(pm, prm_m)
            e_p = _dimer_eval(prm_p)
            e_m = _dimer_eval(prm_m)
            Jd[:, j] = dw * weights * (e_p - e_m) / (2.0 * diff_step) / dd
        Path(prm_p).unlink(missing_ok=True)
        Path(prm_m).unlink(missing_ok=True)
        J = np.vstack([J, Jd])

    # --- Dimer-optimization binding-energy Jacobian row (central FD) ---
    if _config.get("dimeropt_enabled"):
        gw = _config["dimeropt_weight"]
        gd = _config["dimeropt_denom"]
        Jr = np.zeros((1, n_params))
        prm_p = param_file + ".doptfd_p"
        prm_m = param_file + ".doptfd_m"
        for j in range(n_params):
            pp = params.copy(); pp[j] += diff_step
            pm = params.copy(); pm[j] -= diff_step
            write_prm(pp, prm_p)
            write_prm(pm, prm_m)
            b_p = _dimeropt_eval(prm_p)
            b_m = _dimeropt_eval(prm_m)
            Jr[0, j] = gw * (b_p - b_m) / (2.0 * diff_step) / gd
        Path(prm_p).unlink(missing_ok=True)
        Path(prm_m).unlink(missing_ok=True)
        J = np.vstack([J, Jr])

    return J


def _load_tinker_env(tinkerenv):
    """Source tinker.env and merge its exported variables into os.environ."""
    result = subprocess.run(
        ['bash', '-c', f'source {tinkerenv} && env'],
        capture_output=True, text=True,
    )
    for line in result.stdout.splitlines():
        key, sep, val = line.partition('=')
        if sep:
            os.environ[key] = val


def main():
    _setup_logging()

    # Use the object-style API so this works on ruamel.yaml >=0.18,
    # which removed the module-level yaml.load(..., Loader=...) shim.
    yaml_parser = yaml.YAML(typ='safe', pure=True)
    with open('settings.yaml') as f:
        settings = yaml_parser.load(f)

    param_file = settings["parameters"]

    # --- HFE target (optional): enabled when expt_hfe is given ---
    raw_expt_hfe = settings.get("expt_hfe", None)
    hfe_enabled = raw_expt_hfe is not None
    expt_hfe = float(raw_expt_hfe) if hfe_enabled else None
    hfe_weight = float(settings.get("hfe_weight", 1.0))

    raw_opt_params = settings["opt_params"]
    if isinstance(raw_opt_params, str):
        raw_opt_params = [raw_opt_params]

    raw_params_range = settings["params_range"]
    if isinstance(raw_params_range, str):
        raw_params_range = [raw_params_range]

    if len(raw_opt_params) != len(raw_params_range):
        sys.exit(
            f"[Error] opt_params has {len(raw_opt_params)} entry(ies) but "
            f"params_range has {len(raw_params_range)}; they must match."
        )

    # --- Neat-liquid density target (optional): enabled when expt_density /
    # expt_densities is given ---
    raw_densities = settings.get("expt_densities", None)
    if raw_densities is not None:
        expt_densities = [float(d) for d in raw_densities]
    elif settings.get("expt_density", None) is not None:
        expt_densities = [float(settings["expt_density"])]
    else:
        expt_densities = []
    density_enabled = len(expt_densities) > 0

    if density_enabled:
        # Multi-temperature support: "temperatures" (list) takes priority over
        # the single-value "temperature" key for backward compatibility.
        raw_temps = settings.get("temperatures", None)
        if raw_temps is not None:
            temperatures = [float(t) for t in raw_temps]
        else:
            temperatures = [float(settings["temperature"])]

        if len(temperatures) != len(expt_densities):
            sys.exit(
                f"[Error] temperatures has {len(temperatures)} value(s) but "
                f"expt_densities has {len(expt_densities)}; they must match."
            )

        # Density weights: "density_weights" (list) or uniform "density_weight" (scalar).
        # Weights are normalized by the number of temperatures so that adding more
        # temperatures does not inflate the total density contribution to the cost.
        raw_dweights = settings.get("density_weights", None)
        if raw_dweights is not None:
            density_weights = [float(w) for w in raw_dweights]
            if len(density_weights) != len(temperatures):
                sys.exit(
                    f"[Error] density_weights has {len(density_weights)} value(s) but "
                    f"temperatures has {len(temperatures)}; they must match."
                )
        else:
            single_weight = float(settings.get("density_weight", 1.0))
            density_weights = [single_weight] * len(temperatures)

        n_temps = len(temperatures)
        density_weights = [w / n_temps for w in density_weights]

        beta_list = [1.0 / (_KB * T) for T in temperatures]
    else:
        temperatures = []
        density_weights = []
        beta_list = []
        n_temps = 0

    # --- Scale-normalization denominators (ForceBalance-style) ---
    # Defaults: std-dev of expt values when multiple points are available,
    # sqrt(|single value|) when only one point exists (matches FB's convention).
    # Both can be overridden in settings.yaml via hfe_denom / density_denom.
    density_denom = None
    if density_enabled:
        if len(expt_densities) > 1:
            density_denom_default = float(np.std(expt_densities))
        else:
            density_denom_default = float(np.sqrt(abs(expt_densities[0])))
        density_denom = float(settings.get("density_denom", density_denom_default))
        if density_denom <= 0:
            sys.exit("[Error] density_denom must be positive.")

    hfe_denom = None
    if hfe_enabled:
        hfe_denom_default = float(np.sqrt(abs(expt_hfe))) if expt_hfe != 0.0 else 1.0
        hfe_denom = float(settings.get("hfe_denom", hfe_denom_default))
        if hfe_denom <= 0:
            sys.exit("[Error] hfe_denom must be positive.")

    # parmOPT.py lives in <repo>/utils/, so autoBAR.py is one level up.
    autobar_path = str(Path(__file__).resolve().parent.parent / 'autoBAR.py')
    submitexe    = str(Path(__file__).resolve().parent / 'submitTinker.py')
    tinkerenv    = str(Path(__file__).resolve().parent.parent / 'dat' / 'tinker.env')
    _load_tinker_env(tinkerenv)
    nodes        = settings.get("node_list") or []
    check_interval = float(settings.get("checking_time", 60))

    # Parse each opt_params entry: "term_key p1 p2 ..." with a matching params_range "r1 r2 ..."
    # A range of 0 for a parameter marks it as fixed: excluded from the optimizer but
    # still written to the parameter file at its initial value.
    opt_entries = []
    all_initial = []
    all_lb = []
    all_ub = []
    free_start = 0

    for op_str, pr_str in zip(raw_opt_params, raw_params_range):
        s = op_str.split()
        term_idx = s[0]
        entry_params = np.array([float(x) for x in s[1:]])
        n = len(entry_params)
        if n == 0:
            sys.exit(f"[Error] opt_params entry '{op_str}' has no parameter values after the term key.")

        range_vals = [float(v) for v in pr_str.split()]
        if len(range_vals) != n:
            sys.exit(
                f"[Error] params_range entry '{pr_str}' has {len(range_vals)} value(s) but "
                f"opt_params entry '{op_str}' has {n} parameter(s); they must match."
            )

        free_mask = [rv != 0.0 for rv in range_vals]
        n_free = sum(free_mask)
        if n_free == 0:
            log.warning(
                f"All parameters in '{term_idx}' have range=0 and will be fixed; "
                f"the entry is written unchanged but contributes no free variables."
            )

        entry_lb = []
        entry_ub = []
        for ep, rv, is_free in zip(entry_params, range_vals, free_mask):
            if is_free:
                entry_lb.append(ep - rv)
                entry_ub.append(ep + rv)
        entry_lb = np.array(entry_lb)
        entry_ub = np.array(entry_ub)

        # Clamp non-positive lower bounds — letting vdw-style parameters
        # (sigma, epsilon) go through zero sends Tinker into unphysical regimes.
        if len(entry_lb) > 0:
            bad = entry_lb <= 0
            if bad.any():
                log.warning(
                    f"params_range for '{term_idx}' would drive lb <= 0 at free indices "
                    f"{np.where(bad)[0].tolist()}; clamping to {_MIN_LOWER_BOUND}."
                )
                entry_lb = np.where(bad, _MIN_LOWER_BOUND, entry_lb)

        opt_entries.append({
            "term_idx": term_idx,
            "all_params": entry_params.copy(),
            "free_mask": free_mask,
            "n_free": n_free,
            "free_start": free_start,
        })
        free_entry_params = entry_params[np.array(free_mask, dtype=bool)]
        all_initial.extend(free_entry_params)
        all_lb.extend(entry_lb)
        all_ub.extend(entry_ub)
        free_start += n_free

    initial_params = np.array(all_initial)
    n_params = len(initial_params)
    if n_params == 0:
        sys.exit(
            "[Error] No free parameters to optimize. "
            "Set at least one non-zero value in params_range."
        )
    lb = np.array(all_lb)
    ub = np.array(all_ub)
    bounds = (lb, ub)

    # Snapshot the user's pristine parameter file so write_prm can rebuild
    # the in-tree copy on every call without accumulating opt lines.
    snapshot = param_file + ".orig"
    if not os.path.isfile(snapshot):
        if not os.path.isfile(param_file):
            sys.exit(f"[Error] parameter file not found: {param_file}")
        shutil.copy(param_file, snapshot)

    diff_step = 0.0001

    # --- Neat-liquid MD setup (only when the density target is enabled) ---
    if density_enabled:
        liquid_dir = str(Path(settings["liquid_dir"]).resolve())
        # liquid_base: coordinate/trajectory base name inside liquid_dir.
        # Default "neat_liq" avoids collision with autoBAR's ./liquid/ directory.
        liquid_base = settings.get("liquid_base", "neat_liq")
        liquid_key = settings.get("liquid_key", liquid_base)

        # Total system mass from xyz + prm (for density conversion).
        # The xyz is the only file the user must supply in liquid_dir; the .key
        # and .sh are auto-generated below.
        liquid_xyz = str(Path(liquid_dir) / f"{liquid_base}.xyz")
        if not os.path.isfile(liquid_xyz):
            sys.exit(
                f"[Error] Neat-liquid coordinate file not found: {liquid_xyz}\n"
                f"  Place the Tinker .xyz file for the neat liquid in '{liquid_dir}/' "
                f"with the base name '{liquid_base}'.\n"
                f"  The .key and .sh files are auto-generated — only the .xyz is required."
            )
        total_mass = _parse_system_mass(liquid_xyz, param_file)

        # Key file template for neat-liquid MD: liquid_dir/<liquid_key>.key
        liquid_key_file = str((Path(liquid_dir) / liquid_key).with_suffix('.key'))
        if not os.path.isfile(liquid_key_file):
            hfe_key = _find_hfe_liquid_key(settings, autobar_path)
            _derive_liquid_key(hfe_key, liquid_key_file)
            log.info(
                f"Neat-liquid key not found; auto-generated {liquid_key_file} "
                f"from {hfe_key} (removed vdw-annihilate, vdw-lambda, ele-lambda, ligand)"
            )

        # Liquid MD parameters — shared keys with autoBAR HFE liquid settings.
        md_dt = float(settings.get("liquid_md_time_step", 2.0))
        md_t_out = float(settings.get("liquid_md_write_freq", 0.1))
        md_int_type = "4"
        md_pressure = float(settings.get("liquid_md_pressure", 1.0))
        equil_time = float(settings.get("equil_time", 0.02))        # ns
        production_time = float(settings["production_time"])         # ns
        n_equil = round(equil_time * 1000.0 / md_t_out)
        n_production = round(production_time * 1000.0 / md_t_out)

        sh_names = _write_liquid_sh(
            liquid_dir, liquid_base, liquid_key,
            n_equil, n_production, md_dt, md_t_out,
            md_int_type, temperatures, md_pressure, tinkerenv,
        )
        for sh_name in sh_names:
            log.info("Wrote liquid MD script: %s", sh_name)
    else:
        liquid_dir = None
        liquid_base = None
        liquid_key_file = None
        total_mass = None
        md_dt = md_t_out = md_pressure = None
        md_int_type = None
        n_equil = n_production = 0
        sh_names = []

    # Populate the module-level config dict
    _config.update({
        "param_file": param_file,
        "param_file_snapshot": snapshot,
        "hfe_enabled": hfe_enabled,
        "density_enabled": density_enabled,
        "expt_hfe": expt_hfe,
        "expt_densities": expt_densities,
        "hfe_weight": hfe_weight,
        "hfe_denom": hfe_denom,
        "density_weights": density_weights,
        "density_denom": density_denom,
        "temperatures": temperatures,
        "beta_list": beta_list,
        "opt_entries": opt_entries,
        "initial_params": initial_params,
        "autobar_path": autobar_path,
        "diff_step": diff_step,
        "total_mass": total_mass,
        "liquid_dir": liquid_dir,
        "liquid_base": liquid_base,
        "liquid_key_file": liquid_key_file,
        "n_equil": n_equil,
        "n_production": n_production,
        "md_dt": md_dt,
        "md_t_out": md_t_out,
        "md_int_type": md_int_type,
        "md_pressure": md_pressure,
        "rho_frames_list": None,
        "best_cost": np.inf,
        "step": 0,
        "tinkerenv": tinkerenv,
        "submitexe": submitexe,
        "nodes": nodes,
        "check_interval": check_interval,
        "neat_liquid_sh_names": sh_names,
    })

    # Optional 3rd reference: dimer interaction energies (backward compatible).
    _dimer_setup(settings)
    # Optional: dimer-optimization binding energy at AMOEBA's own geometry.
    _dimeropt_setup(settings)

    if not (hfe_enabled or density_enabled
            or _config.get("dimer_enabled") or _config.get("dimeropt_enabled")):
        sys.exit(
            "[Error] No optimization targets enabled. Provide at least one of "
            "the following in settings.yaml: expt_hfe (HFE), "
            "expt_density/expt_densities (neat liquid density), "
            "dimer_data (dimer interaction energy), "
            "dimeropt_start (dimer binding at relaxed geometry)."
        )

    log.info("=== Optimization Settings ===")
    log.info(f'diff_step {diff_step}')
    log.info(f'parameter groups ({len(opt_entries)}):')
    for entry in opt_entries:
        fs = entry["free_start"]
        nf = entry["n_free"]
        fi = 0
        parts = []
        for p, is_free in zip(entry["all_params"], entry["free_mask"]):
            if is_free:
                parts.append(f"{p:.6g} [{lb[fs+fi]:.6g}, {ub[fs+fi]:.6g}]")
                fi += 1
            else:
                parts.append(f"{p:.6g}(fixed)")
        log.info(f'  {entry["term_idx"]}: {", ".join(parts)}')
    if hfe_enabled:
        log.info(f'expt_hfe: {expt_hfe} kcal/mol  hfe_weight: {hfe_weight}  hfe_denom: {hfe_denom:.4g}')
    else:
        log.info('HFE target: disabled (no expt_hfe in settings.yaml)')
    if density_enabled:
        steps_per_frame = round(md_t_out * 1000.0 / md_dt)
        total_md_steps = (n_equil + n_production) * steps_per_frame
        log.info(f'density_denom: {density_denom:.4g} kg/m³  '
                 f'(density weights normalized by n_temps={n_temps})')
        for T, rho_tgt, d_weight in zip(temperatures, expt_densities, density_weights):
            log.info(f'  T={T:.1f} K: expt_density={rho_tgt} kg/m³  '
                     f'density_weight(effective)={d_weight:.6g}')
        log.info(f'total_mass: {total_mass:.4f} g/mol  liquid_dir: {liquid_dir}')
        log.info(f'equil: {equil_time} ns ({n_equil} frames)  '
                 f'production: {production_time} ns ({n_production} frames)  '
                 f'total MD steps per call: {total_md_steps}')
    else:
        log.info('density target: disabled (no expt_density/expt_densities in settings.yaml)')
    if _config.get("dimer_enabled"):
        log.info(f'dimer target: {len(_config["dimer_points"])} points  '
                 f'dimer_weight: {_config["dimer_weight"]}  '
                 f'dimer_denom: {_config["dimer_denom"]:.4g} kcal/mol  '
                 f'workdir: {_config["dimer_workdir"]}')
        for p in _config["dimer_points"]:
            log.info(f'  {p["tag"]}: QM E_int={p["qm"]:.3f} kcal/mol  weight={p["weight"]:g}')
    else:
        log.info('dimer target: disabled (no dimer_data in settings.yaml)')
    if _config.get("dimeropt_enabled"):
        log.info(f'dimeropt binding target: De={_config["dimeropt_target"]:.3f} kcal/mol  '
                 f'weight={_config["dimeropt_weight"]}  denom={_config["dimeropt_denom"]:.4g}  '
                 f'E_mon(cached)={_config["dimeropt_emon"]:.4f}  workdir: {_config["dimeropt_workdir"]}')
    else:
        log.info('dimeropt binding target: disabled')

    # Cleanup covers every perturb_idx the optimizer could ever create this
    # run: the model point (idx=1 on initial, idx=1 reused on trial) plus
    # 2*n_params FD windows, so 2*n_params+1 is always a safe upper bound.
    max_idx = 2 * n_params + 1
    _remove_matching(f'*/{param_file}_??')
    _remove_matching(f'{param_file}_??')
    for i in range(1, max_idx + 1):
        lambda_str = f"{100 + i * 10}"
        _remove_matching(f'*/*e{lambda_str}*')
        _remove_matching(f'*/FEP_{i:02d}')

    # Ensure the in-tree prm starts from the pristine snapshot every run,
    # even if a previous crash left it mutated.
    _restore_param_file()

    result = least_squares(
        fun=model_func,
        x0=initial_params,
        jac=jacobian_fd,
        loss='soft_l1',
        method='trf',
        verbose=2,
        bounds=bounds,
        ftol=0.0001,
        gtol=0.0001,
        xtol=0.0001,
    )

    log.info("=== Optimization Results ===")
    log.info(f"Success: {result.success}")
    log.info(f"Message: {result.message}")
    log.info(f"Optimal parameters: {result.x}")
    log.info(f"Cost (sum of squared residuals): {2 * result.cost}")

    write_final_prm(result.x)

if __name__ == "__main__":
    main()
