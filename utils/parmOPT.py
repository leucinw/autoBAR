#!/usr/bin/env python3
"""
Used to do parameter optimization targeting expt. HFE and neat liquid density
at one or more temperatures. Necessary settings are in the main settings.yaml file

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


def _log_step_table(step, params, hfe, densities):
    """Log a per-step summary table comparing current vs. target properties."""
    expt_hfe = _config["expt_hfe"]
    expt_densities = _config["expt_densities"]
    temperatures = _config["temperatures"]

    col = (24, 12, 12, 12)   # column widths: property, target, current, diff
    hdr = (f"{'Property':<{col[0]}}"
           f"{'Target':>{col[1]}}"
           f"{'Current':>{col[2]}}"
           f"{'Diff':>{col[3]}}")
    sep = "-" * sum(col)

    rows = []
    rows.append(
        f"{'HFE (kcal/mol)':<{col[0]}}"
        f"{expt_hfe:>{col[1]}.4f}"
        f"{hfe:>{col[2]}.4f}"
        f"{hfe - expt_hfe:>{col[3]+1}.4f}"
    )
    for T, rho, rho_tgt in zip(temperatures, densities, expt_densities):
        label = f"Density@{T:.1f}K (kg/m³)"
        rows.append(
            f"{label:<{col[0]}}"
            f"{rho_tgt:>{col[1]}.3f}"
            f"{rho:>{col[2]}.3f}"
            f"{rho - rho_tgt:>{col[3]+1}.3f}"
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
    with open(liquid_key_file, 'w') as f:
        f.write(contents)


def _submit_neat_liquid_to_cluster():
    """Submit per-temperature neat-liquid MD .sh files to the GPU cluster."""
    liquid_dir = _config["liquid_dir"]
    sh_names = _config["neat_liquid_sh_names"]
    submitexe = _config["submitexe"]
    nodes = _config.get("nodes", [])

    nodes_arg = f" -nodes {' '.join(nodes)}" if nodes else ""
    cmd = (f"python {submitexe} -x {' '.join(sh_names)}"
           f" -t GPU{nodes_arg} -p {liquid_dir}")
    log.info("Submitting neat liquid MD to cluster: %s", cmd)
    rc = subprocess.run(cmd, shell=True, cwd=liquid_dir).returncode
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
    """Block until every per-temperature neat-liquid MD log has enough frames."""
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
            log_path = os.path.join(liquid_dir, f"{liquid_base}_{T:.1f}K-md.log")
            n_done = _count_log_frames(log_path)
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
    main_xyz = str(Path(liquid_dir) / f"{liquid_base}.xyz")

    sh_names = []
    for T in temperatures:
        tag = f"_{T:.1f}K"
        xyz_name = f"{liquid_base}{tag}.xyz"
        sh_name  = f"{liquid_base}{tag}.sh"
        log_name = f"{liquid_base}{tag}-md.log"
        xyz_path = str(Path(liquid_dir) / xyz_name)
        sh_path  = str(Path(liquid_dir) / sh_name)

        # Symlink so each temperature has its own coordinate file; Tinker
        # then writes {xyz_name}.arc which keeps arc files separate on parallel runs.
        if not os.path.islink(xyz_path) and not os.path.isfile(xyz_path):
            os.symlink(main_xyz, xyz_path)

        lines = [
            "#!/bin/bash",
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
            f"$DYNAMIC9 {xyz_name} -k {liquid_key}"
            f" {total_steps} {md_dt} {md_t_out}"
            f" {md_int_type} {T} {md_pressure}"
            f" > {log_name}",
            "",
        ]
        with open(sh_path, 'w') as f:
            f.write(f'source {tinkerenv}\n')
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
        tag = f"_{temperature:.1f}K"
        arc_path = os.path.join(liquid_dir, f"{liquid_base}{tag}.arc")
        dyn_tagged = os.path.join(liquid_dir, f"{liquid_base}{tag}.dyn")
        log_path = os.path.join(liquid_dir, f"{liquid_base}{tag}-md.log")
        # Restore per-temperature .dyn so Tinker can continue from it, or
        # remove it to force a fresh start.
        if _config.get("use_dyn", False) and os.path.isfile(dyn_tagged):
            shutil.copy(dyn_tagged, dyn_tinker)
        else:
            Path(dyn_tinker).unlink(missing_ok=True)
    else:
        arc_path = arc_tinker
        dyn_tagged = dyn_tinker
        log_path = os.path.join(liquid_dir, f"{liquid_base}-md.log")
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
    cleaned.append('\n')
    os.makedirs(str(Path(dest_key).parent), exist_ok=True)
    with open(dest_key, 'w') as f:
        f.writelines(cleaned)


def _make_key_with_prm(prm_file, output_dir=None):
    """Create a temp key file from liquid_key_file with PARAMETERS replaced."""
    liquid_key_file = _config["liquid_key_file"]
    with open(liquid_key_file) as f:
        contents = f.read()
    prm_abs = str(Path(prm_file).resolve())
    contents = re.sub(
        r'(?im)^PARAMETERS\s+.*$',
        f'PARAMETERS   {prm_abs}',
        contents
    )
    tmp = tempfile.NamedTemporaryFile(
        mode='w', suffix='.key', dir=output_dir, delete=False
    )
    tmp.write(contents)
    tmp.close()
    return tmp.name


def _parse_analyze_energies(stdout_text, n_equil):
    """Extract per-frame total potential energies from $ANALYZE stdout."""
    energies = []
    for line in stdout_text.splitlines():
        if 'Total Potential Energy' in line:
            m = re.search(r'Total Potential Energy\s*:\s*([-\d.]+)', line)
            if m:
                energies.append(float(m.group(1)))
    if len(energies) <= n_equil:
        raise RuntimeError(
            f"$ANALYZE returned only {len(energies)} energy frames; "
            f"need more than {n_equil} (n_equil)"
        )
    return np.array(energies[n_equil:])


def _spawn_analyze(prm_file, arc_path=None):
    """Launch $ANALYZE as a non-blocking Popen. Returns (proc, tmp_key_path).

    *arc_path* overrides the trajectory file; callers must supply it when
    running at multiple temperatures (each temperature has its own arc).
    """
    analyze_cmd = os.environ.get('ANALYZE')
    if not analyze_cmd:
        raise RuntimeError("$ANALYZE environment variable is not set")

    if arc_path is None:
        raise RuntimeError(
            "_spawn_analyze: arc_path must be provided explicitly"
        )
    arc_abs = str(Path(arc_path).resolve())
    tmp_key = _make_key_with_prm(prm_file)
    proc = subprocess.Popen(
        [analyze_cmd, arc_abs, '-k', tmp_key, 'E'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    return proc, tmp_key


def _collect_analyze(proc, tmp_key):
    """Wait for a spawned analyze process and return per-frame energies."""
    n_equil = _config["n_equil"]
    try:
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(
                f"$ANALYZE failed (rc={proc.returncode}):\n{stderr}"
            )
        return _parse_analyze_energies(stdout, n_equil)
    finally:
        try:
            os.unlink(tmp_key)
        except OSError:
            pass


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
# Core optimizer functions
# ---------------------------------------------------------------------------

def model_func(params):
    """Compute residuals for HFE and neat liquid density at each temperature.

    Residual vector layout: [hfe_residual, density_residual_T0, density_residual_T1, ...]
    """
    param_file = _config["param_file"]
    initial_params = _config["initial_params"]
    expt_hfe = _config["expt_hfe"]
    expt_densities = _config["expt_densities"]
    hfe_weight = _config["hfe_weight"]
    density_weights = _config["density_weights"]
    temperatures = _config["temperatures"]
    n_equil = _config["n_equil"]
    liquid_dir = _config["liquid_dir"]
    liquid_base = _config["liquid_base"]

    Path('result.txt').unlink(missing_ok=True)

    _config["step"] = _config.get("step", 0) + 1
    is_initial = np.array_equal(params, initial_params)

    if is_initial:
        write_prm(params, param_file)
        current_prm = param_file
    else:
        _remove_matching(f'*/{param_file}_01')
        Path(f'{param_file}_01').unlink(missing_ok=True)
        _remove_matching('*/*e110*')
        _remove_matching('*/FEP_01')
        write_prm(params, param_file + "_01")
        current_prm = param_file + "_01"

    # --- Submit HFE (autoBAR) and neat liquid MD in parallel ---
    # Update the shared key file with the current parameter file, then
    # launch autoBAR 'auto' in the background so we can immediately submit
    # the neat-liquid MD jobs to the GPU cluster without waiting for HFE.
    _update_shared_key(current_prm)
    autobar_proc = _run_autobar(background=True)
    _submit_neat_liquid_to_cluster()

    # --- Wait for HFE ---
    rc = autobar_proc.wait()
    if rc != 0:
        raise RuntimeError(f"autoBAR.py exited with code {rc}; aborting optimization")
    if not os.path.isfile('result.txt'):
        raise RuntimeError("autoBAR.py returned 0 but result.txt was not produced")

    # --- Parse HFE result ---
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
    _wait_for_neat_liquid_mds()
    rho_frames_list = []
    densities = []
    for T in temperatures:
        log_path = os.path.join(liquid_dir, f"{liquid_base}_{T:.1f}K-md.log")
        rho_frames = _parse_liquid_trajectory(log_path, n_equil)
        rho_frames_list.append(rho_frames)
        densities.append(rho_frames.mean())

    _config["rho_frames_list"] = rho_frames_list

    _log_step_table(_config["step"], params, hfe, densities)

    residuals = [hfe_weight * (hfe - expt_hfe)]
    for rho, rho_tgt, d_weight in zip(densities, expt_densities, density_weights):
        residuals.append(d_weight * (rho - rho_tgt))
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


def jacobian_fd(params):
    """Compute Jacobian of shape (1 + n_temps, n_params).

    Row 0         = HFE sensitivity (FD via autoBAR FEP).
    Rows 1..n_temps = density sensitivity at each temperature (Eq. 4 of
                      Wang et al. 2013, applied to per-temperature trajectories).

    All $ANALYZE calls for the density rows are spawned in parallel across
    every (param perturbation, temperature) pair and run concurrently with
    the autoBAR call for the HFE row.
    """
    param_file = _config["param_file"]
    initial_params = _config["initial_params"]
    diff_step = _config["diff_step"]
    hfe_weight = _config["hfe_weight"]
    density_weights = _config["density_weights"]
    beta_list = _config["beta_list"]
    rho_frames_list = _config["rho_frames_list"]
    temperatures = _config["temperatures"]
    liquid_dir = _config["liquid_dir"]

    n_temps = len(temperatures)
    n_params = len(params)
    params = np.atleast_1d(params)

    J = np.zeros((1 + n_temps, n_params))
    step = diff_step * np.ones(n_params)

    Path('result.txt').unlink(missing_ok=True)

    is_initial = np.array_equal(params, initial_params)
    perturb_idx = 1 if is_initial else 2

    created_indices = []
    param_perturb_map = {}   # j → (plus_idx, minus_idx)

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

    # Spawn all $ANALYZE processes (density Jacobian) in parallel across every
    # (param perturbation, temperature) pair, then run autoBAR concurrently.
    # Key: (pidx, temp_i) → (proc, tmp_key_path)
    analyze_jobs = {}
    try:
        for j in range(n_params):
            for pidx in param_perturb_map[j]:
                prm_k = param_file + f'_{pidx:02d}'
                for temp_i, T in enumerate(temperatures):
                    liquid_base = _config["liquid_base"]
                    arc_path = os.path.join(liquid_dir, f"{liquid_base}_{T:.1f}K.arc")
                    analyze_jobs[(pidx, temp_i)] = _spawn_analyze(prm_k, arc_path=arc_path)

        _run_autobar()

        E_by_pidx_temp = {}   # (pidx, temp_i) → energies array
        for (pidx, temp_i), (proc, tmp_key) in analyze_jobs.items():
            E_by_pidx_temp[(pidx, temp_i)] = _collect_analyze(proc, tmp_key)
        analyze_jobs.clear()
    finally:
        for pidx_ti, (proc, tmp_key) in analyze_jobs.items():
            try:
                proc.kill()
            except OSError:
                pass
            try:
                os.unlink(tmp_key)
            except OSError:
                pass

    # --- HFE Jacobian (row 0) from result.txt ---
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

    for j in range(n_params):
        if is_initial:
            r_plus = feps[j * 2]
            r_minus = feps[j * 2 + 1]
        else:
            r_plus = feps[j * 2 + 1]
            r_minus = feps[j * 2 + 2]
        J[0, j] = hfe_weight * (r_plus - r_minus) / (2 * diff_step)

    # --- Density Jacobian (rows 1..n_temps) via Eq. 4 ---
    for temp_i, (d_weight, beta, rho_frames) in enumerate(
        zip(density_weights, beta_list, rho_frames_list)
    ):
        for j in range(n_params):
            plus_idx, minus_idx = param_perturb_map[j]
            J[1 + temp_i, j] = d_weight * _density_jacobian_col(
                rho_frames,
                E_by_pidx_temp[(plus_idx, temp_i)],
                E_by_pidx_temp[(minus_idx, temp_i)],
                beta,
                diff_step,
            )

    return J


def main():
    _setup_logging()

    # Use the object-style API so this works on ruamel.yaml >=0.18,
    # which removed the module-level yaml.load(..., Loader=...) shim.
    yaml_parser = yaml.YAML(typ='safe', pure=True)
    with open('settings.yaml') as f:
        settings = yaml_parser.load(f)

    param_file = settings["parameters"]
    expt_hfe = float(settings["expt_hfe"])
    hfe_weight = float(settings.get("hfe_weight", 1.0))
    liquid_dir = settings["liquid_dir"]
    # liquid_base: coordinate/trajectory base name inside liquid_dir.
    # Default "neat_liq" avoids collision with autoBAR's ./liquid/ directory.
    liquid_base = settings.get("liquid_base", "neat_liq")
    liquid_key = settings.get("liquid_key", liquid_base)
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

    # Multi-temperature support: "temperatures" (list) takes priority over
    # the single-value "temperature" key for backward compatibility.
    raw_temps = settings.get("temperatures", None)
    if raw_temps is not None:
        temperatures = [float(t) for t in raw_temps]
    else:
        temperatures = [float(settings["temperature"])]

    # Experimental densities: "expt_densities" (list) or "expt_density" (scalar).
    raw_densities = settings.get("expt_densities", None)
    if raw_densities is not None:
        expt_densities = [float(d) for d in raw_densities]
    else:
        expt_densities = [float(settings["expt_density"])]

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

    # parmOPT.py lives in <repo>/utils/, so autoBAR.py is one level up.
    autobar_path = str(Path(__file__).resolve().parent.parent / 'autoBAR.py')
    submitexe    = str(Path(__file__).resolve().parent / 'submitTinker.py')
    tinkerenv    = str(Path(__file__).resolve().parent.parent / 'dat' / 'tinker.env')
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

    # Populate the module-level config dict
    _config.update({
        "param_file": param_file,
        "param_file_snapshot": snapshot,
        "expt_hfe": expt_hfe,
        "expt_densities": expt_densities,
        "hfe_weight": hfe_weight,
        "density_weights": density_weights,
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

    steps_per_frame = round(md_t_out * 1000.0 / md_dt)
    total_md_frames = n_equil + n_production
    total_md_steps = total_md_frames * steps_per_frame

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
    log.info(f'expt_hfe: {expt_hfe} kcal/mol  hfe_weight: {hfe_weight}')
    log.info(f'density weights normalized by n_temps={n_temps} '
             f'(effective = user_weight / {n_temps})')
    for T, rho_tgt, d_weight in zip(temperatures, expt_densities, density_weights):
        log.info(f'  T={T:.1f} K: expt_density={rho_tgt} kg/m³  '
                 f'density_weight(effective)={d_weight:.6g}')
    log.info(f'total_mass: {total_mass:.4f} g/mol  liquid_dir: {liquid_dir}')
    log.info(f'equil: {equil_time} ns ({n_equil} frames)  '
             f'production: {production_time} ns ({n_production} frames)  '
             f'total MD steps per call: {total_md_steps}')

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

if __name__ == "__main__":
    main()
