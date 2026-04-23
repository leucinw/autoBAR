#!/usr/bin/env python3
"""
Used to do parameter optimization targeting expt. HFE and neat liquid density.
Necessary settings are in the main settings.yaml file

Example usage: python parmOPT.py

- Chengwen Liu
- 2025/11

"""

import glob as globmod
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import ruamel.yaml as yaml
from scipy.optimize import least_squares

# Module-level config dict, replacing scattered global declarations.
# Populated once in main() and read by model_func / write_prm / jacobian_fd.
_config = {}

# Smallest positive lower bound we allow for any parameter — prevents the
# optimizer from driving vdw-style quantities (sigma, epsilon) through zero
# into unphysical territory.
_MIN_LOWER_BOUND = 1e-4

# kB in kcal/(mol·K)
_KB = 0.001987204

# Combined factor for density: ρ (g/cm³) = M (g/mol) / (0.6022140857 × V (Å³))
_DENSITY_FACTOR = 0.6022140857


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


def _run_autobar():
    """Invoke ``autoBAR.py auto`` under the current interpreter.

    autoBAR is blocking: by the time it returns, ``result.txt`` either
    exists (success) or was never written (failure). Raise loudly on
    either condition instead of polling for hours.
    """
    autobar_path = _config["autobar_path"]
    rc = subprocess.run([sys.executable, autobar_path, 'auto']).returncode
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


def _parse_liquid_sh(sh_file):
    """Parse the $DYNAMIC command in liquid.sh and return MD run parameters.

    Returns (dt_fs, t_out_ps, int_type, temperature_K, pressure_atm).
    """
    with open(sh_file) as f:
        content = f.read()
    for line in content.splitlines():
        if '$DYNAMIC' not in line:
            continue
        # Strip shell redirect
        line = re.sub(r'>.*$', '', line).strip()
        parts = line.split()
        parts = parts[1:]           # drop $DYNAMIC itself
        # Remove -k <keyfile> pair
        try:
            k = parts.index('-k')
            parts.pop(k)            # -k
            parts.pop(k)            # keyfile name
        except ValueError:
            pass
        # Remaining positional args (after coord base):
        # coord_base  nsteps  dt_fs  t_out_ps  int_type  temperature  [pressure]
        dt = float(parts[2])
        t_out = float(parts[3])
        int_type = parts[4]
        temperature = float(parts[5])
        pressure = float(parts[6]) if len(parts) > 6 else 1.0
        return dt, t_out, int_type, temperature, pressure
    raise ValueError(f"No $DYNAMIC line found in {sh_file}")


def _run_liquid_md(prm_file):
    """Run liquid NPT MD with prm_file parameters. Returns (arc_path, log_path).

    Runs $DYNAMIC for (n_equil + n_production) frames using the MD parameters
    parsed from liquid.sh, but with a key file that points to prm_file.
    """
    dynamic_cmd = os.environ.get('DYNAMIC')
    if not dynamic_cmd:
        raise RuntimeError("$DYNAMIC environment variable is not set")

    liquid_dir = _config["liquid_dir"]
    n_equil = _config["n_equil"]
    n_production = _config["n_production"]
    md_dt = _config["md_dt"]
    md_t_out = _config["md_t_out"]
    md_int_type = _config["md_int_type"]
    md_temperature = _config["md_temperature"]
    md_pressure = _config["md_pressure"]
    liquid_key_file = _config["liquid_key_file"]

    # steps_per_frame: convert output interval (ps) to steps
    steps_per_frame = round(md_t_out * 1000.0 / md_dt)
    n_steps = (n_equil + n_production) * steps_per_frame

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

    arc_path = os.path.join(liquid_dir, "liquid.arc")
    log_path = os.path.join(liquid_dir, "liquid-md.log")
    Path(arc_path).unlink(missing_ok=True)   # clear stale arc before fresh run

    try:
        with open(log_path, 'w') as log_f:
            rc = subprocess.run(
                [dynamic_cmd, 'liquid', '-k', tmp_key_base,
                 str(n_steps), str(md_dt), str(md_t_out),
                 md_int_type, str(md_temperature), str(md_pressure)],
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
    if not os.path.isfile(arc_path):
        raise RuntimeError("liquid.arc not produced by $DYNAMIC")
    return arc_path, log_path


def _parse_liquid_trajectory(log_file, n_equil):
    """Parse liquid-md.log → rho_frames (g/cm³) after skipping n_equil frames."""
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


def _spawn_analyze(prm_file):
    """Launch $ANALYZE as a non-blocking Popen. Returns (proc, tmp_key_path)."""
    analyze_cmd = os.environ.get('ANALYZE')
    if not analyze_cmd:
        raise RuntimeError("$ANALYZE environment variable is not set")

    liquid_arc = _config["liquid_arc"]
    arc_abs = str(Path(liquid_arc).resolve())
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
    """Compute residuals for HFE and neat liquid density."""
    param_file = _config["param_file"]
    initial_params = _config["initial_params"]
    expt_hfe = _config["expt_hfe"]
    expt_density = _config["expt_density"]
    hfe_weight = _config["hfe_weight"]
    density_weight = _config["density_weight"]
    n_equil = _config["n_equil"]

    Path('result.txt').unlink(missing_ok=True)

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

    _run_autobar()

    # --- HFE ---
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

    # --- Density: fresh MD at current params ---
    _, log_path = _run_liquid_md(current_prm)
    rho_frames = _parse_liquid_trajectory(log_path, n_equil)
    _config["rho_frames"] = rho_frames
    density = rho_frames.mean()

    print(f'Current params: {params}')
    print(f'Current HFE: {hfe:.4f} kcal/mol   Density: {density:.5f} g/cm³')

    return np.array([
        hfe_weight * (hfe - expt_hfe),
        density_weight * (density - expt_density),
    ])


def write_prm(params, fname):
    """Write the current opt line into *fname*.

    The in-tree ``param_file`` is rebuilt from the pristine snapshot before
    every write, so optimizer iterations never accumulate duplicate opt
    lines. Sidecar files (``param_file_NN``) are also derived from the
    snapshot, not from the (possibly modified) in-tree file.
    """
    param_file = _config["param_file"]
    opt_term_idx = _config["opt_term_idx"]
    snapshot = _config["param_file_snapshot"]

    if not os.path.isfile(snapshot):
        raise FileNotFoundError(f"Parameter snapshot not found: {snapshot}")

    shutil.copy(snapshot, fname)
    with open(fname, 'a') as f:
        line = (opt_term_idx.replace('-', '   ')
                + '  ' + '  '.join(str(p) for p in params) + '\n')
        f.write(line)


def jacobian_fd(params):
    """Compute Jacobian: row 0 = HFE (FD via autoBAR), row 1 = density (Eq. 4).

    All $ANALYZE calls for the density Jacobian are spawned in parallel and
    run concurrently with the autoBAR call for the HFE Jacobian.
    The trajectory used is from the preceding model_func call (current params).
    """
    param_file = _config["param_file"]
    initial_params = _config["initial_params"]
    diff_step = _config["diff_step"]
    hfe_weight = _config["hfe_weight"]
    density_weight = _config["density_weight"]
    beta = _config["beta"]
    rho_frames = _config["rho_frames"]

    n_params = len(params)
    params = np.atleast_1d(params)

    J = np.zeros((2, n_params))
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

    # Spawn all $ANALYZE processes (density Jacobian) in parallel,
    # then run autoBAR (HFE Jacobian). Both execute concurrently.
    analyze_jobs = {}   # perturb_idx → (proc, tmp_key)
    try:
        for j in range(n_params):
            for pidx in param_perturb_map[j]:
                prm_k = param_file + f'_{pidx:02d}'
                analyze_jobs[pidx] = _spawn_analyze(prm_k)

        _run_autobar()

        E_by_idx = {}
        for pidx, (proc, tmp_key) in analyze_jobs.items():
            E_by_idx[pidx] = _collect_analyze(proc, tmp_key)
        analyze_jobs.clear()
    finally:
        for pidx, (proc, tmp_key) in analyze_jobs.items():
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

    # --- Density Jacobian (row 1) via Eq. 4 ---
    for j in range(n_params):
        plus_idx, minus_idx = param_perturb_map[j]
        J[1, j] = density_weight * _density_jacobian_col(
            rho_frames,
            E_by_idx[plus_idx],
            E_by_idx[minus_idx],
            beta,
            diff_step,
        )

    return J


def main():
    # Use the object-style API so this works on ruamel.yaml >=0.18,
    # which removed the module-level yaml.load(..., Loader=...) shim.
    yaml_parser = yaml.YAML(typ='safe', pure=True)
    with open('settings.yaml') as f:
        settings = yaml_parser.load(f)

    param_file = settings["parameters"]
    expt_hfe = float(settings["expt_hfe"])
    expt_density = float(settings["expt_density"])
    hfe_weight = float(settings.get("hfe_weight", 1.0))
    density_weight = float(settings.get("density_weight", 1.0))
    liquid_dir = settings["liquid_dir"]
    liquid_key = settings.get("liquid_key", "liquid")
    n_equil = int(settings.get("n_equil", 200))
    n_production = int(settings["n_production"])
    temperature = float(settings["temperature"])
    opt_params = settings["opt_params"]
    params_range = settings["params_range"]

    # parmOPT.py lives in <repo>/utils/, so autoBAR.py is one level up.
    autobar_path = str(Path(__file__).resolve().parent.parent / 'autoBAR.py')

    # opt_params format: "vdwpair-401-402 3.8 0.05"
    s = opt_params.split()
    opt_term_idx = s[0]
    initial_params = np.atleast_1d([float(x) for x in s[1:]])
    n_params = len(initial_params)

    # form the upper and lower bound
    range_values = params_range.split()
    if len(range_values) != n_params:
        sys.exit(
            f"[Error] params_range has {len(range_values)} value(s) but "
            f"opt_params has {n_params} parameter(s); they must match."
        )
    range_values = np.array([float(v) for v in range_values])
    lb = initial_params - range_values
    ub = initial_params + range_values

    # Clamp non-positive lower bounds — letting vdw-style parameters
    # (sigma, epsilon) go through zero sends Tinker into unphysical regimes.
    bad = lb <= 0
    if bad.any():
        print(
            f"[Warning] params_range would drive lb <= 0 for indices "
            f"{np.where(bad)[0].tolist()}; clamping to {_MIN_LOWER_BOUND}."
        )
        lb = np.where(bad, _MIN_LOWER_BOUND, lb)
    bounds = (lb, ub)

    # Snapshot the user's pristine parameter file so write_prm can rebuild
    # the in-tree copy on every call without accumulating opt lines.
    snapshot = param_file + ".orig"
    if not os.path.isfile(snapshot):
        if not os.path.isfile(param_file):
            sys.exit(f"[Error] parameter file not found: {param_file}")
        shutil.copy(param_file, snapshot)

    diff_step = 0.0001
    beta = 1.0 / (_KB * temperature)

    # Total system mass from xyz + prm (for density conversion)
    liquid_xyz = str(Path(liquid_dir) / "liquid.xyz")
    total_mass = _parse_system_mass(liquid_xyz, param_file)

    # Key file template for liquid MD: liquid_dir/liquid_key.key
    liquid_key_file = str((Path(liquid_dir) / liquid_key).with_suffix('.key'))
    if not os.path.isfile(liquid_key_file):
        sys.exit(f"[Error] liquid key file not found: {liquid_key_file}")

    # Parse MD run parameters from liquid.sh
    sh_file = str(Path(liquid_dir) / "liquid.sh")
    md_dt, md_t_out, md_int_type, md_temperature, md_pressure = _parse_liquid_sh(sh_file)

    liquid_arc = str(Path(liquid_dir) / "liquid.arc")

    # Populate the module-level config dict
    _config.update({
        "param_file": param_file,
        "param_file_snapshot": snapshot,
        "expt_hfe": expt_hfe,
        "expt_density": expt_density,
        "hfe_weight": hfe_weight,
        "density_weight": density_weight,
        "opt_term_idx": opt_term_idx,
        "initial_params": initial_params,
        "autobar_path": autobar_path,
        "diff_step": diff_step,
        "beta": beta,
        "total_mass": total_mass,
        "liquid_dir": liquid_dir,
        "liquid_key_file": liquid_key_file,
        "liquid_arc": liquid_arc,
        "n_equil": n_equil,
        "n_production": n_production,
        "md_dt": md_dt,
        "md_t_out": md_t_out,
        "md_int_type": md_int_type,
        "md_temperature": md_temperature,
        "md_pressure": md_pressure,
        "rho_frames": None,
    })

    steps_per_frame = round(md_t_out * 1000.0 / md_dt)
    total_md_frames = n_equil + n_production
    total_md_steps = total_md_frames * steps_per_frame

    print("\n=== Optimization Settings ===")
    print('bounds', bounds)
    print('initial_params', initial_params)
    print('diff_step', diff_step)
    print(f'expt_hfe: {expt_hfe} kcal/mol  expt_density: {expt_density} g/cm³')
    print(f'hfe_weight: {hfe_weight}  density_weight: {density_weight}')
    print(f'temperature: {temperature} K  beta: {beta:.4f} mol/kcal')
    print(f'total_mass: {total_mass:.4f} g/mol  liquid_dir: {liquid_dir}')
    print(f'n_equil: {n_equil}  n_production: {n_production}  '
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

    print("\n=== Optimization Results ===")
    print("Success:", result.success)
    print("Message:", result.message)
    print("Optimal parameters:", result.x)
    print("Cost (sum of squared residuals):", 2 * result.cost)

if __name__ == "__main__":
    main()
