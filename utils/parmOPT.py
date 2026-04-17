#!/usr/bin/env python3
"""
Used to do parameter optimization targeting expt. HFE
Necessary settings are in the main settings.yaml file

Example usage: python parmOPT.py

- Chengwen Liu
- 2025/11

"""

import glob as globmod
import os
import shutil
import subprocess
import sys
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


def model_func(params):
    """
    Compute residuals (difference between predicted and experimental values).
    """
    param_file = _config["param_file"]
    initial_params = _config["initial_params"]
    expt_hfe = _config["expt_hfe"]

    Path('result.txt').unlink(missing_ok=True)

    if np.array_equal(params, initial_params):
        # Reference point: rebuild the in-tree prm from the pristine snapshot
        # plus one opt line. Never append to an already-modified file.
        write_prm(params, param_file)
    else:
        _remove_matching(f'*/{param_file}_01')
        Path(f'{param_file}_01').unlink(missing_ok=True)
        _remove_matching('*/*e110*')
        _remove_matching('*/FEP_01')
        write_prm(params, param_file + "_01")

    _run_autobar()

    # read result.txt file
    with open('result.txt') as f:
        lines = f.readlines()

    fe0 = None
    fe1 = None
    for line in lines:
        if 'SUM OF ' in line:
            fe0 = float(line.split()[-2])
        if 'FEP_001' in line:
            fe1 = float(line.split()[-1])

    print(f'Current params: {params}')
    if np.array_equal(params, initial_params):
        if fe0 is None:
            raise RuntimeError("result.txt missing 'SUM OF ' line for reference point")
        print(f'Current HFE: {fe0}')
        return np.atleast_1d([fe0 - expt_hfe])
    else:
        if fe1 is None:
            raise RuntimeError("result.txt missing 'FEP_001' line for trial point")
        print(f'Current HFE: {fe1}')
        return np.atleast_1d([fe1 - expt_hfe])


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

    # Rebuild the target file from the pristine snapshot every time.
    shutil.copy(snapshot, fname)
    with open(fname, 'a') as f:
        line = (opt_term_idx.replace('-', '   ')
                + '  ' + '  '.join(str(p) for p in params) + '\n')
        f.write(line)


def jacobian_fd(params):
    """
    Compute Jacobian matrix numerically using central finite difference.
    The dynamic/bar simulations are running parallelly so there will be
    `2*len(params)` speedup comparing to just have model_func with jac='2-point'.
    """
    param_file = _config["param_file"]
    initial_params = _config["initial_params"]
    diff_step = _config["diff_step"]

    n_params = len(params)
    params = np.atleast_1d(params)

    # only ONE HFE data
    J = np.zeros((1, n_params))
    step = diff_step * np.ones(n_params)

    Path('result.txt').unlink(missing_ok=True)

    is_initial = np.array_equal(params, initial_params)
    perturb_idx = 1 if is_initial else 2

    # Perturbation windows we create this call — used below to sanity-check
    # the FEP row count in result.txt.
    created_indices = []

    for j in range(n_params):

        # deep copy the numpy array
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
        perturb_idx += 1

    _run_autobar()

    # read result.txt file
    with open('result.txt') as f:
        lines = f.readlines()

    feps = []
    for line in lines:
        if 'FEP_' in line:
            fep = float(line.split()[-1])
            feps.append(fep)

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
        J[:, j] = (r_plus - r_minus) / (2 * diff_step)

    return J


def main():
    # Use the object-style API so this works on ruamel.yaml >=0.18,
    # which removed the module-level yaml.load(..., Loader=...) shim.
    yaml_parser = yaml.YAML(typ='safe', pure=True)
    with open('settings.yaml') as f:
        settings = yaml_parser.load(f)

    param_file = settings["parameters"]
    expt_hfe = float(settings["expt_hfe"])
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
    # An existing snapshot is authoritative — treat it as the base going
    # forward (the in-tree file may already carry an opt line from a
    # previous run).

    # x*diff_step — only used inside jacobian_fd; NB least_squares ignores
    # its own diff_step when a callable `jac` is supplied.
    diff_step = 0.0001

    # Populate the module-level config dict (replaces scattered globals)
    _config.update({
        "param_file": param_file,
        "param_file_snapshot": snapshot,
        "expt_hfe": expt_hfe,
        "opt_term_idx": opt_term_idx,
        "initial_params": initial_params,
        "autobar_path": autobar_path,
        "diff_step": diff_step,
    })

    print("\n=== Optimization Settings ===")
    print('bounds', bounds)
    print('initial_params', initial_params)
    print('diff_step', diff_step)

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

    # Run optimization. Leave x_scale at scipy's default so the trust
    # region scales with the Jacobian instead of treating sigma~3.8 and
    # epsilon~0.05 as numerically equivalent.
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
