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
import time
from pathlib import Path

import numpy as np
import ruamel.yaml as yaml
from scipy.optimize import least_squares

# Module-level config dict, replacing scattered global declarations.
# Populated once in main() and read by model_func / write_prm / jacobian_fd.
_config = {}

# Maximum time (seconds) to wait for result.txt before raising an error
_MAX_WAIT_SECONDS = 86400  # 24 hours


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


def _wait_for_result(poll_interval=60.0, timeout=_MAX_WAIT_SECONDS):
    """Block until ``result.txt`` appears, raising on timeout."""
    elapsed = 0.0
    while not os.path.isfile('result.txt'):
        if elapsed >= timeout:
            raise TimeoutError(
                f"result.txt not found after waiting {timeout} seconds"
            )
        time.sleep(poll_interval)
        elapsed += poll_interval


def model_func(params):
    """
    Compute residuals (difference between predicted and experimental values).
    """
    param_file = _config["param_file"]
    initial_params = _config["initial_params"]
    autobar_path = _config["autobar_path"]
    expt_hfe = _config["expt_hfe"]

    Path('result.txt').unlink(missing_ok=True)

    if np.array_equal(params, initial_params):
        write_prm(params, param_file)
    else:
        _remove_matching(f'*/{param_file}_01')
        Path(f'{param_file}_01').unlink(missing_ok=True)
        _remove_matching('*/*e110*')
        _remove_matching('*/FEP_01')
        write_prm(params, param_file + "_01")

    subprocess.run(['python', autobar_path, 'auto'], check=False)

    _wait_for_result()

    # read result.txt file
    with open('result.txt') as f:
        lines = f.readlines()

    fe0 = 0.0
    fe1 = 0.0
    for line in lines:
        if 'SUM OF ' in line:
            fe0 = float(line.split()[-2])
        if 'FEP_001' in line:
            fe1 = float(line.split()[-1])

    print(f'Current params: {params}')
    if np.array_equal(params, initial_params):
        print(f'Current HFE: {fe0}')
        return np.atleast_1d([fe0 - expt_hfe])
    else:
        print(f'Current HFE: {fe1}')
        return np.atleast_1d([fe1 - expt_hfe])


def write_prm(params, fname):
    param_file = _config["param_file"]
    opt_term_idx = _config["opt_term_idx"]

    if not os.path.isfile(param_file):
        raise FileNotFoundError(f"Parameter file not found: {param_file}")
    if param_file != fname:
        shutil.copy(param_file, fname)
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
    autobar_path = _config["autobar_path"]
    diff_step = _config["diff_step"]

    n_params = len(params)
    params = np.atleast_1d(params)

    # only ONE HFE data
    J = np.zeros((1, n_params))
    step = diff_step * np.ones(n_params)

    Path('result.txt').unlink(missing_ok=True)

    perturb_idx = 1
    if not np.array_equal(params, initial_params):
        perturb_idx = 2

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
        perturb_idx += 1

        # minus finite difference
        param_file_p = param_file + f'_{perturb_idx:02d}'
        _remove_matching(f'*/{param_file_p}')
        Path(param_file_p).unlink(missing_ok=True)
        _remove_matching(f'*/FEP_{perturb_idx:02d}')
        _remove_matching(f'*/*e{lambda_str}*')
        params_minus -= dp
        write_prm(params_minus, param_file_p)
        perturb_idx += 1

    subprocess.run(['python', autobar_path, 'auto'], check=False)

    _wait_for_result()

    # read result.txt file
    with open('result.txt') as f:
        lines = f.readlines()

    feps = []
    for line in lines:
        if 'FEP_' in line:
            fep = float(line.split()[-1])
            feps.append(fep)

    for j in range(n_params):
        if np.array_equal(params, initial_params):
            r_plus = feps[j * 2]
            r_minus = feps[j * 2 + 1]
        else:
            r_plus = feps[j * 2 + 1]
            r_minus = feps[j * 2 + 2]
        J[:, j] = (r_plus - r_minus) / (2 * diff_step)

    return J


def main():
    with open('settings.yaml') as f:
        settings = yaml.load(f, Loader=yaml.Loader)

    param_file = settings["parameters"]
    expt_hfe = float(settings["expt_hfe"])
    opt_params = settings["opt_params"]
    params_range = settings["params_range"]

    script_dir = os.path.dirname(os.path.abspath(__file__))
    autobar_path = os.path.join(script_dir.split('utils')[0], 'autoBAR.py')

    # opt_params format: "vdwpair-401-402 3.8 0.05"
    s = opt_params.split()
    opt_term_idx = s[0]
    initial_params = np.atleast_1d([float(x) for x in s[1:]])

    # x*diff_step
    diff_step = 0.0001

    # Populate the module-level config dict (replaces scattered globals)
    _config.update({
        "param_file": param_file,
        "expt_hfe": expt_hfe,
        "opt_term_idx": opt_term_idx,
        "initial_params": initial_params,
        "autobar_path": autobar_path,
        "diff_step": diff_step,
    })

    # form the upper and lower bound
    lb = np.zeros(len(initial_params))
    ub = np.zeros(len(initial_params))
    range_values = params_range.split()

    for i in range(len(initial_params)):
        lb[i] = initial_params[i] - float(range_values[i])
        ub[i] = initial_params[i] + float(range_values[i])

    bounds = (lb, ub)

    # xs = x/x_scale
    x_scale = np.ones(len(initial_params))

    print("\n=== Optimization Settings ===")
    print('bounds', bounds)
    print('initial_params', initial_params)
    print('diff_step', diff_step)
    print('x_scale', x_scale)

    # clean up the FEP_01...09 files when starting
    _remove_matching(f'*/{param_file}_0?')
    _remove_matching(f'{param_file}_0?')
    for i in range(1, 10):
        _remove_matching(f'*/*e1{i}0*')
        _remove_matching(f'*/FEP_0{i}')

    # Run optimization
    result = least_squares(
        fun=model_func,
        x0=initial_params,
        jac=jacobian_fd,
        diff_step=diff_step,
        x_scale=x_scale,
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
