#!/usr/bin/env python3
"""
Used to do parameter optimization targeting expt. HFE
Necessary settings are in the main settings.yaml file

Example usage: python parmOPT.py

- Chengwen Liu
- 2025/11

"""

import os
import time
import shutil
import numpy as np
import ruamel.yaml as yaml
from scipy.optimize import least_squares

def model_func(params):
    """
    Compute residuals (difference between predicted and experimental values).
    """
    if os.path.isfile('result.txt'):
      os.system('rm -f result.txt')

    if np.all(params == initial_params):
      write_prm(params, param_file)
    else:
      os.system(f'rm -f */{param_file}_01')
      os.system(f'rm -f {param_file}_01')
      os.system('rm -f */*e110*')
      os.system('rm -rf */FEP_01')
      write_prm(params, param_file+"_01")
    
    os.system(f'python {autobar_path} auto')
    
    # wait for result.txt file
    while True:
      if os.path.isfile('result.txt'):
        break
      else:
        time.sleep(60.0)
    
    # read result.txt file
    lines = open('result.txt').readlines()
    fe0 = 0.0
    fe1 = 0.0
    for line in lines:
      if 'SUM OF ' in line:
        fe0 = float(line.split()[-2])
      if 'FEP_001' in line:
        fe1 = float(line.split()[-1])

    print(f'Current params: {params}')
    if np.all(params == initial_params):
      print(f'Current HFE: {fe0}')
      return np.array([fe0 - expt_hfe])
    else:
      print(f'Current HFE: {fe1}')
      return np.array([fe1 - expt_hfe])

def write_prm(params, fname):
    lines = open(param_file).readlines()
    if not (param_file == fname):
      shutil.copy(param_file, fname)
    with open(fname, 'a') as f:
      line = opt_term_idx.replace('-', '   ') + '  ' + '  '.join([str(p) for p in params]) + '\n'
      f.write(line)  
    return

def jacobian_fd(params):
    """
    Compute Jacobian matrix numerically using central finite difference.
    """
    n_params = len(params)
    residuals = np.atleast_1d(model_func(params))
    m = len(residuals)

    J = np.zeros((m, n_params))
    step = diff_step * np.ones(n_params)
    
    if os.path.isfile('result.txt'):
      os.system('rm -f result.txt')
    
    perturb_idx = 1
    if not np.all(params == initial_params):
      perturb_idx = 2 
    
    for j in range(n_params):
      # deep copy the numpy array
      params_plus = params.copy()
      params_minus = params.copy()
      dp = np.zeros_like(params)
      dp[j] = step[j]
      
      # plus finite difference   
      lambda_str = f"{100 + perturb_idx*10}"
      param_file_p = param_file + f'_{perturb_idx:02d}' 
      os.system(f'rm -f */{param_file_p}')
      os.system(f'rm -f {param_file_p}')
      os.system(f'rm -rf */FEP_{perturb_idx:02d}')
      os.system(f'rm -f */*e{lambda_str}*')
      params_plus += dp
      write_prm(params_plus, param_file_p)
      perturb_idx += 1
     
      # minus finite difference
      param_file_p = param_file + f'_{perturb_idx:02d}' 
      os.system(f'rm -f */{param_file_p}')
      os.system(f'rm -f {param_file_p}')
      os.system(f'rm -rf */FEP_{perturb_idx:02d}')
      os.system(f'rm -f */*e{lambda_str}*')
      params_minus -=  dp
      write_prm(params_minus, param_file_p)
      perturb_idx += 1
      
    os.system(f'python {autobar_path} auto')
    
    # wait for result.txt file
    while True:
      if os.path.isfile('result.txt'):
        break
      else:
        time.sleep(60.0)
    
    # read result.txt file
    lines = open('result.txt').readlines()
    
    # reset starting idx
   
    feps = []
    for line in lines:
      if 'FEP_' in line:
        fep = float(line.split()[-1])
        feps.append(fep)

    for j in range(n_params):
      if np.all(params == initial_params):
        r_plus  = fep[j*2] 
        r_minus = fep[j*2 + 1] 
      else: 
        r_plus  = fep[j*2 + 1] 
        r_minus = fep[j*2 + 2] 
      J[:, j] = (r_plus - r_minus) / (2 * diff_step)

    return J

def main():
    with open('settings.yaml') as f:
      FEsimsettings = yaml.load(f,Loader=yaml.Loader)
    
    global param_file, expt_hfe, opt_params
    param_file = FEsimsettings["parameters"]
    expt_hfe = float(FEsimsettings["expt_hfe"])

    opt_params = FEsimsettings["opt_params"]
    params_range = FEsimsettings["params_range"]
   
    global autobar_path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    autobar_path = os.path.join(script_dir.split('utils')[0], 'autoBAR.py')

    # opt_params format: "vdwpair-401-402 3.8 0.05"
    s = opt_params.split()
    global opt_term_idx
    opt_term_idx = s[0]
    global initial_params
    initial_params = np.array([float(x) for x in s[1:]])
    
    # form the upper and lower bound
    lb = np.zeros(len(initial_params)) 
    ub = np.zeros(len(initial_params)) 

    for i in range(len(initial_params)):
      lb[i] = initial_params[i] - float(params_range.split()[i])
      ub[i] = initial_params[i] + float(params_range.split()[i])
   
    bounds = (lb, ub)
    
    # x*diff_step
    global diff_step
    diff_step = 0.0001
    
    # xs = x/x_scale
    x_scale = np.ones(len(initial_params))
    
    print("\n=== Optimization Settings ===")
    print('bounds', bounds) 
    print('initial_params', initial_params) 
    print('diff_step', diff_step) 
    print('x_scale', x_scale) 
    
    # clean up the FEP_01...09 files when starting
    os.system(f'rm -f */{param_file}_0?')
    os.system(f'rm -f {param_file}_0?')
    os.system('rm -f */*e1{10..90}*')
    os.system('rm -rf */FEP_{01..09}')
    
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
        bounds = bounds,
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