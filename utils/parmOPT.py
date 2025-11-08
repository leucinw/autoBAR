#!/usr/bin/env python3
"""
Used to do parameter optimization targeting expt. HFE
Necessary settings are in the main settings.yaml file

Example usage: python parmOPT.py

- Chengwen Liu
- 2025/11

"""

import os
import sys
import time
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
      if 'FEP_01' in line:
        fe1 = float(line.split()[-1])

    if np.all(params == initial_params):
      return np.array([fe0 - expt_hfe])
    else:
      return np.array([fe1 - expt_hfe])

def write_prm(params, fname):
    lines = open(param_file).readlines()
    if not np.all(params == initial_params):
      os.system(f"cp {param_file} {fname}")
    with open(fname, 'a') as f:
      line = opt_term_idx.replace('-', '   ') + '  ' + '  '.join([str(p) for p in params]) + '\n'
      f.write(line)  
    return

# callback function, require scipy>1.16!!
def print_step(params):
    print(f"Current parameters: {params}")
          
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
    
    # diff_step
    diff_step = 0.0001
    # xs = x/x_scale
    x_scale = np.array([10.0, 1.0])
    
    print("\n=== Optimization Settings ===")
    print('bounds', bounds) 
    print('initial_params', initial_params) 
    print('diff_step', diff_step) 
    print('x_scale', x_scale) 
    
    # clean up the FEP_01 files when starting
    os.system(f'rm -f */{param_file}_01')
    os.system(f'rm -f {param_file}_01')
    os.system('rm -f */*e110*')
    os.system('rm -rf */FEP_01')
    
    # Run optimization
    result = least_squares(
        fun=model_func,
        x0=initial_params,
        jac='2-point',
        diff_step=diff_step,
        x_scale=x_scale,
        loss='soft_l1',
        method='trf', 
        verbose=2,
        callback=print_step,
        bounds = bounds,
        ftol=0.001,
        gtol=0.001,
        xtol=0.001,
    )

    print("\n=== Optimization Results ===")
    print("Success:", result.success)
    print("Message:", result.message)
    print("Optimal parameters:", result.x)
    print("Cost (sum of squared residuals):", 2 * result.cost)

if __name__ == "__main__":
    main()