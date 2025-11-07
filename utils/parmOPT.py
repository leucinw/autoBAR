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
      os.system('rm -f {param_file}_01')
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
    os.system(f"cp {param_file} {fname}")
    with open(fname, 'a') as f:
      line = opt_term_idx.replace('-', '   ') + '  ' + '  '.join([str(p) for p in params]) + '\n'
      f.write(line)  
    return

# callback argument only valid for scipy > 1.16
def print_step(params):
    print(f"Current parameters: {params}")
          
def main():
    with open('settings.yaml') as f:
      FEsimsettings = yaml.load(f,Loader=yaml.Loader)
    
    global param_file, expt_hfe, opt_params
    param_file = FEsimsettings["parameters"]
    expt_hfe = float(FEsimsettings["expt_hfe"])

    opt_params = FEsimsettings["opt_params"]
   
    global autobar_path
    autobar_path = "/home/liuchw/Documents/Github.leucinw/autoBAR/autoBAR.py"

    # opt_params format: "vdwpair-401-402 3.8 0.05"
    s = opt_params.split()
    global opt_term_idx
    opt_term_idx = s[0]
    global initial_params
    initial_params = np.array([float(x) for x in s[1:]])
    
    # step size used to do finite difference
    diff_step = np.array([0.0005, 0.0001])

    # Run optimization
    result = least_squares(
        fun=model_func,
        x0=initial_params,
        jac='2-point',
        diff_step=diff_step,
        method='trf', 
        verbose=2,
        callback=print_step,
    )

    print("\n=== Optimization Results ===")
    print("Success:", result.success)
    print("Message:", result.message)
    print("Optimal parameters:", result.x)
    print("Cost (sum of squared residuals):", 2 * result.cost)

if __name__ == "__main__":
    main()