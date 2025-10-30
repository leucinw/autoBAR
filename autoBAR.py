
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import time
import argparse
import subprocess
import numpy as np
import ruamel.yaml as yaml
from datetime import datetime
from utils.checkautobar import *
from utils.elescale import *
from scipy.optimize import least_squares

def setup():
  for phase in phases:
    xyzfile = phase_xyz[phase]
    keylines = phase_key[phase]
    with open(f"{phase}.key", 'w') as f:
      for line in keylines:
        if 'parameters' in line.lower():
          line = f'parameters     {prm}\n'
        f.write(line)
    gasminsh = 'gas-min.sh'
    if phase == 'gas':
      with open(gasminsh, 'w') as f:
        f.write(f'source {tinkerenv}\n')
        f.write(f'{phase_minimize[phase]} {xyzfile} -key gas.key 0.1 > gas-min.log \n')
        f.write(f'wait\nmv {phase_xyz[phase]}_2 {phase_xyz[phase]}\n')
      if not os.path.isfile("gas-min.log"):
        os.system(f"sh {gasminsh}")
    liquidminsh = 'liquid-min.sh'
    if phase == 'liquid':
      with open(liquidminsh, 'w') as f:
        f.write(f'source {tinkerenv}\n')
        f.write(f'{phase_minimize[phase]} {xyzfile} -key liquid.key 0.2 > liquid-min.log\n')
        f.write(f'wait\nmv {phase_xyz[phase]}_2 {phase_xyz[phase]}\n')
      if not os.path.isfile("liquid-min.log"):
        os.system(f"sh {liquidminsh}")
    os.system(f"rm -f {phase}/*.xyz {phase}/*.key")
    if not (ignoregas == 1 and phase == 'gas'):
      prmfiles = []
      for i in range(len(orderparams)):
        elb, vlb = orderparams[i]
        fname = "%s-e%s-v%s"%(phase, "%03d"%round(elb*100), "%03d"%round(vlb*100))
        with open(fname + ".key", 'w') as fw:
          for line in keylines:
            if 'parameters' in line.lower():
              if (elb*vlb > 1.0):
                assert elb == vlb, RED + f" Error: lambdas greater than 1 but not the same: {elb}, {vlb} " + ENDC
                idx = int(round(elb*100)/10) - 10
                perturbprm = f"{prm}_{idx:02d}"
                line = f'parameters     {homedir}/{phase}/{perturbprm}\n'
                if perturbprm not in prmfiles:
                  prmfiles.append(perturbprm)
              else:
                line = f'parameters     {homedir}/{phase}/{prm}\n'
                if prm not in prmfiles:
                  prmfiles.append(prm)
            fw.write(line)
          fw.write('\n')
          fw.write(f'ligand -1 {natom}\n')
          if (elb*vlb > 1.0):
            elb = 1.0
            vlb = 1.0
          
          if manualelescale:
            fw.write('ele-lambda 1.0\n')
          else:
            fw.write(f'ele-lambda {elb}\n')
          fw.write(f'vdw-lambda {vlb}\n')
          
          if manualelescale:
            scaledprms = scaledownele(phase_xyz['gas'], prm, elb)
            fw.write(f'\n# scale down electrostatic related parameters by {elb}\n')
            for s in scaledprms:
              fw.write(f'{s}\n')
            
        linkxyz = f"ln -sf {homedir}/{phase_xyz[phase]} {homedir}/{phase}/{fname}.xyz"
        movekey = f"mv {fname}.key ./{phase}"
        subprocess.run(linkxyz, shell=True)
        subprocess.run(movekey, shell=True)
      
      # make copy of prmfile
      for prmfile in prmfiles:
        copyprm = f"cp {prmfile}  ./{phase}"
        subprocess.run(copyprm, shell=True)
  print(GREEN + ' [GOOD] BAR simulation files generated!' + ENDC)
  return

def dynamic():
  liquidshs = []
  gasshs = []
  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    os.system(f"rm -rf {phasedir}/*.sh")
    os.chdir(phasedir)
    for i in range(len(orderparams)):
      elb, vlb = orderparams[i]
      fname = "%s-e%s-v%s"%(phase, "%03d"%round(elb*100), "%03d"%round(vlb*100))
      if (elb*vlb > 1.0) and copyarcforperturb:
        fname0 = "%s-e100-v100"%(phase)
        cmd = f"ln -sf {fname0}.arc {fname}.arc"
        os.system(cmd)
        print(GREEN + f" {cmd} " + ENDC)
      
      xyzfile = fname + ".xyz"
      keyfile = xyzfile.replace("xyz", "key")
      logfile = xyzfile.replace("xyz", "log")
      if os.path.isfile(logfile):
        if verbose > 0:
          print(YELLOW + f" [Warning] {logfile} exists in {phasedir} folder for {fname}" + ENDC)
      else:
        if phase == 'liquid':
          liquidsh = fname + '.sh'
          liquidshs.append(liquidsh)
          with open(liquidsh, 'w') as f:
            f.write(f'source {tinkerenv}\n')
            if liquidensemble == "NPT":
              f.write(f'{phase_dynamic[phase]} {xyzfile} -key {keyfile} {liquidtotalstep} {liquidtimestep} {liquidwriteout} 4 {liquidT} {liquidP} > {logfile}\n')
            if liquidensemble == "NVT":
              f.write(f'{phase_dynamic[phase]} {xyzfile} -key {keyfile} {liquidtotalstep} {liquidtimestep} {liquidwriteout} 2 {liquidT} > {logfile}\n')
        if phase == 'gas':
          gassh = fname + '.sh'
          gasshs.append(gassh)
          with open(gassh, 'w') as f:
            f.write(f'source {tinkerenv}\n')
            f.write(f"{phase_dynamic[phase]} {xyzfile} -key {keyfile} {gastotalstep} {gastimestep} {gaswriteout} 2 {gastemperature} > {logfile}\n")
  
  # submit jobs to clusters 
  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    os.chdir(phasedir)
    hoststr = os.getenv('HOSTNAME').split('.')[0]
    timestr = str(time.time()).replace('.', '')
    jobpooldir = "/home/liuchw/bin/JobPool/"
    scriptfile = os.path.join(jobpooldir, f"{hoststr}-{timestr}.sh")
    if (phase == 'gas') and (gasshs != []):
      if os.getlogin() == 'liuchw':
        shstr = f"python /home/liuchw/bin/TinkerGPU2022/submitTinker.py -x {' '.join(gasshs)} -t CPU -n 4 -p {phasedir}"
        print(GREEN + ' ' + shstr + ENDC)
        with open(scriptfile, 'w') as f:
          f.write(shstr)
      else:  
        shstr = f"python {submitexe} -x {' '.join(gasshs)} -t CPU -nodes {' '.join(nodes)} -n 4 -p {phasedir}"
        print(GREEN + ' ' + shstr + ENDC)
        os.system(shstr)
    if (phase == 'liquid') and (liquidshs != []):
      if os.getlogin() == 'liuchw':
        shstr = f"python /home/liuchw/bin/TinkerGPU2022/submitTinker.py -x {' '.join(liquidshs)} -t GPU -p {phasedir}"
        print(GREEN + ' ' + shstr + ENDC)
        with open(scriptfile, 'w') as f:
          f.write(shstr)
      else:  
        shstr = f"python {submitexe} -x {' '.join(liquidshs)} -t GPU -nodes {' '.join(nodes)} -p {phasedir}"
        print(GREEN + ' ' + shstr + ENDC)
        os.system(shstr)
  return

def bar():
  print(YELLOW + " Checking the completeness of the MD trajectories, please wait... " + ENDC)
  
  if copyarcforperturb:
    checkorderparams = orderparams[:-1]
  else:
    checkorderparams = orderparams
  if inputaction == "auto":
    proceed = False
    while not proceed: 
      proceed, phase_simsnapshot = checkdynamic(liquidtotalsnapshot, gastotalsnapshot, phases, checkorderparams, homedir, verbose)
      if proceed:
        break
      now = datetime.now().strftime("%b %d %Y %H:%M:%S")
      print(YELLOW + f" [{now}] Waiting for dynamic jobs to finish ..." + ENDC)
      time.sleep(checkingtime)
  else:
    proceed, phase_simsnapshot = checkdynamic(liquidtotalsnapshot, gastotalsnapshot, phases, checkorderparams, homedir, verbose)
    
  if proceed:
    liquidshs = []
    gasshs = []
    for phase in phases:
      phasedir = os.path.join(homedir, phase)
      os.chdir(phasedir)
      for i in range(len(orderparams)-1):
        elb0, vlb0 = orderparams[i]
        elb1, vlb1 = orderparams[i+1]
        if elb1 > 1.0:
          elb0, vlb0 = 1.0, 1.0 
        elb0 = "%03d"%round(elb0*100)
        elb1 = "%03d"%round(elb1*100)
        vlb0 = "%03d"%round(vlb0*100)
        vlb1 = "%03d"%round(vlb1*100)
        fname0 = "%s-e%s-v%s"%(phase, elb0, vlb0)
        fname1 = "%s-e%s-v%s"%(phase, elb1, vlb1)
        arcfile0 = fname0 + ".arc"
        arcfile1 = fname1 + ".arc"
        liquidbarname = f"bar_e{elb0}-v{vlb0}_e{elb1}-v{vlb1}.sh"
        gasbarname = f"bar_e{elb1}-v{vlb1}_e{elb0}-v{vlb0}.sh"
        barfiledir = os.path.join(homedir, phase)
        # put the FEP files in a separate folder
        if (int(elb1) > 100 and (int(vlb1) > 100)):
          idx = int(int(elb1)/10 - 10)
          fepdir = f"FEP_{idx:02d}"
          barfiledir = os.path.join(homedir, phase, fepdir)
          os.system(f"mkdir -p {homedir}/{phase}/{fepdir}") 
          linkarc = f"ln -sf {homedir}/{phase}/{arcfile0} {homedir}/{phase}/{fepdir}/{arcfile0}"
          os.system(linkarc)
          key0 = f"{phase}-e{elb0}-v{vlb0}.key"
          key1 = f"{phase}-e{elb1}-v{vlb1}.key"
          linkkey = f"ln -sf {homedir}/{phase}/{key0} {homedir}/{phase}/{fepdir}/{key0}"
          os.system(linkkey)
          linkkey = f"ln -sf {homedir}/{phase}/{key1} {homedir}/{phase}/{fepdir}/{key1}"
          os.system(linkkey)
          if copyarcforperturb:
            linkarc = f"ln -sf {homedir}/{phase}/{arcfile0} {homedir}/{phase}/{fepdir}/{arcfile1}"
            os.system(linkarc)
          else: 
            linkarc = f"ln -sf {homedir}/{phase}/{arcfile1} {homedir}/{phase}/{fepdir}/{arcfile1}"
            os.system(linkarc)
        
        if phase == 'liquid':
          printname = fname0
          outfile = fname0 + ".out"
          barfile = fname0 + ".bar"
          enefile = fname0 + ".ene"
          startsnapshot = int(liquidtotalsnapshot/5.0) + 1
          
          if (not os.path.isfile(barfiledir + "/" + barfile)):
            with open(liquidbarname, 'w') as f:
              f.write(f"source {tinkerenv}\n")
              f.write(f"{liquidbarexe} 1 {arcfile0} {liquidT} {arcfile1} {liquidT} N > {outfile} && \n")
              f.write(f"{liquidbarexe} 2 {barfile} {startsnapshot} {liquidtotalsnapshot} 1 {startsnapshot} {liquidtotalsnapshot} 1 > {enefile} \n")
            if [liquidbarname, barfiledir] not in liquidshs:
              liquidshs.append([liquidbarname, barfiledir])
            if (int(elb1) > 100):
              os.system(f"mv {liquidbarname} {homedir}/{phase}/{fepdir}")
          else:
            if verbose > 0:
              print(YELLOW + f" [Warning] {barfile} exists in {barfiledir} folder! " + ENDC)
        
        if phase == 'gas':
          printname = fname1
          outfile = fname1 + ".out"
          barfile = fname1 + ".bar"
          enefile = fname1 + ".ene"
          startsnapshot = int(gastotalsnapshot/5.0) + 1
          
          if gastotaltime != 0:
            if (not os.path.isfile(barfiledir + '/' + barfile)):
              with open(gasbarname, 'w') as f:
                f.write(f"source {tinkerenv}\n")
                f.write(f"{gasbarexe} 1 {arcfile1} {gastemperature} {arcfile0} {gastemperature} N > {outfile} && \n")
                f.write(f"{gasbarexe} 2 {barfile} {startsnapshot} {gastotalsnapshot} 1 {startsnapshot} {gastotalsnapshot} 1 > {enefile} \n")
              if [gasbarname, barfiledir] not in gasshs:
                gasshs.append([gasbarname, barfiledir])
              if (int(elb1) > 100):
                os.system(f"mv {gasbarname} {homedir}/{phase}/{fepdir}")
            else:
              if verbose > 0:
                print(YELLOW + f" [Warning] {barfile} exists in {barfiledir} folder!" + ENDC)
    # submit jobs to clusters 
    for phase in phases:
      phasedir = os.path.join(homedir, phase)
      os.chdir(phasedir)
      if (phase == 'gas') and (gasshs != []):
        if os.getlogin() == 'liuchw':
          os.system(f"python /home/liuchw/bin/TinkerGPU2022/submitTinker.py -x {' '.join([x[0] for x in gasshs])} -t CPU -n 4 -p {' '.join([x[1] for x in gasshs])}")
        else:
          os.system(f"python {submitexe} -x {' '.join([x[0] for x in gasshs])} -t CPU -n 4 -p {' '.join([x[1] for x in gasshs])} -nodes {' '.join(nodes)}")
      if (phase == 'liquid') and (liquidshs != []):
        if os.getlogin() == 'liuchw':
          os.system(f"python /home/liuchw/bin/TinkerGPU2022/submitTinker.py -x {' '.join([x[0] for x in liquidshs])} -t GPU -p {' '.join([x[1] for x in liquidshs])}")
        else:
          os.system(f"python {submitexe} -x {' '.join([x[0] for x in liquidshs])} -t GPU -p {' '.join([x[1] for x in liquidshs])} -nodes {' '.join(nodes)}")
  return

def result():
  print(YELLOW + " Checking the completeness of the BAR analysis ..." + ENDC)
  proceed, gasperturbsteps, gasenes, liquidperturbsteps, liquidenes, fep_gasperturbsteps, fep_gasenes, fep_liquidperturbsteps, fep_liquidenes = checkbar(phases, orderparams, homedir, ignoregas, verbose)
  
  if proceed:
    FEgas = []
    Errgas = []
    FEliquid = []
    Errliquid = []
    fep_FEgas = []
    fep_Errgas = []
    fep_FEliquid = []
    fep_Errliquid = []
    
    for gasene in gasenes:
      find = False
      for line in open(gasene).readlines():
        if "Free Energy via BAR Iteration" in line:
          find = True
          FEgas.append(float(line.split()[-4]))
          Errgas.append(float(line.split()[-2]))
        if (not find) and ("Free Energy via BAR Bootstrap" in line):
          find = True
          FEgas.append(float(line.split()[-4]))
          Errgas.append(float(line.split()[-2]))
    
    for gasene in fep_gasenes:
      find = False
      for line in open(gasene).readlines():
        if "Free Energy via BAR Iteration" in line:
          find = True
          fep_FEgas.append(float(line.split()[-4]))
          fep_Errgas.append(float(line.split()[-2]))
        if (not find) and ("Free Energy via BAR Bootstrap" in line):
          find = True
          fep_FEgas.append(float(line.split()[-4]))
          fep_Errgas.append(float(line.split()[-2]))
    
    for liquidene in liquidenes:
      find = False
      for line in open(liquidene).readlines():
        if "Free Energy via BAR Iteration" in line:
          find = True
          FEliquid.append(float(line.split()[-4]))
          Errliquid.append(float(line.split()[-2]))
        if (not find) and ("Free Energy via BAR Bootstrap" in line):
          find = True
          FEliquid.append(float(line.split()[-4]))
          Errliquid.append(float(line.split()[-2]))
    
    for liquidene in fep_liquidenes:
      find = False
      for line in open(liquidene).readlines():
        if "Free Energy via BAR Iteration" in line:
          find = True
          fep_FEliquid.append(float(line.split()[-4]))
          fep_Errliquid.append(float(line.split()[-2]))
        if (not find) and ("Free Energy via BAR Bootstrap" in line):
          find = True
          fep_FEliquid.append(float(line.split()[-4]))
          fep_Errliquid.append(float(line.split()[-2]))
    
    FEall = FEgas + FEliquid
    Errall = Errgas + Errliquid
    totFE = np.array(FEall).sum()
    totErr = np.sqrt(np.square(np.array(Errall)).sum())
    
    fo = open(os.path.join(homedir, "result.txt"), "w")
    print(GREEN + "%20s%20s%30s%20s"%("StateA", "StateB", "FreeEnergy(kcal/mol)", "Error(kcal/mol)") + ENDC)
    fo.write("%20s%20s%30s%20s\n"%("StateA", "StateB", "FreeEnergy(kcal/mol)", "Error(kcal/mol)"))
    for i in range(len(FEgas)-1,-1,-1):
      state0, state1 = gasperturbsteps[i]
      print("%20s%20s%25.4f%20.4f"%(state0, state1, FEgas[i], Errgas[i]))
      fo.write("%20s%20s%25.4f%20.4f\n"%(state0, state1, FEgas[i], Errgas[i]))
    for i in range(len(FEliquid)):
      state0, state1 = liquidperturbsteps[i]
      print("%20s%20s%25.4f%20.4f"%(state0, state1, FEliquid[i], Errliquid[i]))
      fo.write("%20s%20s%25.4f%20.4f\n"%(state0, state1, FEliquid[i], Errliquid[i]))
    print(GREEN + "%40s%25.4f%20.4f"%("SUM OF THE TOTAL FREE ENERGY (FE0)", totFE, totErr) + ENDC)
    fo.write("%40s%25.4f%20.4f\n"%("SUM OF THE TOTAL FREE ENERGY (FE0)", totFE, totErr))
    
    print(YELLOW + "    %20s%20s%27s%20s"%("GAS", "LIQUID", "GAS+LIQUID", "GAS+LIQUID+FE0") + ENDC)
    fo.write("    %20s%20s%27s%20s\n"%("GAS", "LIQUID", "GAS+LIQUID", "GAS+LIQUID+FE0"))
    for i in range(len(fep_FEliquid)): 
      feg = 0.0
      if ignoregas == 0:
        feg = fep_FEgas[i]
      fel = fep_FEliquid[i]
      print(f"    FEP_{i+1:03d}{feg:14.4f}{fel:19.4f}{feg+fel:25.4f}{totFE+feg+fel:20.4f}")
      fo.write(f"    FEP_{i+1:03d}{feg:14.4f}{fel:19.4f}{feg+fel:25.4f}{totFE+feg+fel:20.4f}\n")
    
    fo.close()
  return

def opt():
  print(GREEN + " Optimizing parameters ..." + ENDC)
  p_init = []
  terms = []
  types = []
  indices = []
  deltas = []
  for t in tuning_parms:
    s = t.split("_")
    terms.append(s[0])
    types.append(s[1])
    indices.append(int(s[2]))
    deltas.append(float(s[3]))
  
  lines = open(prm).readlines()
  for line in lines:
    ss = line.split()
    if len(ss) > 3:
      for i in range(len(terms)):
        if (ss[0] == terms[i]) and (ss[1] == types[i]):
          p_init.append(float(ss[indices[i]]))
  p_init = np.array(p_init)

  def write_ff(p):
    with open(prm + "_10", 'w') as f:
      k = 0
      for line in lines:
        tuning = False
        ss = line.split()
        if len(ss) > 3:
          for i in range(len(terms)):
            if (ss[0] == terms[i]) and (ss[1] == types[i]):
              ss[indices[i]] = f"{p[k]:10.5f}"
              k += 1
              tuning = True
        if tuning:
          line = '  '.join(ss) + "\n"
        f.write(line)
    return
  
  def sp_fe(p):
    restraint_factor = 1.0
    os.system(f"rm -f result.txt {prm}_10")
    if not np.array_equal(p, p_init):
      write_ff(p)
      os.system("rm -f result.txt */*200* liquid/liquid-e100-v100.{bar,ene}")
    
    cmdstr = f"python {os.path.join(rootdir, 'autoBAR.py')} auto"
    os.system(cmdstr)
    
    while True:
      if not os.path.isfile('result.txt'):
        time.sleep(30.0)
      else:
        calc_fe = float(open('result.txt').readlines()[-2].split()[-2])
        break 
    return np.array([calc_fe - expt_fe] + list((p-p_init)*restraint_factor))
 
  upper_bound = []
  lower_bound = []
  for i in range(len(p_init)):
    upper_bound.append(p_init[i] + deltas[i])
    lower_bound.append(p_init[i] - deltas[i])
  ret = least_squares(sp_fe, p_init, loss='soft_l1', bounds = (lower_bound, upper_bound), verbose=2, diff_step=1e-5, ftol=1e-3, gtol=1e-3, xtol=1e-3)
  return

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('act', help = "Actions to take.", choices = ['setup', 'dynamic', 'bar', 'result', 'auto', 'opt'], type = str.lower) 
  
  # colors for screen output
  global RED, ENDC, GREEN, YELLOW
  RED = '\033[91m'
  ENDC = '\033[0m'
  GREEN = '\033[92m'
  YELLOW = '\33[93m'

  # global settings 
  global inputaction
  inputaction = vars(parser.parse_args())['act']
  global rootdir, homedir
  rootdir = os.path.join(os.path.split(__file__)[0])
  homedir = os.getcwd()
  if not os.path.isfile(os.path.join(homedir, "settings.yaml")):
    sys.exit(RED + "Please provide 'settings.yaml' file; " + ENDC + GREEN + f"An example is here for you {os.path.join(rootdir, 'dat', 'settings.yaml')}" + ENDC)
  else:
    with open('settings.yaml') as f:
      FEsimsettings = yaml.load(f,Loader=yaml.Loader)
  
  global verbose
  verbose = 1
  if 'verbose' in FEsimsettings.keys():
    verbose = int(FEsimsettings['verbose'])
  
  global prm, checkingtime
  lig = FEsimsettings['gas_xyz'] 
  box = FEsimsettings['box_xyz'] 
  prm = FEsimsettings["parameters"]
  checkingtime = FEsimsettings["checking_time"]
  global natom
  natom = int(open(lig).readlines()[0].split()[0])
  global nodes
  try:
    nodes = FEsimsettings["node_list"]
  except:
    if os.getlogin() != 'liuchw':
      sys.exit(RED + "node_list must be provided" + ENDC)
    else:
      pass
 
  if inputaction == 'opt':
    global expt_fe, tuning_parms 
    if "expt_fe" not in FEsimsettings.keys():
      sys.exit(RED + " I could not find expt_fe in settings.yaml!" + ENDC)
    if "tuning_parms" not in FEsimsettings.keys():
      sys.exit(RED + " I could not find tuning_parms in settings.yaml!" + ENDC)
    expt_fe = FEsimsettings["expt_fe"]
    tuning_parms = FEsimsettings["tuning_parms"]
    
  # special list for liuchw
  if os.getlogin() == 'liuchw':
    nodes = []
    lines = open('/home/liuchw/bin/TinkerGPU2022/nodes.dat').readlines()
    for line in lines:
      if ("GPU " in line) and ("#" not in line[0]):
        d = line.split()
        if d[1] not in nodes:
          nodes.append(d[1])
  
  if nodes == None:
    sys.exit(RED + "[Error] node_list must be provided" + ENDC) 
  global submitexe
  submitexe = os.path.join(rootdir, "utils", "submitTinker.py")
  global tinkerenv 
  tinkerenv = os.path.join(rootdir, "dat", "tinker.env")
  global orderparams
  orderparams = []
  lambdawindow = FEsimsettings["lambda_window"].upper()
  if lambdawindow == "COURSER":
    orderprmfile = os.path.join(rootdir, "dat", "orderparams_courser")
  else:
    orderprmfile = os.path.join(rootdir, "dat", "orderparams_default")
  # user specified order parameter file 
  try:
    orderprmfile = FEsimsettings['lambda_window_file']
    print(YELLOW + f" [Warning] You are responsible for your {orderprmfile}" + ENDC)
  except:
    print(GREEN + f" [GOOD] You are using lambda window settings from {orderprmfile}" + ENDC) 

  global manualelescale
  try:
    manualelescale = FEsimsettings['manual_ele_scale']
  except:
    manualelescale = False 
  
  for line in open(orderprmfile).readlines():
    if "#" not in line:
      d = line.split()
      orderparams.append([float(d[0]), float(d[1])])
  
  global copyarcforperturb
  copyarcforperturb = False
  if "copy_arc_for_perturb" in FEsimsettings.keys():
    copyarcforperturb = bool(FEsimsettings["copy_arc_for_perturb"])
    print(YELLOW + f" [Warning] Copy ARC file from e100-v100: {bool(copyarcforperturb)}" + ENDC) 
   
  for i in range(1,100,1): # allow 99 files at a time to do FEP
    if os.path.isfile(os.path.join(homedir, f"{prm}_{i:02d}")):
      fakelambda = round(1.0 + i*0.1, 1)
      orderparams.append([fakelambda, fakelambda])
      print(GREEN + f" [GOOD] You are doing one-step FEP for parameter file {prm}_{i:02d}" + ENDC) 


  # liquid phase specific settings 
  global phases, liquidkeylines
  phases = ['liquid']
  try:
    liquidkey = FEsimsettings['liquid_key']
    print(YELLOW + f"[Warning] You are responsible for your {liquidkey}" + ENDC)
  except:
    liquidkey = os.path.join(rootdir, "dat", "liquid.key")
  liquidkeylines = open(liquidkey).readlines()
  # tinker.env
  global liquidmdexe, liquidbarexe
  liquidmdexe = '$DYNAMIC9' 
  liquidminexe = '$MINIMIZE9' 
  liquidbarexe = '$BAR9' 
  global phase_xyz, phase_key, phase_dynamic, phase_minimize
  phase_xyz = {'liquid':box}
  phase_key = {'liquid':liquidkeylines}
  global liquidtotaltime, liquidtimestep, liquidwriteout, liquidtotalstep, liquidtotalsnapshot
  global liquidT, liquidP, liquidensemble
  liquidtotaltime = FEsimsettings["liquid_md_total_time"] 
  liquidtimestep = FEsimsettings["liquid_md_time_step"] 
  liquidwriteout = FEsimsettings["liquid_md_write_freq"] 
  liquidT = FEsimsettings["liquid_md_temperature"] 
  liquidP = FEsimsettings["liquid_md_pressure"]
  liquidensemble = FEsimsettings["liquid_md_ensemble"].upper()
  liquidtotalstep = int((1000000.0*liquidtotaltime)/liquidtimestep)
  liquidtotalsnapshot = int(1000*liquidtotaltime/liquidwriteout)
  if not os.path.isdir('liquid'):
    os.system("mkdir liquid")

  # gas phase specific settings
  global ignoregas
  global gaskeylines
  try:
    gaskey = FEsimsettings['gas_key']
    print(YELLOW + f"[Warning] You are responsible for your {gaskey}" + ENDC)
  except:
    gaskey = os.path.join(rootdir, "dat", "gas.key")
  
  gaskeylines = open(gaskey).readlines()
  # tinker.env
  global gasmdexe, gasbarexe
  gasmdexe = '$DYNAMIC8' 
  gasminexe = '$MINIMIZE8' 
  gasbarexe = '$BAR8' 
  global gastotaltime, gastimestep, gaswriteout, gastemperature, gastotalstep, gastotalsnapshot
  gastotaltime = FEsimsettings["gas_md_total_time"] 
  gastimestep = FEsimsettings["gas_md_time_step"] 
  gastemperature = FEsimsettings["gas_md_temperature"] 
  gaswriteout = FEsimsettings["gas_md_write_freq"]
  gastotalstep = int((1000000.0*gastotaltime)/gastimestep)
  gastotalsnapshot = int(1000*gastotaltime/gaswriteout)
  phase_xyz = {'liquid':box, 'gas':lig}
  phase_key = {'liquid':liquidkeylines, 'gas':gaskeylines}
  phase_dynamic = {'liquid':liquidmdexe, 'gas':gasmdexe}
  phase_minimize = {'liquid':liquidminexe, 'gas':gasminexe}
  
  # check xyz files
  lines = open(box).readlines()
  natomliquid = int(lines[0].split()[0])
  if natomliquid != len(lines)-2:
    print(YELLOW + f"[Warning] No box info in {box}. Can be in liquid.key instead." + ENDC)
  else:
    [a,b,c] = lines[1].split()[0:3]
    if min([float(a), float(b), float(c)]) < 30.0:
      sys.exit(RED + f"[Error] Please provide a bigger box (>30*30*30)" + ENDC)
  lines = open(lig).readlines()
  natomgas = int(lines[0].split()[0])
  if natomgas < 5:
    gastotaltime = 0.0
    print(YELLOW + f" [Warning] I set the simulation time to 0 since it only contains {natomgas} atoms" + ENDC)
  if gastotaltime == 0.0:
    ignoregas = 1
  else:
    ignoregas = 0
  if (ignoregas == 0):
    phases.insert(0,'gas')
    if not os.path.isdir('gas'):
      os.system("mkdir gas")
  
  actions = {'setup':setup, 'dynamic':dynamic, 'bar':bar, 'result':result, 'opt':opt}
  if inputaction in actions.keys():
    actions[inputaction]()
  if inputaction == 'auto': 
    actions['setup']()
    actions['dynamic']()
    
    if copyarcforperturb:
      checkorderparams = orderparams[:-1]
    else:
      checkorderparams = orderparams
    
    # check if dynamic is complete
    dynamic_good = False
    while not dynamic_good: 
      time.sleep(30.0)
      dynamic_good, _ = checkdynamic(liquidtotalsnapshot, gastotalsnapshot, phases, checkorderparams, homedir, verbose)
    
    actions['bar']()
    # check if bar is complete  
    bar_good = False
    while not bar_good:
      time.sleep(30.0)
      bar_good, _, _, _, _, _, _, _, _ = checkbar(phases, orderparams, homedir, ignoregas, verbose)
    
    actions['result']()
