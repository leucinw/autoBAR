
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import time
import yaml
import argparse
import subprocess
import numpy as np
from datetime import datetime
from utils.checkautobar import *

# color
RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\33[93m'

def setup():
  for phase in phases:
    for i in range(len(orderparams)):
      elb, vlb = orderparams[i]
      fname = "%s-e%s-v%s"%(phase, "%03d"%int(elb*100), "%03d"%int(vlb*100))
      with open(fname + ".key", 'w') as fw:
        for line in phase_key[phase]:
          if 'parameters' in line.lower():
            line = f'parameters     ../{prm}\n'
          fw.write(line)
        fw.write('\n')
        fw.write(f'ligand -1 {natom}\n')
        fw.write(f'ele-lambda {elb}\n')
        fw.write(f'vdw-lambda {vlb}\n')
      linkxyz = f"ln -sf {homedir}/{phase_xyz[phase]} {homedir}/{phase}/{fname}.xyz"
      movekey = f"mv {fname}.key ./{phase}"
      subprocess.run(linkxyz, shell=True)
      subprocess.run(movekey, shell=True)
  print(GREEN + ' [GOOD] BAR simulation files generated!' + ENDC)
  return

def dynamic():
  phase_dynamic = {'liquid':liquidmdexe, 'gas':gasmdexe}
  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    os.chdir(phasedir)
    for i in range(len(orderparams)):
      elb, vlb = orderparams[i]
      fname = "%s-e%s-v%s"%(phase, "%03d"%int(elb*100), "%03d"%int(vlb*100))
      xyzfile = fname + ".xyz"
      keyfile = xyzfile.replace("xyz", "key")
      logfile = xyzfile.replace("xyz", "log")
      if (not os.path.isfile(logfile)):
        if phase == 'liquid':
          if liquidensemble == "NPT":
            dynamiccmd = f"{phase_dynamic[phase]} {xyzfile} -key {keyfile} {liquidtotalstep} {liquidtimestep} {liquidwriteout} 4 {liquidT} {liquidP} > {logfile}"
          elif liquidensemble == "NVT":
            dynamiccmd = f"{phase_dynamic[phase]} {xyzfile} -key {keyfile} {liquidtotalstep} {liquidtimestep} {liquidwriteout} 2 {liquidT} > {logfile}"
          else:
            sys.exit(RED + "Error: only NPT or NVT ensemble is supported for Free energy simulations" + ENDC)
          if nodes == []:
            subprocess.run(dynamiccmd, shell=True)
            print(GREEN + ' [GOOD] Dynamic jobs submitted for {fname}' + ENDC)
          else:
            subprocess.run(f"nohup python {submitexe} -c '{dynamiccmd}' -n {nodes[i]} 2>err &", shell=True)
            print(GREEN + ' [GOOD] Dynamic jobs submitted for {fname}' + ENDC)
        if phase == 'gas':
          dynamiccmd = f"{phase_dynamic[phase]} {xyzfile} -key {keyfile} {gastotalstep} {gastimestep} {gaswriteout} 2 {gastemperature} > {logfile}"
          if nodes == []:
            subprocess.run(dynamiccmd, shell=True)
            print(GREEN + ' [GOOD] Dynamic jobs submitted for {fname}' + ENDC)
          else:
            subprocess.run(f"nohup python {submitexe} -c '{dynamiccmd}' -n {nodes[i]} 2>err &", shell=True)
            print(GREEN + ' [GOOD] Dynamic jobs submitted for {fname}' + ENDC)
      else:
        print(YELLOW + f" [Warning] {logfile} exists in {phase} folder for {fname}" + ENDC)
  return

def bar():
  print(YELLOW + " Checking the completeness of the MD trajectories, please wait... " + ENDC)
  if inputaction == "auto":
    proceed = False
    while not proceed: 
      proceed, phase_simtime = checkdynamic(liquidtotaltime, gastotaltime, phases, orderparams, homedir)
      if proceed:
        break
      now = datetime.now().strftime("%b %d %Y %H:%M:%S")
      print(YELLOW + f" [{now}] Waiting for dynamic jobs to finish ..." + ENDC)
      time.sleep(checktime)
  else:
    proceed, phase_simtime = checkdynamic(liquidtotaltime, gastotaltime, phases, orderparams, homedir)
  
  if proceed:
    for phase in phases:
      for i in range(0,len(orderparams)-1,1):
        elb0, vlb0 = orderparams[i]
        elb1, vlb1 = orderparams[i+1]
        fname0 = "%s-e%s-v%s"%(phase, "%03d"%int(elb0*100), "%03d"%int(vlb0*100))
        fname1 = "%s-e%s-v%s"%(phase, "%03d"%int(elb1*100), "%03d"%int(vlb1*100))
        arcfile0 = fname0 + ".arc"
        arcfile1 = fname1 + ".arc"
        if phase == 'liquid':
          printname = fname0
          outfile = fname0 + ".out"
          barfile = fname0 + ".bar"
          enefile = fname0 + ".ene"
          barcmd1 = f"{liquidbarexe} 1 {arcfile0} {liquidT} {arcfile1} {liquidT} N > {outfile} && "
          totalsnapshot = int(phase_simtime[phase]/liquidwriteout)
          startsnapshot = int(totalsnapshot/5.0) + 1
          barcmd2 = f"{liquidbarexe} 2 {barfile} {startsnapshot} {totalsnapshot} 1 {startsnapshot} {totalsnapshot} 1 > {enefile} "
          barstr = barcmd1 + barcmd2
          if nodes == []:
            if (not os.path.isfile(os.path.join(homedir, phase, barfile))):
              subprocess.run(cmdstr, shell=True)
            else:
              print(GREEN + f" [Warning] {barfile} exist in {phase} folder! If you want to overwrite, please set `overwritebar=True`" + ENDC)
          else:
            if (not os.path.isfile(os.path.join(homedir, phase, barfile))):
              phasedir = os.path.join(homedir, phase)
              os.chdir(phasedir)
              cmdstr = f"nohup python {submitexe} -c '{barstr}' -n {nodes[i]} 2>err &"
              subprocess.run(cmdstr, shell=True)
              print(GREEN + f" submitted {printname} on {nodes[i]}" + ENDC)
            else:
              print(GREEN + f" [Warning] {barfile} exist in {phase} folder! If you want to overwrite, please set `overwritebar:True` in settings.yaml" + ENDC)
        if phase == 'gas':
          printname = fname1
          outfile = fname1 + ".out"
          barfile = fname1 + ".bar"
          enefile = fname1 + ".ene"
          barcmd1 = f"{gasbarexe} 1 {arcfile1} {gastemperature} {arcfile0} {gastemperature} N > {outfile} && "
          totalsnapshot = int(phase_simtime[phase]/liquidwriteout)
          startsnapshot = int(totalsnapshot/5.0) + 1
          barcmd2 = f"{gasbarexe} 2 {barfile} {startsnapshot} {totalsnapshot} 1 {startsnapshot} {totalsnapshot} 1 > {enefile} "
          barstr = barcmd1 + barcmd2
          if nodes == []:
            if (not os.path.isfile(os.path.join(homedir, phase, barfile))):
              subprocess.run(barstr, shell=True)
            else:
              print(GREEN + f" [Warning] {barfile} exist in {phase} folder! If you want to overwrite, please set `overwritebar:True` in settings.yaml" + ENDC)
          else:
            if (not os.path.isfile(os.path.join(homedir, phase, barfile))):
              phasedir = os.path.join(homedir, phase)
              os.chdir(phasedir)
              cmdstr = f"nohup python {submitexe} -c '{barstr}' -n {nodes[i]} 2>err &"
              subprocess.run(cmdstr, shell=True)
              print(GREEN + f" submitted {printname} on {nodes[i]}" + ENDC)
            else:
              print(GREEN + f" [Warning] {barfile} exist in {phase} folder! If you want to overwrite, please set `overwritebar:True` in settings.yaml" + ENDC)
  return

def result():
  print(YELLOW + " Checking the completeness of the BAR analysis" + ENDC)
  if inputaction == 'auto':
    proceed = False
    while not proceed: 
      proceed, gasperturbsteps, gasenes, liquidperturbsteps, liquidenes = checkbar(phases, orderparams, homedir)
      if proceed:
        break
      now = datetime.now().strftime("%b %d %Y %H:%M:%S")
      print(YELLOW + f" [{now}] Waiting for bar jobs to finish ..." + ENDC)
      time.sleep(checktime)
  else:
    proceed, gasperturbsteps, gasenes, liquidperturbsteps, liquidenes = checkbar(phases, orderparams, homedir)
  if proceed:
    FEgas = []
    Errgas = []
    FEliquid = []
    Errliquid = []
    for gasene in gasenes:
      find = False
      for line in open(os.path.join(homedir, 'gas', gasene)).readlines():
        if "Free Energy via BAR Iteration" in line:
          find = True
          FEgas.append(float(line.split()[-4]))
          Errgas.append(float(line.split()[-2]))
        if (not find) and ("Free Energy via BAR Bootstrap" in line):
          find = True
          FEgas.append(float(line.split()[-4]))
          Errgas.append(float(line.split()[-2]))
    
    for liquidene in liquidenes:
      find = False
      for line in open(os.path.join(homedir, 'liquid', liquidene)).readlines():
        if "Free Energy via BAR Iteration" in line:
          find = True
          FEliquid.append(float(line.split()[-4]))
          Errliquid.append(float(line.split()[-2]))
        if (not find) and ("Free Energy via BAR Bootstrap" in line):
          find = True
          FEliquid.append(float(line.split()[-4]))
          Errliquid.append(float(line.split()[-2]))
    
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
    print(GREEN + "%40s%25.4f%20.4f"%("SUMMARY OF THE TOTAL FREE ENERGY", totFE, totErr) + ENDC)
    fo.write("%40s%25.4f%20.4f\n"%("SUMMARY OF THE TOTAL FREE ENERGY", totFE, totErr))
    fo.close()
  return

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('act', help = "Actions to take.", choices = ['setup', 'dynamic', 'bar', 'result', 'auto'], type = str.lower) 
  
  global inputaction
  inputaction = vars(parser.parse_args())['act']
  
  global rootdir, homedir
  rootdir = os.path.join(os.path.split(__file__)[0])
  homedir = os.getcwd()

  if not os.path.isfile(os.path.join(homedir, "settings.yaml")):
    sys.exit(RED + "Please provide 'settings.yaml' file; " + ENDC + GREEN + f"An example is here for you {os.path.join(rootdir, 'dat', 'settings.yaml')}" + ENDC)
  else:
    with open('settings.yaml') as f:
      FEsimsettings = yaml.load(f, Loader=yaml.FullLoader)
  
  global prm, checktime
  lig = FEsimsettings['gas_xyz'] 
  box = FEsimsettings['box_xyz'] 
  prm = FEsimsettings["tinker_prm"]
  hfe = FEsimsettings["hydration"] 
  nodelist = FEsimsettings["node_list"]
  checktime = FEsimsettings["checktime"]

  global orderparams
  orderparams = []
  lambdawindow = FEsimsettings["lambda_window"]
  if lambdawindow == "COURSER":
    orderprmfile = os.path.join(rootdir, "dat", "orderparams_courser")
  else:
    orderprmfile = os.path.join(rootdir, "dat", "orderparams")
  for line in open(orderprmfile).readlines():
    if "#" not in line:
      d = line.split()
      orderparams.append([float(d[0]), float(d[1])])
  
  global liquidkeylines, gaskeylines
  liquidkeylines = open(os.path.join(rootdir, "dat", "liquid.key")).readlines()
  gaskeylines = open(os.path.join(rootdir, "dat", "gas.key")).readlines()
  
  global phases
  phases = ['liquid']
  if hfe:
    phases.append('gas')
  for phase in phases:
    if not os.path.isdir(f"{phase}"):
      os.system(f"mkdir {phase}")
  
  global submitexe
  submitexe = os.path.join(rootdir, "utils", "submit.py")
  
  global nodes
  nodes = FEsimsettings["node_list"]
  if nodes == None:
    nodes = []

  global natom, phase_xyz, phase_key
  natom = int(open(lig).readlines()[0].split()[0])
  phase_xyz = {'liquid':box, 'gas':lig}
  phase_key = {'liquid':liquidkeylines, 'gas':gaskeylines}

  global liquidtotaltime, liquidtimestep, liquidwriteout, liquidtotalstep
  global liquidT, liquidP, liquidensemble
  liquidtotaltime = FEsimsettings["liquid_md_total_time"] 
  liquidtimestep = FEsimsettings["liquid_md_time_step"] 
  liquidwriteout = FEsimsettings["liquid_md_write_freq"] 
  liquidT = FEsimsettings["liquid_md_temperature"] 
  liquidP = FEsimsettings["liquid_md_pressure"]
  liquidensemble = FEsimsettings["liquid_md_ensemble"].upper()
  liquidtotalstep = int((1000000.0*liquidtotaltime)/liquidtimestep)

  global gastotaltime, gastimestep, gaswriteout, gastemperature, gastotalstep
  gastotaltime = FEsimsettings["gas_md_total_time"] 
  gastimestep = FEsimsettings["gas_md_time_step"] 
  gastemperature = FEsimsettings["gas_md_temperature"] 
  gaswriteout = FEsimsettings["gas_md_write_freq"]
  gastotalstep = int((1000000.0*gastotaltime)/gastimestep)

  global liquidmdexe, liquidbarexe, gasmdexe, gasbarexe
  liquidmdexe = "/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda10.2/dynamic9.sh"
  liquidbarexe = "/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda10.2/bar9.sh"
  gasmdexe = "/home/liuchw/Softwares/tinkers/Tinker-latest/source-C8/dynamic.x" 
  gasbarexe = "/home/liuchw/Softwares/tinkers/Tinker-latest/source-C8/bar.x"
  
  overwritebar = FEsimsettings["overwritebar"]
  if overwritebar:
    os.system(f"rm -f */*.bar */*.ene")

  actions = {'setup':setup, 'dynamic':dynamic, 'bar':bar, 'result':result}
  if inputaction in actions.keys():
    actions[inputaction]()
  else:
    for action in actions.keys():
      actions[action]()
  return

if __name__ == "__main__":
  main()
