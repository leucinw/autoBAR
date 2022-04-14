
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
    if not (ignoregas == 1 and phase == 'gas'):
      for i in range(len(orderparams)):
        elb, vlb = orderparams[i]
        fname = "%s-e%s-v%s"%(phase, "%03d"%int(elb*100), "%03d"%int(vlb*100))
        with open(fname + ".key", 'w') as fw:
          for line in phase_key[phase]:
            if 'parameters' in line.lower():
              if (elb*vlb > 1.0):
                line = f'parameters     ../{perturbprm}\n'
              else:
                line = f'parameters     ../{prm}\n'
            if 'polar-eps ' in line.lower():
              line = f'polar-eps {polareps} \n'
            fw.write(line)
          fw.write('\n')
          fw.write(f'ligand -1 {natom}\n')
          if (elb*vlb > 1.0):
            elb = 1.0
            vlb = 1.0
          fw.write(f'ele-lambda {elb}\n')
          fw.write(f'vdw-lambda {vlb}\n')
            
        linkxyz = f"ln -sf {homedir}/{phase_xyz[phase]} {homedir}/{phase}/{fname}.xyz"
        movekey = f"mv {fname}.key ./{phase}"
        subprocess.run(linkxyz, shell=True)
        subprocess.run(movekey, shell=True)
  print(GREEN + ' [GOOD] BAR simulation files generated!' + ENDC)
  return

def dynamic():
  liquidcmds = []
  gascmds = []
  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    os.system(f"rm -rf {phasedir}/*.sh")
    os.chdir(phasedir)
    for i in range(len(orderparams)):
      elb, vlb = orderparams[i]
      fname = "%s-e%s-v%s"%(phase, "%03d"%int(elb*100), "%03d"%int(vlb*100))
      if (elb*vlb > 1.0) and copyarcforperturb:
        print(YELLOW + f" Skipping dynamic for {fname} " + ENDC)
        break
      xyzfile = fname + ".xyz"
      keyfile = xyzfile.replace("xyz", "key")
      logfile = xyzfile.replace("xyz", "log")
      if os.path.isfile(logfile):
        print(YELLOW + f" [Warning] {logfile} exists in {phase} folder for {fname}" + ENDC)
      else:
        if phase == 'liquid':
          if liquidensemble == "NPT":
            dynamiccmd = f" 'cd {phasedir}; {phase_dynamic[phase]} {xyzfile} -key {keyfile} {liquidtotalstep} {liquidtimestep} {liquidwriteout} 4 {liquidT} {liquidP} > {logfile}' "
            liquidcmds.append(dynamiccmd)
          if liquidensemble == "NVT":
            dynamiccmd = f" 'cd {phasedir}; {phase_dynamic[phase]} {xyzfile} -key {keyfile} {liquidtotalstep} {liquidtimestep} {liquidwriteout} 2 {liquidT} > {logfile}' "
            liquidcmds.append(dynamiccmd)
        if phase == 'gas':
          dynamiccmd = f" 'cd {phasedir}; {phase_dynamic[phase]} {xyzfile} -key {keyfile} {gastotalstep} {gastimestep} {gaswriteout} 2 {gastemperature} > {logfile}' "
          gascmds.append(dynamiccmd)
    
  # submit jobs to clusters 
  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    os.chdir(phasedir)
    if phase == 'gas':
      cmdstr = " ".join(gascmds)
      os.system(f"python {submitexe} -c {cmdstr} -t CPU -nodes {' '.join(nodes)} -n 4")
    if phase == 'liquid':
      cmdstr = " ".join(liquidcmds)
      os.system(f"python {submitexe} -c {cmdstr} -t GPU -nodes {' '.join(nodes)}")
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
      proceed, phase_simtime = checkdynamic(liquidtotaltime, gastotaltime, phases, checkorderparams, homedir)
      if proceed:
        break
      now = datetime.now().strftime("%b %d %Y %H:%M:%S")
      print(YELLOW + f" [{now}] Waiting for dynamic jobs to finish ..." + ENDC)
      time.sleep(checkingtime)
  else:
    proceed, phase_simtime = checkdynamic(liquidtotaltime, gastotaltime, phases, checkorderparams, homedir)
    
  if proceed:
    liquidcmds = []
    gascmds = []
    for phase in phases:
      phasedir = os.path.join(homedir, phase)
      os.chdir(phasedir)
      for i in range(0,len(orderparams)-1,1):
        elb0, vlb0 = orderparams[i]
        elb1, vlb1 = orderparams[i+1]
        fname0 = "%s-e%s-v%s"%(phase, "%03d"%int(elb0*100), "%03d"%int(vlb0*100))
        fname1 = "%s-e%s-v%s"%(phase, "%03d"%int(elb1*100), "%03d"%int(vlb1*100))
        arcfile0 = fname0 + ".arc"
        arcfile1 = fname1 + ".arc"
        if (elb1*vlb1 > 1.0) and copyarcforperturb:
          linkarc = f"ln -sf {arcfile0} {arcfile1}"
          os.system(linkarc)
        if phase == 'liquid':
          printname = fname0
          outfile = fname0 + ".out"
          barfile = fname0 + ".bar"
          enefile = fname0 + ".ene"
          barcmd1 = f"cd {phasedir}; {liquidbarexe} 1 {arcfile0} {liquidT} {arcfile1} {liquidT} N > {outfile} && "
          totalsnapshot = int(phase_simtime[phase]/liquidwriteout)
          startsnapshot = int(totalsnapshot/5.0) + 1
          barcmd2 = f"{liquidbarexe} 2 {barfile} {startsnapshot} {totalsnapshot} 1 {startsnapshot} {totalsnapshot} 1 > {enefile} "
          barstr = f" '{barcmd1} {barcmd2} ' "
          if (not os.path.isfile(os.path.join(homedir, phase, barfile))):
            liquidcmds.append(barstr)
          else:
            print(GREEN + f" [Warning] {barfile} exist in {phase} folder! " + ENDC)
        if phase == 'gas':
          printname = fname1
          outfile = fname1 + ".out"
          barfile = fname1 + ".bar"
          enefile = fname1 + ".ene"
          barcmd1 = f"cd {phasedir}; {gasbarexe} 1 {arcfile1} {gastemperature} {arcfile0} {gastemperature} N > {outfile} && "
          totalsnapshot = int(phase_simtime[phase]/liquidwriteout)
          startsnapshot = int(totalsnapshot/5.0) + 1
          barcmd2 = f"{gasbarexe} 2 {barfile} {startsnapshot} {totalsnapshot} 1 {startsnapshot} {totalsnapshot} 1 > {enefile} "
          barstr = f" '{barcmd1} {barcmd2} ' "
          if gastotaltime != 0:
            if (not os.path.isfile(os.path.join(homedir, phase, barfile))):
              gascmds.append(barstr)
            else:
              print(GREEN + f" [Warning] {barfile} exist in {phase} folder!" + ENDC)
  
    # submit jobs to clusters 
    for phase in phases:
      phasedir = os.path.join(homedir, phase)
      os.chdir(phasedir)
      if phase == 'gas':
        cmdstr = " ".join(gascmds)
        os.system(f"python {submitexe} -c {cmdstr} -t CPU -nodes {' '.join(nodes)} -n 4")
      if phase == 'liquid':
        cmdstr = " ".join(liquidcmds)
        os.system(f"python {submitexe} -c {cmdstr} -t GPU -nodes {' '.join(nodes)}")
  return

def result():
  print(YELLOW + " Checking the completeness of the BAR analysis" + ENDC)
  if inputaction == 'auto':
    proceed = False
    while not proceed:
      proceed, gasperturbsteps, gasenes, liquidperturbsteps, liquidenes = checkbar(phases, orderparams, homedir, ignoregas)
      if proceed:
        break
      now = datetime.now().strftime("%b %d %Y %H:%M:%S")
      print(YELLOW + f" [{now}] Waiting for bar jobs to finish ..." + ENDC)
      time.sleep(checkingtime)
  else:
    proceed, gasperturbsteps, gasenes, liquidperturbsteps, liquidenes = checkbar(phases, orderparams, homedir, ignoregas)
  
  if proceed:
    if copyarcforperturb:
      for phase in phases:
        rmcmd = f"rm -f {phase}/*200.arc"
        os.system(rmcmd)
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
      FEsimsettings = yaml.load(f, Loader=yaml.FullLoader)
  global prm, checkingtime
  lig = FEsimsettings['gas_xyz'] 
  box = FEsimsettings['box_xyz'] 
  prm = FEsimsettings["parameters"]
  checkingtime = FEsimsettings["checking_time"]
  global natom
  natom = int(open(lig).readlines()[0].split()[0])
  global polareps 
  polareps = FEsimsettings["polar_eps"]
  global nodes
  nodes = FEsimsettings["node_list"]
  if nodes == None:
    sys.exit(RED + "[Error] node_list must be provided" + ENDC) 
  global submitexe
  submitexe = os.path.join(rootdir, "utils", "submitTinker.py")
  global orderparams
  orderparams = []
  lambdawindow = FEsimsettings["lambda_window"].upper()
  if lambdawindow == "COURSER":
    orderprmfile = os.path.join(rootdir, "dat", "orderparams_courser")
  else:
    orderprmfile = os.path.join(rootdir, "dat", "orderparams_default")
  for line in open(orderprmfile).readlines():
    if "#" not in line:
      d = line.split()
      orderparams.append([float(d[0]), float(d[1])])
  
  global perturbprm
  global copyarcforperturb
  copyarcforperturb = False
  perturbprm = prm + "_"
  if os.path.isfile(os.path.join(homedir, perturbprm)):
    orderparams.append([2.0, 2.0])
    copyarcforperturb = FEsimsettings["copy_arc_for_perturb"]

  # liquid phase specific settings 
  global phases, liquidkeylines
  phases = ['liquid']
  liquidkeylines = open(os.path.join(rootdir, "dat", "liquid.key")).readlines()
  global liquidmdexe, liquidbarexe
  liquidmdexe = "/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda10.2/dynamic9.sh"
  liquidbarexe = "/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda10.2/bar9.sh"
  global phase_xyz, phase_key, phase_dynamic
  phase_xyz = {'liquid':box}
  phase_key = {'liquid':liquidkeylines}
  phase_dynamic = {'liquid':liquidmdexe}
  global liquidtotaltime, liquidtimestep, liquidwriteout, liquidtotalstep
  global liquidT, liquidP, liquidensemble
  liquidtotaltime = FEsimsettings["liquid_md_total_time"] 
  liquidtimestep = FEsimsettings["liquid_md_time_step"] 
  liquidwriteout = FEsimsettings["liquid_md_write_freq"] 
  liquidT = FEsimsettings["liquid_md_temperature"] 
  liquidP = FEsimsettings["liquid_md_pressure"]
  liquidensemble = FEsimsettings["liquid_md_ensemble"].upper()
  liquidtotalstep = int((1000000.0*liquidtotaltime)/liquidtimestep)
  if not os.path.isdir('liquid'):
    os.system("mkdir liquid")

  # gas phase specific settings
  global ignoregas
  global gaskeylines,gasmdexe,gasbarexe
  gaskeylines = open(os.path.join(rootdir, "dat", "gas.key")).readlines()
  gasmdexe = "/home/liuchw/Softwares/tinkers/Tinker-latest/source-C8/dynamic.x" 
  gasbarexe = "/home/liuchw/Softwares/tinkers/Tinker-latest/source-C8/bar.x"
  global gastotaltime, gastimestep, gaswriteout, gastemperature, gastotalstep
  gastotaltime = FEsimsettings["gas_md_total_time"] 
  gastimestep = FEsimsettings["gas_md_time_step"] 
  gastemperature = FEsimsettings["gas_md_temperature"] 
  gaswriteout = FEsimsettings["gas_md_write_freq"]
  gastotalstep = int((1000000.0*gastotaltime)/gastimestep)
  phase_xyz = {'liquid':box, 'gas':lig}
  phase_key = {'liquid':liquidkeylines, 'gas':gaskeylines}
  phase_dynamic = {'liquid':liquidmdexe, 'gas':gasmdexe}
  
  # check xyz files
  lines = open(box).readlines()
  natomliquid = int(lines[0].split()[0])
  if natomliquid != len(lines)-2:
    sys.exit(RED + f"[Error] Please provide box info in {box}" + ENDC)
  [a,b,c] = lines[1].split()[0:3]
  if min([float(a), float(b), float(c)]) < 30.0:
    sys.exit(RED + f"[Error] Please provide a bigger box (>30*30*30)" + ENDC)

  natomgas = int(open(lig).readlines()[0].split()[0])
  if natomgas == 1:
    gastotaltime = 0.0
    print(YELLOW + "[Warning] I set the simulation time to 0 since it is a single ion/atom" + ENDC)
  if gastotaltime == 0.0:
    ignoregas = 1
  else:
    ignoregas = 0
  if (not os.path.isdir('gas')) and (ignoregas == 0):
    os.system("mkdir gas")
    phases.append('gas')
  
  actions = {'setup':setup, 'dynamic':dynamic, 'bar':bar, 'result':result}
  if inputaction in actions.keys():
    actions[inputaction]()
  elif inputaction == 'auto': 
    for action in actions.keys():
      actions[action]()
  else:
    sys.exit(RED + f"[Error] {inputaction} does not supported!" + ENDC)
  return

if __name__ == "__main__":
  main()
