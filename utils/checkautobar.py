
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import subprocess

RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\33[93m'

def checkdynamic(liquidtotalsnapshot, gastotalsnapshot, phases, orderparams, homedir):
  statuslist = []
  phase_simsnapshot = {'liquid': liquidtotalsnapshot, 'gas':gastotalsnapshot}
  for phase in phases:
    for i in range(len(orderparams)):
      elb, vlb = orderparams[i]
      fname = "%s-e%s-v%s"%(phase, "%03d"%int(elb*100), "%03d"%int(vlb*100))
      arcfile = fname + ".arc"
      if os.path.isfile(os.path.join(homedir, phase, arcfile)):
        simsnapshot = 0
        # get the simulation snapshots here
        cmd = f"head -n1 {os.path.join(homedir, phase, arcfile)} "
        checkstr = subprocess.check_output(cmd, shell=True).decode("utf-8")
        firstline = '^' + checkstr.replace('\n', '$')
        cmd = f"LC_ALL=C fgrep '{firstline}' {os.path.join(homedir, phase, arcfile)} | wc -l"
        checkstr = subprocess.check_output(cmd, shell=True).decode("utf-8")
        simsnapshot = int(checkstr)
        if (simsnapshot == phase_simsnapshot[phase]):
          per = int(simsnapshot/phase_simsnapshot[phase]*100)
          print(GREEN + f"{fname:>20s}: " + u'\u2584'*(int(per/2))  + f" [{per:>3d}%]" + ENDC)
          statuslist.append(True)
        else:
          per = int(simsnapshot/phase_simsnapshot[phase]*100)
          print(YELLOW + f"{fname:>20s}: " + u'\u2584'*(int(per/2))  + f" [{per:>3d}%]" + ENDC)
          statuslist.append(False) 
      else:
        if gastotalsnapshot == 0:
          statuslist.append(True)
        else:
          print(RED + f"{arcfile} does not exist" + ENDC)
          statuslist.append(False)
  completed = all(statuslist)
  return completed, phase_simsnapshot

def checkbar(phases, orderparams, homedir, ignoregas):
  statuslist = []
  gasenes = []
  liquidenes = []
  gasperturbsteps = []
  liquidperturbsteps = []
  checkphases = phases
  if ignoregas == 1:
    checkphases = ['liquid']
  
  for phase in checkphases:
    for i in range(0,len(orderparams)-1,1):
      elb0, vlb0 = orderparams[i]
      elb1, vlb1 = orderparams[i+1]
      fname0 = "%s-e%s-v%s"%(phase, "%03d"%int(elb0*100), "%03d"%int(vlb0*100))
      fname1 = "%s-e%s-v%s"%(phase, "%03d"%int(elb1*100), "%03d"%int(vlb1*100))
      if phase == 'gas':
        enefile = fname1 + ".ene"
        fname = fname1
        gasenes.append(enefile)
        gasperturbsteps.append([fname1, fname0])
      if phase == 'liquid':
        enefile = fname0 + ".ene"
        fname = fname0
        liquidenes.append(enefile)
        liquidperturbsteps.append([fname0, fname1])
      if not os.path.isfile(os.path.join(homedir, phase, enefile)):
        print(RED + " " + fname + f": free energy file (.ene) not found!" + ENDC)
        statuslist.append(False) 
      else:
        lines = open(os.path.join(homedir, phase, enefile)).readlines()
        barfinished = False
        for line in lines[-5:]:
          if "BAR Estimate of -T*dS" in line:
            barfinished = True
        statuslist.append(barfinished)
  finished = statuslist.count(True)
  targeted = len(statuslist)
  if finished != targeted:
    print(YELLOW + f" Finished {finished} out of {targeted} bar analysis steps" + ENDC)
  else:
    print(GREEN + f" Finished {finished} out of {targeted} bar analysis steps" + ENDC)
  completed = all(statuslist)
  return completed, gasperturbsteps, gasenes, liquidperturbsteps, liquidenes
