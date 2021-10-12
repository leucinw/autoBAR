
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os

RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\33[93m'

def checkdynamic(liquidtotaltime, gastotaltime, phases, orderparams, homedir):
  statuslist = [] 
  phase_simtime = {'liquid':1000.0*liquidtotaltime, 'gas':1000.0*gastotaltime}
  for phase in phases:
    for i in range(len(orderparams)):
      elb, vlb = orderparams[i]
      fname = "%s-e%s-v%s"%(phase, "%03d"%int(elb*100), "%03d"%int(vlb*100))
      logfile = fname + ".log"
      if os.path.isfile(os.path.join(homedir, phase, logfile)):
        lines = open(os.path.join(homedir, phase, logfile)).readlines()
        simtime = 0.0
        for line in lines[-20:]:
          if "Current Time" in line:
            simtime = float(line.split()[2])
        if (simtime == phase_simtime[phase]):
          print(GREEN + "  [" + fname + f"]: {simtime}/{phase_simtime[phase]} ps!" + ENDC)
          statuslist.append(True)
        else:
          print(RED + "  [" + fname + f"]: {simtime}/{phase_simtime[phase]} ps!" + ENDC)
          statuslist.append(False) 
      else:
        statuslist.append(False)
  completed = all(statuslist)
  return completed, phase_simtime

def checkbar(phases, orderparams, homedir):
  statuslist = []
  gasenes = []
  liquidenes = []
  gasperturbsteps = []
  liquidperturbsteps = []
  for phase in phases:
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
