
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

# color
RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\033[93m'

def check_cpu_avail(node, nproc_required):
  
  # assume fully occupied
  tot_occ = 64.0 
  
  try:
    # occupied nproc
    cmd = f'ssh {node} "top -n1 -b" | grep " R \| S " '
    sp_ret = subprocess.check_output(cmd, timeout=10.0, shell=True).decode("utf-8").split('\n')[:-1]
    tot_occ = 0 
    for r in sp_ret:
      if 'R'  in r:
        occ = r.split('R')[-1].split()[0]
        if occ.replace('.', '', 1).isdigit():
          tot_occ += float(occ)/100.0
      if 'S'  in r:
        occ = r.split('S')[-1].split()[0]
        if occ.replace('.', '', 1).isdigit():
          tot_occ += float(occ)/100.0
    tot_occ = round(tot_occ)
  
  except:
    pass
  
  # assume no CPUs 
  nproc = 0
  try:
    # total nproc
    cmd = f'ssh {node} "nproc" '
    sp_ret = subprocess.check_output(cmd, timeout=10.0, shell=True).decode("utf-8").split('\n')[0]
    if sp_ret != '':
      nproc = int(sp_ret) 
  except:
    pass

  # limit the usage to 80%
  nproc = int(nproc*0.8)
  
  # available nproc
  avail = False
  avail_nproc = nproc - tot_occ
  if avail_nproc > nproc_required:
    avail = True
  return avail

def check_gpu_avail(node):
  sp_ret = '  '
  try:
    # occupied nproc
    cmd = f'ssh {node} "nvidia-smi" 2>/dev/null'
    sp_ret = subprocess.check_output(cmd, timeout=10.0, shell=True).decode("utf-8").split('\n')[:-1]
  except:
    pass

  tot_cards = []
  occ_cards = []

  for r in sp_ret:
    if 'N/A' in r:
      if r.split()[1].isdigit():
        tot_cards.append(r.split()[1])
    # tinker9/dynamic9/bar9 is for Tinker9
    # dynamic is for openmm
    # gmx is for gromacs
    if ('tinker9' in r) or ('dynamic' in r) or ('bar9' in r) or ('dynamic9' in r) or ('gmx' in r):
      occ_cards.append(r.split()[1])
  tot_cards = list(set(tot_cards)) 
  ava_cards = tot_cards
  if occ_cards != []:
    ava_cards = []
    for c in tot_cards:
      if c not in occ_cards:
        ava_cards.append(c)
  return ava_cards 

def submit_jobs(jobcmds, jobtype):
  rtx30s = ['node36', 'node54', 'node102', 'node193', 'node206', \
            'black', 'bme-black', 'pluto', 'bme-pluto', 'bme-mercury', 'mercury', \
            'bme-mars', 'mars', 'bme-venus', 'venus', 'node103', 'node104', 'node105', \
            'bme-rna', 'rna', 'bme-jupiter', 'jupiter', 'bme-earth', 'earth']
  njob_pointer = 0
  if jobtype == "CPU":
    for i in range(len(cpu_node_list)):
      if njob_pointer >= len(jobcmds): break
      if check_cpu_avail(cpu_node_list[i], nproc):
        cmdstr = f'ssh {cpu_node_list[i]} "' + jobcmds[njob_pointer] + ' " &'
        subprocess.run(cmdstr, shell=True)
        jobcmds[njob_pointer] = 'x'
        print(f"{cmdstr}")
        njob_pointer += 1
  else:
    for i in range(len(gpu_node_list)): 
      if njob_pointer >= len(jobcmds): break
      ava_cards = check_gpu_avail(gpu_node_list[i]) 
      if ava_cards != []:
        for card in ava_cards:
          cuda_device = f'export CUDA_VISIBLE_DEVICES="{card}"'
          pci_bus_id = 'export CUDA_DEVICE_ORDER=PCI_BUS_ID'
          if njob_pointer < len(jobcmds):
            cmdstr = f"ssh {gpu_node_list[i]} '{pci_bus_id}; {cuda_device}; {jobcmds[njob_pointer]} ' &"
            if gpu_node_list[i] in rtx30s:
              if 'cuda10.2' in cmdstr:
                cmdstr = cmdstr.replace('cuda10.2', 'cuda11.2')
            subprocess.run(cmdstr, shell=True)
            jobcmds[njob_pointer] = 'x'
            print(f"   --> {cmdstr}")
            njob_pointer += 1
    # wait for 15 sec. to let job appear on a node
    # i.e., shown by nvidia-smi command
    # to avoid submitting multiple jobs on one GPU card
    time.sleep(15.0)
     
  # return the remainig jobcmds
  tmp = [] 
  for jobcmd in jobcmds:
    if jobcmd != 'x':
      tmp.append(jobcmd)
  jobcmds = tmp
  return jobcmds

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-x', dest = 'jobshs',  nargs='+', help = "Scripts to run", default = []) 
  parser.add_argument('-p', dest = 'path',  help = "Working directory", default=None) 
  parser.add_argument('-c', dest = 'jobcmds',  nargs='+', help = "Commands to run", default = []) 
  parser.add_argument('-t', dest = 'type',  help = "Job type", choices =['CPU', 'GPU'], required=True, type = str.upper) 
  parser.add_argument('-n', dest = 'nproc',  help = "Nproc requested", default=2, type=int) 
  parser.add_argument('-nodes', dest = 'nodes',  nargs='+', help = "node list") 
  args = vars(parser.parse_args())
  jobshs = args["jobshs"]
  jobcmds = args["jobcmds"]
  jtyp = args["type"]
  global nproc,path
  path = args["path"]
  nproc = args["nproc"]
  nodes = args["nodes"]  

  global gpu_node_list
  global cpu_node_list
  
  workingdir = os.getcwd()
  if path != None:
    workingdir = path
  if jobcmds != []:
    for i in range(len(jobcmds)):
      jobcmds[i] = f"cd {workingdir}; {jobcmds[i]}"
  if jobshs != []:
    jobcmds = []
    for jobsh in jobshs:
      jobcmds.append(f'cd {workingdir}; sh {jobsh}')
  
  print(GREEN + f"   === Submitting {jtyp} Jobs to Ren Lab Clusters === " + ENDC)
  gpu_node_list = nodes
  cpu_node_list = nodes
  jobcmds = submit_jobs(jobcmds, jtyp)
  while jobcmds != []:
    time.sleep(5.0)
    jobcmds = submit_jobs(jobcmds, jtyp)