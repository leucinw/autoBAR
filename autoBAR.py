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
from pathlib import Path
from dataclasses import dataclass, field
from typing import List
import numpy as np
import ruamel.yaml as yaml
from datetime import datetime
from utils.checkautobar import *
from utils.elescale import *

# ANSI color codes for terminal output
RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\33[93m'


@dataclass
class Config:
  """Holds all run-time configuration for an autoBAR simulation."""
  # paths
  rootdir: str = ""
  homedir: str = ""
  prm: str = ""
  tinkerenv: str = ""
  submitexe: str = ""
  # phase data
  phases: List[str] = field(default_factory=list)
  phase_xyz: dict = field(default_factory=dict)
  phase_key: dict = field(default_factory=dict)
  phase_dynamic: dict = field(default_factory=dict)
  phase_minimize: dict = field(default_factory=dict)
  # lambda windows
  orderparams: List[list] = field(default_factory=list)
  copyarcforperturb: bool = False
  manualelescale: bool = False
  # liquid MD
  liquidtotaltime: float = 0.0
  liquidtimestep: float = 0.0
  liquidwriteout: float = 0.0
  liquidtotalstep: int = 0
  liquidtotalsnapshot: int = 0
  liquidT: float = 298.0
  liquidP: float = 1.0
  liquidensemble: str = "NPT"
  liquidmdexe: str = "$DYNAMIC9"
  liquidbarexe: str = "$BAR9"
  # gas MD
  gastotaltime: float = 0.0
  gastimestep: float = 0.0
  gaswriteout: float = 0.0
  gastemperature: float = 298.0
  gastotalstep: int = 0
  gastotalsnapshot: int = 0
  gasmdexe: str = "$DYNAMIC8"
  gasbarexe: str = "$BAR8"
  # misc
  natom: int = 0
  nodes: List[str] = field(default_factory=list)
  verbose: int = 1
  ignoregas: int = 0
  checkingtime: float = 60.0
  inputaction: str = ""


def parse_ene_file(path: str) -> tuple:
  """Return (free_energy, error) from a Tinker .ene file.

  Tries 'BAR Iteration' first, then falls back to 'BAR Bootstrap'.
  Raises ValueError if neither is found.
  """
  iteration_result = None
  bootstrap_result = None
  with open(path) as fh:
    for line in fh:
      if "Free Energy via BAR Iteration" in line:
        parts = line.split()
        iteration_result = (float(parts[-4]), float(parts[-2]))
      elif bootstrap_result is None and "Free Energy via BAR Bootstrap" in line:
        parts = line.split()
        bootstrap_result = (float(parts[-4]), float(parts[-2]))
  if iteration_result is not None:
    return iteration_result
  if bootstrap_result is not None:
    return bootstrap_result
  raise ValueError(f"No BAR free energy found in {path}")


def setup(cfg):
  """Generate Tinker input files (XYZ, KEY) for all lambda windows."""
  for phase in cfg.phases:
    xyzfile = cfg.phase_xyz[phase]
    keylines = cfg.phase_key[phase]
    with open(f"{phase}.key", 'w') as f:
      for line in keylines:
        if 'parameters' in line.lower():
          line = f'parameters     {cfg.prm}\n'
        f.write(line)
    gasminsh = 'gas-min.sh'
    if phase == 'gas':
      with open(gasminsh, 'w') as f:
        f.write(f'source {cfg.tinkerenv}\n')
        f.write(f'{cfg.phase_minimize[phase]} {xyzfile} -key gas.key 0.1 > gas-min.log \n')
        f.write(f'wait\nmv {cfg.phase_xyz[phase]}_2 {cfg.phase_xyz[phase]}\n')
      if not Path("gas-min.log").is_file():
        subprocess.run(f"sh {gasminsh}", shell=True)
    liquidminsh = 'liquid-min.sh'
    if phase == 'liquid':
      with open(liquidminsh, 'w') as f:
        f.write(f'source {cfg.tinkerenv}\n')
        f.write(f'{cfg.phase_minimize[phase]} {xyzfile} -key liquid.key 0.2 > liquid-min.log\n')
        f.write(f'wait\nmv {cfg.phase_xyz[phase]}_2 {cfg.phase_xyz[phase]}\n')
      if not Path("liquid-min.log").is_file():
        subprocess.run(f"sh {liquidminsh}", shell=True)
    subprocess.run(f"rm -f {phase}/*.xyz {phase}/*.key", shell=True)
    if not (cfg.ignoregas == 1 and phase == 'gas'):
      prmfiles = []
      for i in range(len(cfg.orderparams)):
        elb, vlb = cfg.orderparams[i]
        fname = lambda_fname(phase, elb, vlb)
        with open(fname + ".key", 'w') as fw:
          for line in keylines:
            if 'parameters' in line.lower():
              if (elb*vlb > 1.0):
                assert elb == vlb, RED + f" Error: lambdas greater than 1 but not the same: {elb}, {vlb} " + ENDC
                idx = int(round(elb*100)/10) - 10
                perturbprm = f"{cfg.prm}_{idx:02d}"
                line = f'parameters     {cfg.homedir}/{phase}/{perturbprm}\n'
                if perturbprm not in prmfiles:
                  prmfiles.append(perturbprm)
              else:
                line = f'parameters     {cfg.homedir}/{phase}/{cfg.prm}\n'
                if cfg.prm not in prmfiles:
                  prmfiles.append(cfg.prm)
            fw.write(line)
          fw.write('\n')
          fw.write(f'ligand -1 {cfg.natom}\n')
          if (elb*vlb > 1.0):
            elb = 1.0
            vlb = 1.0

          if cfg.manualelescale:
            fw.write('ele-lambda 1.0\n')
          else:
            fw.write(f'ele-lambda {elb}\n')
          fw.write(f'vdw-lambda {vlb}\n')

          if cfg.manualelescale:
            scaledprms = scaledownele(cfg.phase_xyz['gas'], cfg.prm, elb)
            fw.write(f'\n# scale down electrostatic related parameters by {elb}\n')
            for s in scaledprms:
              fw.write(f'{s}\n')

        linkxyz = f"ln -sf {cfg.homedir}/{cfg.phase_xyz[phase]} {cfg.homedir}/{phase}/{fname}.xyz"
        movekey = f"mv {fname}.key ./{phase}"
        subprocess.run(linkxyz, shell=True)
        subprocess.run(movekey, shell=True)

      # make copy of prmfile
      for prmfile in prmfiles:
        copyprm = f"cp {prmfile}  ./{phase}"
        subprocess.run(copyprm, shell=True)
  print(GREEN + ' [GOOD] BAR simulation files generated!' + ENDC)
  return


def dynamic(cfg):
  """Write MD shell scripts and submit them to the compute cluster."""
  liquidshs = []
  gasshs = []
  for phase in cfg.phases:
    phasedir = Path(cfg.homedir) / phase
    subprocess.run(f"rm -rf {phasedir}/*.sh", shell=True)
    os.chdir(phasedir)
    for i in range(len(cfg.orderparams)):
      elb, vlb = cfg.orderparams[i]
      fname = lambda_fname(phase, elb, vlb)
      if (elb*vlb > 1.0) and cfg.copyarcforperturb:
        fname0 = f"{phase}-e100-v100"
        cmd = f"ln -sf {fname0}.arc {fname}.arc"
        subprocess.run(cmd, shell=True)
        print(GREEN + f" {cmd} " + ENDC)

      xyzfile = fname + ".xyz"
      keyfile = xyzfile.replace("xyz", "key")
      logfile = xyzfile.replace("xyz", "log")
      if (phasedir / logfile).is_file():
        if cfg.verbose > 0:
          print(YELLOW + f" [Warning] {logfile} exists in {phasedir} folder for {fname}" + ENDC)
      else:
        if phase == 'liquid':
          liquidsh = fname + '.sh'
          liquidshs.append(liquidsh)
          with open(liquidsh, 'w') as f:
            f.write(f'source {cfg.tinkerenv}\n')
            if cfg.liquidensemble == "NPT":
              f.write(f'{cfg.phase_dynamic[phase]} {xyzfile} -key {keyfile} {cfg.liquidtotalstep} {cfg.liquidtimestep} {cfg.liquidwriteout} 4 {cfg.liquidT} {cfg.liquidP} > {logfile}\n')
            if cfg.liquidensemble == "NVT":
              f.write(f'{cfg.phase_dynamic[phase]} {xyzfile} -key {keyfile} {cfg.liquidtotalstep} {cfg.liquidtimestep} {cfg.liquidwriteout} 2 {cfg.liquidT} > {logfile}\n')
        if phase == 'gas':
          gassh = fname + '.sh'
          gasshs.append(gassh)
          with open(gassh, 'w') as f:
            f.write(f'source {cfg.tinkerenv}\n')
            f.write(f"{cfg.phase_dynamic[phase]} {xyzfile} -key {keyfile} {cfg.gastotalstep} {cfg.gastimestep} {cfg.gaswriteout} 2 {cfg.gastemperature} > {logfile}\n")

  # submit jobs to clusters
  for phase in cfg.phases:
    phasedir = Path(cfg.homedir) / phase
    os.chdir(phasedir)
    hoststr = os.getenv('HOSTNAME').split('.')[0]
    timestr = str(time.time()).replace('.', '')
    jobpooldir = "/home/liuchw/bin/JobPool/"
    scriptfile = Path(jobpooldir) / f"{hoststr}-{timestr}.sh"
    if (phase == 'gas') and (gasshs != []):
      if os.getlogin() == 'liuchw':  # developer-only
        shstr = f"python /home/liuchw/bin/TinkerGPU2022/submitTinker.py -x {' '.join(gasshs)} -t CPU -n 4 -p {phasedir}"
        print(GREEN + ' ' + shstr + ENDC)
        with open(scriptfile, 'w') as f:
          f.write(shstr)
      else:
        shstr = f"python {cfg.submitexe} -x {' '.join(gasshs)} -t CPU -nodes {' '.join(cfg.nodes)} -n 4 -p {phasedir}"
        print(GREEN + ' ' + shstr + ENDC)
        subprocess.run(shstr, shell=True)
    if (phase == 'liquid') and (liquidshs != []):
      if os.getlogin() == 'liuchw':  # developer-only
        shstr = f"python /home/liuchw/bin/TinkerGPU2022/submitTinker.py -x {' '.join(liquidshs)} -t GPU -p {phasedir}"
        print(GREEN + ' ' + shstr + ENDC)
        with open(scriptfile, 'w') as f:
          f.write(shstr)
      else:
        shstr = f"python {cfg.submitexe} -x {' '.join(liquidshs)} -t GPU -nodes {' '.join(cfg.nodes)} -p {phasedir}"
        print(GREEN + ' ' + shstr + ENDC)
        subprocess.run(shstr, shell=True)
  return


def bar(cfg):
  """Write BAR analysis shell scripts and submit them after MD completes."""
  print(YELLOW + " Checking the completeness of the MD trajectories, please wait... " + ENDC)

  if cfg.copyarcforperturb:
    checkorderparams = cfg.orderparams[:-1]
  else:
    checkorderparams = cfg.orderparams
  if cfg.inputaction == "auto":
    proceed = False
    while not proceed:
      proceed, phase_simsnapshot = checkdynamic(cfg.liquidtotalsnapshot, cfg.gastotalsnapshot, cfg.phases, checkorderparams, cfg.homedir, cfg.verbose)
      if proceed:
        break
      now = datetime.now().strftime("%b %d %Y %H:%M:%S")
      print(YELLOW + f" [{now}] Waiting for dynamic jobs to finish ..." + ENDC)
      time.sleep(cfg.checkingtime)
  else:
    proceed, phase_simsnapshot = checkdynamic(cfg.liquidtotalsnapshot, cfg.gastotalsnapshot, cfg.phases, checkorderparams, cfg.homedir, cfg.verbose)

  if proceed:
    liquidshs = []
    gasshs = []
    for phase in cfg.phases:
      phasedir = Path(cfg.homedir) / phase
      os.chdir(phasedir)
      for i in range(len(cfg.orderparams)-1):
        elb0, vlb0 = cfg.orderparams[i]
        elb1, vlb1 = cfg.orderparams[i+1]
        if elb1 > 1.0:
          elb0, vlb0 = 1.0, 1.0
        elb0s = f"{round(elb0*100):03d}"
        elb1s = f"{round(elb1*100):03d}"
        vlb0s = f"{round(vlb0*100):03d}"
        vlb1s = f"{round(vlb1*100):03d}"
        fname0 = f"{phase}-e{elb0s}-v{vlb0s}"
        fname1 = f"{phase}-e{elb1s}-v{vlb1s}"
        arcfile0 = fname0 + ".arc"
        arcfile1 = fname1 + ".arc"
        liquidbarname = f"bar_e{elb0s}-v{vlb0s}_e{elb1s}-v{vlb1s}.sh"
        gasbarname = f"bar_e{elb1s}-v{vlb1s}_e{elb0s}-v{vlb0s}.sh"
        barfiledir = Path(cfg.homedir) / phase
        # put the FEP files in a separate folder
        if (int(elb1s) > 100) and (int(vlb1s) > 100):
          idx = int(int(elb1s)/10 - 10)
          fepdir = f"FEP_{idx:02d}"
          barfiledir = Path(cfg.homedir) / phase / fepdir
          subprocess.run(f"mkdir -p {cfg.homedir}/{phase}/{fepdir}", shell=True)
          linkarc = f"ln -sf {cfg.homedir}/{phase}/{arcfile0} {cfg.homedir}/{phase}/{fepdir}/{arcfile0}"
          subprocess.run(linkarc, shell=True)
          key0 = f"{phase}-e{elb0s}-v{vlb0s}.key"
          key1 = f"{phase}-e{elb1s}-v{vlb1s}.key"
          linkkey = f"ln -sf {cfg.homedir}/{phase}/{key0} {cfg.homedir}/{phase}/{fepdir}/{key0}"
          subprocess.run(linkkey, shell=True)
          linkkey = f"ln -sf {cfg.homedir}/{phase}/{key1} {cfg.homedir}/{phase}/{fepdir}/{key1}"
          subprocess.run(linkkey, shell=True)
          if cfg.copyarcforperturb:
            linkarc = f"ln -sf {cfg.homedir}/{phase}/{arcfile0} {cfg.homedir}/{phase}/{fepdir}/{arcfile1}"
            subprocess.run(linkarc, shell=True)
          else:
            linkarc = f"ln -sf {cfg.homedir}/{phase}/{arcfile1} {cfg.homedir}/{phase}/{fepdir}/{arcfile1}"
            subprocess.run(linkarc, shell=True)

        if phase == 'liquid':
          outfile = fname0 + ".out"
          barfile = fname0 + ".bar"
          enefile = fname0 + ".ene"
          startsnapshot = int(cfg.liquidtotalsnapshot/5.0) + 1

          if not (barfiledir / barfile).is_file():
            with open(liquidbarname, 'w') as f:
              f.write(f"source {cfg.tinkerenv}\n")
              f.write(f"{cfg.liquidbarexe} 1 {arcfile0} {cfg.liquidT} {arcfile1} {cfg.liquidT} N > {outfile} && \n")
              f.write(f"{cfg.liquidbarexe} 2 {barfile} {startsnapshot} {cfg.liquidtotalsnapshot} 1 {startsnapshot} {cfg.liquidtotalsnapshot} 1 > {enefile} \n")
            if [liquidbarname, str(barfiledir)] not in liquidshs:
              liquidshs.append([liquidbarname, str(barfiledir)])
            if int(elb1s) > 100:
              subprocess.run(f"mv {liquidbarname} {cfg.homedir}/{phase}/{fepdir}", shell=True)
          else:
            if cfg.verbose > 0:
              print(YELLOW + f" [Warning] {barfile} exists in {barfiledir} folder! " + ENDC)

        if phase == 'gas':
          outfile = fname1 + ".out"
          barfile = fname1 + ".bar"
          enefile = fname1 + ".ene"
          startsnapshot = int(cfg.gastotalsnapshot/5.0) + 1

          if cfg.gastotaltime != 0:
            if not (barfiledir / barfile).is_file():
              with open(gasbarname, 'w') as f:
                f.write(f"source {cfg.tinkerenv}\n")
                f.write(f"{cfg.gasbarexe} 1 {arcfile1} {cfg.gastemperature} {arcfile0} {cfg.gastemperature} N > {outfile} && \n")
                f.write(f"{cfg.gasbarexe} 2 {barfile} {startsnapshot} {cfg.gastotalsnapshot} 1 {startsnapshot} {cfg.gastotalsnapshot} 1 > {enefile} \n")
              if [gasbarname, str(barfiledir)] not in gasshs:
                gasshs.append([gasbarname, str(barfiledir)])
              if int(elb1s) > 100:
                subprocess.run(f"mv {gasbarname} {cfg.homedir}/{phase}/{fepdir}", shell=True)
            else:
              if cfg.verbose > 0:
                print(YELLOW + f" [Warning] {barfile} exists in {barfiledir} folder!" + ENDC)
    # submit jobs to clusters
    for phase in cfg.phases:
      phasedir = Path(cfg.homedir) / phase
      os.chdir(phasedir)
      if (phase == 'gas') and (gasshs != []):
        if os.getlogin() == 'liuchw':  # developer-only
          subprocess.run(f"python /home/liuchw/bin/TinkerGPU2022/submitTinker.py -x {' '.join([x[0] for x in gasshs])} -t CPU -n 4 -p {' '.join([x[1] for x in gasshs])}", shell=True)
        else:
          subprocess.run(f"python {cfg.submitexe} -x {' '.join([x[0] for x in gasshs])} -t CPU -n 4 -p {' '.join([x[1] for x in gasshs])} -nodes {' '.join(cfg.nodes)}", shell=True)
      if (phase == 'liquid') and (liquidshs != []):
        if os.getlogin() == 'liuchw':  # developer-only
          subprocess.run(f"python /home/liuchw/bin/TinkerGPU2022/submitTinker.py -x {' '.join([x[0] for x in liquidshs])} -t GPU -p {' '.join([x[1] for x in liquidshs])}", shell=True)
        else:
          subprocess.run(f"python {cfg.submitexe} -x {' '.join([x[0] for x in liquidshs])} -t GPU -p {' '.join([x[1] for x in liquidshs])} -nodes {' '.join(cfg.nodes)}", shell=True)
  return


def result(cfg):
  """Parse BAR .ene files and write the final free energy summary to result.txt."""
  print(YELLOW + " Checking the completeness of the BAR analysis ..." + ENDC)
  proceed, gasperturbsteps, gasenes, liquidperturbsteps, liquidenes, fep_gasperturbsteps, fep_gasenes, fep_liquidperturbsteps, fep_liquidenes = checkbar(cfg.phases, cfg.orderparams, cfg.homedir, cfg.ignoregas, cfg.verbose)

  if proceed:
    if gasenes:
      FEgas, Errgas = map(list, zip(*[parse_ene_file(f) for f in gasenes]))
    else:
      FEgas, Errgas = [], []
    if fep_gasenes:
      fep_FEgas, fep_Errgas = map(list, zip(*[parse_ene_file(f) for f in fep_gasenes]))
    else:
      fep_FEgas, fep_Errgas = [], []
    if liquidenes:
      FEliquid, Errliquid = map(list, zip(*[parse_ene_file(f) for f in liquidenes]))
    else:
      FEliquid, Errliquid = [], []
    if fep_liquidenes:
      fep_FEliquid, fep_Errliquid = map(list, zip(*[parse_ene_file(f) for f in fep_liquidenes]))
    else:
      fep_FEliquid, fep_Errliquid = [], []

    FEall = FEgas + FEliquid
    Errall = Errgas + Errliquid
    totFE = np.array(FEall).sum()
    totErr = np.sqrt(np.square(np.array(Errall)).sum())

    fo = open(Path(cfg.homedir) / "result.txt", "w")
    print(GREEN + f"{'StateA':>20s}{'StateB':>20s}{'FreeEnergy(kcal/mol)':>30s}{'Error(kcal/mol)':>20s}" + ENDC)
    fo.write(f"{'StateA':>20s}{'StateB':>20s}{'FreeEnergy(kcal/mol)':>30s}{'Error(kcal/mol)':>20s}\n")
    for i in range(len(FEgas)-1, -1, -1):
      state0, state1 = gasperturbsteps[i]
      print(f"{state0:>20s}{state1:>20s}{FEgas[i]:>25.4f}{Errgas[i]:>20.4f}")
      fo.write(f"{state0:>20s}{state1:>20s}{FEgas[i]:>25.4f}{Errgas[i]:>20.4f}\n")
    for i in range(len(FEliquid)):
      state0, state1 = liquidperturbsteps[i]
      print(f"{state0:>20s}{state1:>20s}{FEliquid[i]:>25.4f}{Errliquid[i]:>20.4f}")
      fo.write(f"{state0:>20s}{state1:>20s}{FEliquid[i]:>25.4f}{Errliquid[i]:>20.4f}\n")
    print(GREEN + f"{'SUM OF THE TOTAL FREE ENERGY (FE0)':>40s}{totFE:>25.4f}{totErr:>20.4f}" + ENDC)
    fo.write(f"{'SUM OF THE TOTAL FREE ENERGY (FE0)':>40s}{totFE:>25.4f}{totErr:>20.4f}\n")

    print(YELLOW + f"    {'GAS':>20s}{'LIQUID':>20s}{'GAS+LIQUID':>27s}{'GAS+LIQUID+FE0':>20s}" + ENDC)
    fo.write(f"    {'GAS':>20s}{'LIQUID':>20s}{'GAS+LIQUID':>27s}{'GAS+LIQUID+FE0':>20s}\n")
    for i in range(len(fep_FEliquid)):
      feg = 0.0
      if cfg.ignoregas == 0:
        feg = fep_FEgas[i]
      fel = fep_FEliquid[i]
      print(f"    FEP_{i+1:03d}{feg:14.4f}{fel:19.4f}{feg+fel:25.4f}{totFE+feg+fel:20.4f}")
      fo.write(f"    FEP_{i+1:03d}{feg:14.4f}{fel:19.4f}{feg+fel:25.4f}{totFE+feg+fel:20.4f}\n")

    fo.close()
  return


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('act', help="Actions to take.", choices=['setup', 'dynamic', 'bar', 'result', 'auto'], type=str.lower)

  cfg = Config()
  cfg.inputaction = vars(parser.parse_args())['act']
  cfg.rootdir = str(Path(__file__).parent)
  cfg.homedir = str(Path.cwd())

  if not Path(cfg.homedir, "settings.yaml").is_file():
    sys.exit(RED + "Please provide 'settings.yaml' file; " + ENDC + GREEN + f"An example is here for you {Path(cfg.rootdir, 'dat', 'settings.yaml')}" + ENDC)
  else:
    with open('settings.yaml') as f:
      FEsimsettings = yaml.load(f, Loader=yaml.Loader)

  cfg.verbose = int(FEsimsettings.get('verbose', 1))

  lig = FEsimsettings['gas_xyz']
  box = FEsimsettings['box_xyz']
  cfg.prm = FEsimsettings["parameters"]
  cfg.checkingtime = FEsimsettings["checking_time"]
  cfg.natom = int(open(lig).readlines()[0].split()[0])

  cfg.nodes = FEsimsettings.get("node_list", None)
  if cfg.nodes is None:
    if os.getlogin() != 'liuchw':
      sys.exit(RED + "node_list must be provided" + ENDC)

  # developer-only: special node list for liuchw
  if os.getlogin() == 'liuchw':
    cfg.nodes = []
    lines = open('/home/liuchw/bin/TinkerGPU2022/nodes.dat').readlines()
    for line in lines:
      if ("GPU " in line) and ("#" not in line[0]):
        d = line.split()
        if d[1] not in cfg.nodes:
          cfg.nodes.append(d[1])

  if cfg.nodes is None:
    sys.exit(RED + "[Error] node_list must be provided" + ENDC)

  cfg.submitexe = str(Path(cfg.rootdir, "utils", "submitTinker.py"))
  cfg.tinkerenv = str(Path(cfg.rootdir, "dat", "tinker.env"))

  lambdawindow = FEsimsettings["lambda_window"].upper()
  # Accept both spellings: COARSER (correct) and COURSER (legacy)
  if lambdawindow in ("COURSER", "COARSER"):
    orderprmfile = str(Path(cfg.rootdir, "dat", "orderparams_courser"))
  else:
    orderprmfile = str(Path(cfg.rootdir, "dat", "orderparams_default"))
  # user specified order parameter file
  if 'lambda_window_file' in FEsimsettings:
    orderprmfile = FEsimsettings['lambda_window_file']
    print(YELLOW + f" [Warning] You are responsible for your {orderprmfile}" + ENDC)
  else:
    print(GREEN + f" [GOOD] You are using lambda window settings from {orderprmfile}" + ENDC)

  cfg.manualelescale = FEsimsettings.get('manual_ele_scale', False)

  for line in open(orderprmfile).readlines():
    if "#" not in line:
      d = line.split()
      cfg.orderparams.append([float(d[0]), float(d[1])])

  cfg.copyarcforperturb = False
  if "copy_arc_for_perturb" in FEsimsettings:
    cfg.copyarcforperturb = bool(FEsimsettings["copy_arc_for_perturb"])
    print(YELLOW + f" [Warning] Copy ARC file from e100-v100: {bool(cfg.copyarcforperturb)}" + ENDC)

  for i in range(1, 100, 1):  # allow 99 files at a time to do FEP
    if Path(cfg.homedir, f"{cfg.prm}_{i:02d}").is_file():
      fakelambda = round(1.0 + i*0.1, 1)
      cfg.orderparams.append([fakelambda, fakelambda])
      print(GREEN + f" [GOOD] You are doing one-step FEP for parameter file {cfg.prm}_{i:02d}" + ENDC)

  # liquid phase specific settings
  cfg.phases = ['liquid']
  liquidkey = FEsimsettings.get('liquid_key', str(Path(cfg.rootdir, "dat", "liquid.key")))
  if 'liquid_key' in FEsimsettings:
    print(YELLOW + f"[Warning] You are responsible for your {liquidkey}" + ENDC)
  liquidkeylines = open(liquidkey).readlines()

  liquidmdexe = '$DYNAMIC9'
  liquidminexe = '$MINIMIZE9'
  cfg.liquidbarexe = '$BAR9'
  cfg.phase_xyz = {'liquid': box}
  cfg.phase_key = {'liquid': liquidkeylines}

  cfg.liquidtotaltime = FEsimsettings["liquid_md_total_time"]
  cfg.liquidtimestep = FEsimsettings["liquid_md_time_step"]
  cfg.liquidwriteout = FEsimsettings["liquid_md_write_freq"]
  cfg.liquidT = FEsimsettings["liquid_md_temperature"]
  cfg.liquidP = FEsimsettings["liquid_md_pressure"]
  cfg.liquidensemble = FEsimsettings["liquid_md_ensemble"].upper()
  cfg.liquidtotalstep = int((1000000.0*cfg.liquidtotaltime)/cfg.liquidtimestep)
  cfg.liquidtotalsnapshot = int(1000*cfg.liquidtotaltime/cfg.liquidwriteout)
  if not Path('liquid').is_dir():
    subprocess.run("mkdir liquid", shell=True)

  # gas phase specific settings
  gaskey = FEsimsettings.get('gas_key', str(Path(cfg.rootdir, "dat", "gas.key")))
  if 'gas_key' in FEsimsettings:
    print(YELLOW + f"[Warning] You are responsible for your {gaskey}" + ENDC)
  gaskeylines = open(gaskey).readlines()

  gasmdexe = '$DYNAMIC8'
  gasminexe = '$MINIMIZE8'
  cfg.gasbarexe = '$BAR8'

  cfg.gastotaltime = FEsimsettings["gas_md_total_time"]
  cfg.gastimestep = FEsimsettings["gas_md_time_step"]
  cfg.gastemperature = FEsimsettings["gas_md_temperature"]
  cfg.gaswriteout = FEsimsettings["gas_md_write_freq"]
  cfg.gastotalstep = int((1000000.0*cfg.gastotaltime)/cfg.gastimestep)
  cfg.gastotalsnapshot = int(1000*cfg.gastotaltime/cfg.gaswriteout)
  cfg.phase_xyz = {'liquid': box, 'gas': lig}
  cfg.phase_key = {'liquid': liquidkeylines, 'gas': gaskeylines}
  cfg.phase_dynamic = {'liquid': liquidmdexe, 'gas': gasmdexe}
  cfg.phase_minimize = {'liquid': liquidminexe, 'gas': gasminexe}

  # check xyz files
  lines = open(box).readlines()
  natomliquid = int(lines[0].split()[0])
  if natomliquid != len(lines)-2:
    print(YELLOW + f"[Warning] No box info in {box}. Can be in liquid.key instead." + ENDC)
  else:
    [a, b, c] = lines[1].split()[0:3]
    if min([float(a), float(b), float(c)]) < 30.0:
      sys.exit(RED + f"[Error] Please provide a bigger box (>30*30*30)" + ENDC)
  lines = open(lig).readlines()
  natomgas = int(lines[0].split()[0])
  if natomgas < 5:
    cfg.gastotaltime = 0.0
    print(YELLOW + f" [Warning] I set the simulation time to 0 since it only contains {natomgas} atoms" + ENDC)
  if cfg.gastotaltime == 0.0:
    cfg.ignoregas = 1
  else:
    cfg.ignoregas = 0
  if cfg.ignoregas == 0:
    cfg.phases.insert(0, 'gas')
    if not Path('gas').is_dir():
      subprocess.run("mkdir gas", shell=True)

  actions = {'setup': setup, 'dynamic': dynamic, 'bar': bar, 'result': result}
  if cfg.inputaction in actions:
    actions[cfg.inputaction](cfg)
  if cfg.inputaction == 'auto':
    actions['setup'](cfg)
    actions['dynamic'](cfg)

    if cfg.copyarcforperturb:
      checkorderparams = cfg.orderparams[:-1]
    else:
      checkorderparams = cfg.orderparams

    # check if dynamic is complete
    dynamic_good = False
    while not dynamic_good:
      time.sleep(30.0)
      dynamic_good, _ = checkdynamic(cfg.liquidtotalsnapshot, cfg.gastotalsnapshot, cfg.phases, checkorderparams, cfg.homedir, cfg.verbose)

    actions['bar'](cfg)
    # check if bar is complete
    bar_good = False
    while not bar_good:
      time.sleep(30.0)
      bar_good, _, _, _, _, _, _, _, _ = checkbar(cfg.phases, cfg.orderparams, cfg.homedir, cfg.ignoregas, cfg.verbose)

    actions['result'](cfg)
