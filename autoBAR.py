
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import glob
import time
import shutil
import getpass
import argparse
import subprocess
import numpy as np
import ruamel.yaml as yaml
from datetime import datetime
from utils.checkautobar import (
    RED, ENDC, GREEN, YELLOW,
    format_lambda_name, checkdynamic, checkbar, count_arc_snapshots, _read_last_lines,
    _bar_sh_steps_match,
    _bar_file_snapshot_count,
)
from utils.elescale import scaledownele


def _force_symlink(src, dst):
  """Create a symlink, removing any existing file at dst."""
  if os.path.islink(dst) or os.path.exists(dst):
    os.remove(dst)
  os.symlink(src, dst)


def _read_free_energy(enefile):
  """Read free energy and error from a BAR .ene file.

  Prefers 'BAR Iteration' over 'BAR Bootstrap'.
  Returns (free_energy, error) or None if not found.
  """
  with open(enefile) as f:
    for line in f:
      if "Free Energy via BAR Iteration" in line:
        tokens = line.split()
        return float(tokens[-4]), float(tokens[-2])
    # rewind and try bootstrap
    f.seek(0)
    for line in f:
      if "Free Energy via BAR Bootstrap" in line:
        tokens = line.split()
        return float(tokens[-4]), float(tokens[-2])
  return None


def _remove_phase_files(phase_dir, extensions):
  """Remove files matching given extensions from a directory."""
  for ext in extensions:
    for fpath in glob.glob(os.path.join(phase_dir, f"*.{ext}")):
      os.remove(fpath)


def _ene_complete(enepath):
  """Return True if enepath exists and contains the BAR convergence line."""
  if not os.path.isfile(enepath):
    return False
  return any("BAR Estimate of -T*dS" in line for line in _read_last_lines(enepath))



def _remaining_steps(arcpath, total_steps, steps_per_snapshot):
  """Return (remaining_steps, existing_snapshots) given an existing arc file.

  If the arc is absent or unreadable, returns (total_steps, 0).
  If the arc already meets or exceeds the target, returns (0, existing).
  Tinker restarts from the .dyn checkpoint and appends to the arc, so
  we only need to request the outstanding steps.
  """
  if not os.path.isfile(arcpath):
    return total_steps, 0
  existing = count_arc_snapshots(arcpath)
  if not isinstance(existing, int):   # count_arc_snapshots returns str on error
    return total_steps, 0
  remaining = total_steps - existing * steps_per_snapshot
  return max(0, remaining), existing


def setup():
  for phase in phases:
    xyzfile = phase_xyz[phase]
    keylines = phase_key[phase]
    with open(f"{phase}.key", 'w') as f:
      for line in keylines:
        if 'parameters' in line.lower():
          line = f'parameters     {prm}\n'
        f.write(line)

    if phase == 'gas':
      gasminsh = 'gas-min.sh'
      with open(gasminsh, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write(f'source {tinkerenv}\n')
        f.write(f'{phase_minimize[phase]} {xyzfile} -key gas.key 0.1 > gas-min.log \n')
        f.write(f'wait\nmv {phase_xyz[phase]}_2 {phase_xyz[phase]}\n')
      if not os.path.isfile("gas-min.log"):
        rc = subprocess.run(['bash', gasminsh]).returncode
        if rc != 0 or not os.path.isfile(phase_xyz[phase]):
          sys.exit(RED + f"[Error] gas minimization failed (rc={rc}); see gas-min.log" + ENDC)
    elif phase == 'liquid':
      liquidminsh = 'liquid-min.sh'
      with open(liquidminsh, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write(f'source {tinkerenv}\n')
        f.write(f'{phase_minimize[phase]} {xyzfile} -key liquid.key 0.2 > liquid-min.log\n')
        f.write(f'wait\nmv {phase_xyz[phase]}_2 {phase_xyz[phase]}\n')
      if not os.path.isfile("liquid-min.log"):
        rc = subprocess.run(['bash', liquidminsh]).returncode
        if rc != 0 or not os.path.isfile(phase_xyz[phase]):
          sys.exit(RED + f"[Error] liquid minimization failed (rc={rc}); see liquid-min.log" + ENDC)

    _remove_phase_files(phase, ['xyz', 'key'])

    if not (ignoregas == 1 and phase == 'gas'):
      prmfiles = []
      for elb, vlb in orderparams:
        fname = format_lambda_name(phase, elb, vlb)
        with open(fname + ".key", 'w') as fw:
          for line in keylines:
            if 'parameters' in line.lower():
              if (elb * vlb > 1.0):
                assert elb == vlb, RED + f" Error: lambdas greater than 1 but not the same: {elb}, {vlb} " + ENDC
                idx = int(round(elb * 100) / 10) - 10
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
          eff_elb = 1.0 if (elb * vlb > 1.0) else elb
          eff_vlb = 1.0 if (elb * vlb > 1.0) else vlb

          if manualelescale:
            fw.write('ele-lambda 1.0\n')
          else:
            fw.write(f'ele-lambda {eff_elb}\n')
          fw.write(f'vdw-lambda {eff_vlb}\n')

          if manualelescale:
            scaledprms = scaledownele(phase_xyz['gas'], prm, eff_elb)
            fw.write(f'\n# scale down electrostatic related parameters by {eff_elb}\n')
            for s in scaledprms:
              fw.write(f'{s}\n')

        xyz_src = os.path.join(homedir, phase_xyz[phase])
        xyz_dst = os.path.join(homedir, phase, f"{fname}.xyz")
        _force_symlink(xyz_src, xyz_dst)
        shutil.move(fname + ".key", os.path.join(phase, fname + ".key"))

      # make copy of prmfile
      for prmfile in prmfiles:
        if os.path.isfile(prmfile):
          shutil.copy2(prmfile, os.path.join(phase, prmfile))
  print(GREEN + ' [GOOD] BAR simulation files generated!' + ENDC)
  return


def dynamic():
  liquidshs = []
  gasshs = []
  # steps_per_snapshot: how many MD steps produce one frame written to the arc
  liquid_steps_per_snap = int(round(liquidwriteout / liquidtimestep * 1000))
  gas_steps_per_snap    = int(round(gaswriteout    / gastimestep    * 1000))

  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    for shfile in glob.glob(os.path.join(phasedir, "*.sh")):
      if not os.path.basename(shfile).startswith('bar_'):
        os.remove(shfile)
    os.chdir(phasedir)
    for elb, vlb in orderparams:
      fname = format_lambda_name(phase, elb, vlb)
      if (elb * vlb > 1.0) and copyarcforperturb:
        fname0 = f"{phase}-e100-v100"
        src_arc = os.path.join(phasedir, f"{fname0}.arc")
        dst_arc = os.path.join(phasedir, f"{fname}.arc")
        _force_symlink(src_arc, dst_arc)
        print(GREEN + f" ln -sf {fname0}.arc {fname}.arc " + ENDC)
        continue

      xyzfile = fname + ".xyz"
      keyfile = fname + ".key"
      logfile = fname + ".log"
      arcpath = os.path.join(phasedir, fname + ".arc")

      if phase == 'liquid':
        remaining, existing_snaps = _remaining_steps(arcpath, liquidtotalstep, liquid_steps_per_snap)
        total_snap = liquidtotalsnapshot
      else:
        remaining, existing_snaps = _remaining_steps(arcpath, gastotalstep, gas_steps_per_snap)
        total_snap = gastotalsnapshot

      if remaining <= 0:
        if verbose > 0:
          print(GREEN + f" [Done]  {fname}: arc already has {total_snap} snapshots, skipping" + ENDC)
        continue

      dynpath = os.path.join(phasedir, fname + ".dyn")
      errpath = os.path.join(phasedir, fname + ".err")
      is_continuation = existing_snaps > 0 and os.path.isfile(dynpath)
      if is_continuation and os.path.isfile(errpath):
        print(RED + f" [Error]  {fname}: .err file detected, skipping resume."
              f" Please check {errpath} before resubmitting." + ENDC)
        continue
      log_redirect = ">>" if is_continuation else ">"
      if is_continuation and verbose > 0:
        print(YELLOW + f" [Resume] {fname}: {existing_snaps}/{total_snap} snapshots done,"
              f" running {remaining} more steps" + ENDC)
      elif existing_snaps > 0 and not os.path.isfile(dynpath) and verbose > 0:
        print(YELLOW + f" [Restart] {fname}: arc has {existing_snaps} snapshots but .dyn missing,"
              f" restarting from scratch" + ENDC)

      if phase == 'liquid':
        liquidsh = fname + '.sh'
        liquidshs.append(liquidsh)
        with open(liquidsh, 'w') as f:
          f.write(f'source {tinkerenv}\n')
          if liquidensemble == "NPT":
            f.write(f'{phase_dynamic[phase]} {xyzfile} -key {keyfile} {remaining} {liquidtimestep} {liquidwriteout} 4 {liquidT} {liquidP} {log_redirect} {logfile}\n')
          elif liquidensemble == "NVT":
            f.write(f'{phase_dynamic[phase]} {xyzfile} -key {keyfile} {remaining} {liquidtimestep} {liquidwriteout} 2 {liquidT} {log_redirect} {logfile}\n')
      elif phase == 'gas':
        gassh = fname + '.sh'
        gasshs.append(gassh)
        with open(gassh, 'w') as f:
          f.write(f'source {tinkerenv}\n')
          f.write(f"{phase_dynamic[phase]} {xyzfile} -key {keyfile} {remaining} {gastimestep} {gaswriteout} 2 {gastemperature} {log_redirect} {logfile}\n")

  # submit jobs to clusters
  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    os.chdir(phasedir)
    if phase == 'gas' and gasshs:
      nodes_arg = f" -nodes {' '.join(nodes)}" if nodes else ""
      shstr = f"python {submitexe} -x {' '.join(gasshs)} -t CPU{nodes_arg} -n 4 -p {phasedir}"
      print(GREEN + ' ' + shstr + ENDC)
      rc = subprocess.run(shstr, shell=True).returncode
      if rc != 0:
        print(RED + f"[Warning] gas submitTinker exited with code {rc}; some jobs may not have been dispatched" + ENDC)
    if phase == 'liquid' and liquidshs:
      nodes_arg = f" -nodes {' '.join(nodes)}" if nodes else ""
      shstr = f"python {submitexe} -x {' '.join(liquidshs)} -t GPU{nodes_arg} -p {phasedir}"
      print(GREEN + ' ' + shstr + ENDC)
      rc = subprocess.run(shstr, shell=True).returncode
      if rc != 0:
        print(RED + f"[Warning] liquid submitTinker exited with code {rc}; some jobs may not have been dispatched" + ENDC)
  return


def bar():
  print(YELLOW + " Checking the completeness of the MD trajectories, please wait... " + ENDC)

  checkorderparams = orderparams[:-1] if copyarcforperturb else orderparams
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

  if not proceed:
    return

  liquidshs = []
  gasshs = []
  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    os.chdir(phasedir)
    for i in range(len(orderparams) - 1):
      elb0, vlb0 = orderparams[i]
      elb1, vlb1 = orderparams[i + 1]
      if elb1 > 1.0:
        elb0, vlb0 = 1.0, 1.0
      e0 = f"{round(elb0 * 100):03d}"
      e1 = f"{round(elb1 * 100):03d}"
      v0 = f"{round(vlb0 * 100):03d}"
      v1 = f"{round(vlb1 * 100):03d}"
      fname0 = f"{phase}-e{e0}-v{v0}"
      fname1 = f"{phase}-e{e1}-v{v1}"
      arcfile0 = fname0 + ".arc"
      arcfile1 = fname1 + ".arc"
      liquidbarname = f"bar_e{e0}-v{v0}_e{e1}-v{v1}.sh"
      gasbarname = f"bar_e{e1}-v{v1}_e{e0}-v{v0}.sh"
      barfiledir = os.path.join(homedir, phase)
      # put the FEP files in a separate folder
      if int(e1) > 100 and int(v1) > 100:
        idx = int(int(e1) / 10 - 10)
        fepdir = f"FEP_{idx:02d}"
        barfiledir = os.path.join(homedir, phase, fepdir)
        os.makedirs(barfiledir, exist_ok=True)
        _force_symlink(
          os.path.join(homedir, phase, arcfile0),
          os.path.join(barfiledir, arcfile0)
        )
        for keyname in [f"{phase}-e{e0}-v{v0}.key", f"{phase}-e{e1}-v{v1}.key"]:
          _force_symlink(
            os.path.join(homedir, phase, keyname),
            os.path.join(barfiledir, keyname)
          )
        arc_src = arcfile0 if copyarcforperturb else arcfile1
        _force_symlink(
          os.path.join(homedir, phase, arc_src),
          os.path.join(barfiledir, arcfile1)
        )

      if phase == 'liquid':
        outfile = fname0 + ".out"
        barfile = fname0 + ".bar"
        enefile = fname0 + ".ene"
        startsnapshot = int(liquidtotalsnapshot / 5.0) + 1

        barpath = os.path.join(barfiledir, barfile)
        enepath = os.path.join(barfiledir, enefile)
        shpath  = os.path.join(barfiledir, liquidbarname)
        if os.path.isfile(barpath) and os.path.isfile(shpath) and not _bar_sh_steps_match(shpath, startsnapshot, liquidtotalsnapshot):
          print(YELLOW + f" [Regen] {barfile}: snapshot range changed, removing stale .bar/.ene" + ENDC)
          for p in [barpath, enepath]:
            if os.path.isfile(p): os.remove(p)

        if os.path.isfile(barpath) and _bar_file_snapshot_count(barpath) != liquidtotalsnapshot:
          print(YELLOW + f" [Regen] {barfile}: .bar is from a different run ({_bar_file_snapshot_count(barpath)} vs {liquidtotalsnapshot} snapshots), removing stale .bar/.ene" + ENDC)
          for p in [barpath, enepath]:
            if os.path.isfile(p): os.remove(p)

        bar_done = os.path.isfile(barpath) and _ene_complete(enepath) and os.path.isfile(shpath)
        if not bar_done:
          if os.path.isfile(barpath) and verbose > 0:
            print(YELLOW + f" [Rerun] {barfile}: .bar exists but .ene missing/incomplete, resubmitting" + ENDC)
          with open(liquidbarname, 'w') as f:
            f.write(f"source {tinkerenv}\n")
            f.write(f"{liquidbarexe} 1 {arcfile0} {liquidT} {arcfile1} {liquidT} N > {outfile} && \n")
            f.write(f"{liquidbarexe} 2 {barfile} {startsnapshot} {liquidtotalsnapshot} 1 {startsnapshot} {liquidtotalsnapshot} 1 > {enefile} \n")
          if [liquidbarname, barfiledir] not in liquidshs:
            liquidshs.append([liquidbarname, barfiledir])
          if int(e1) > 100:
            shutil.move(liquidbarname, os.path.join(barfiledir, liquidbarname))
        else:
          if verbose > 0:
            print(YELLOW + f" [Warning] {barfile} exists in {barfiledir} folder! " + ENDC)

      elif phase == 'gas':
        outfile = fname1 + ".out"
        barfile = fname1 + ".bar"
        enefile = fname1 + ".ene"
        startsnapshot = int(gastotalsnapshot / 5.0) + 1

        if gastotaltime != 0:
          barpath = os.path.join(barfiledir, barfile)
          enepath = os.path.join(barfiledir, enefile)
          shpath  = os.path.join(barfiledir, gasbarname)
          if os.path.isfile(barpath) and os.path.isfile(shpath) and not _bar_sh_steps_match(shpath, startsnapshot, gastotalsnapshot):
            print(YELLOW + f" [Regen] {barfile}: snapshot range changed, removing stale .bar/.ene" + ENDC)
            for p in [barpath, enepath]:
              if os.path.isfile(p): os.remove(p)

          if os.path.isfile(barpath) and _bar_file_snapshot_count(barpath) != gastotalsnapshot:
            print(YELLOW + f" [Regen] {barfile}: .bar is from a different run ({_bar_file_snapshot_count(barpath)} vs {gastotalsnapshot} snapshots), removing stale .bar/.ene" + ENDC)
            for p in [barpath, enepath]:
              if os.path.isfile(p): os.remove(p)

          bar_done = os.path.isfile(barpath) and _ene_complete(enepath) and os.path.isfile(shpath)
          if not bar_done:
            if os.path.isfile(barpath) and verbose > 0:
              print(YELLOW + f" [Rerun] {barfile}: .bar exists but .ene missing/incomplete, resubmitting" + ENDC)
            with open(gasbarname, 'w') as f:
              f.write(f"source {tinkerenv}\n")
              f.write(f"{gasbarexe} 1 {arcfile1} {gastemperature} {arcfile0} {gastemperature} N > {outfile} && \n")
              f.write(f"{gasbarexe} 2 {barfile} {startsnapshot} {gastotalsnapshot} 1 {startsnapshot} {gastotalsnapshot} 1 > {enefile} \n")
            if [gasbarname, barfiledir] not in gasshs:
              gasshs.append([gasbarname, barfiledir])
            if int(e1) > 100:
              shutil.move(gasbarname, os.path.join(barfiledir, gasbarname))
          else:
            if verbose > 0:
              print(YELLOW + f" [Warning] {barfile} exists in {barfiledir} folder!" + ENDC)

  # submit jobs to clusters
  for phase in phases:
    phasedir = os.path.join(homedir, phase)
    os.chdir(phasedir)
    if phase == 'gas' and gasshs:
      nodes_arg = f" -nodes {' '.join(nodes)}" if nodes else ""
      cmd = f"python {submitexe} -x {' '.join([x[0] for x in gasshs])} -t CPU -n 4 -p {' '.join([x[1] for x in gasshs])}{nodes_arg}"
      rc = subprocess.run(cmd, shell=True).returncode
      if rc != 0:
        print(RED + f"[Warning] gas BAR submitTinker exited with code {rc}; some jobs may not have been dispatched" + ENDC)
    if phase == 'liquid' and liquidshs:
      nodes_arg = f" -nodes {' '.join(nodes)}" if nodes else ""
      cmd = f"python {submitexe} -x {' '.join([x[0] for x in liquidshs])} -t GPU -p {' '.join([x[1] for x in liquidshs])}{nodes_arg}"
      rc = subprocess.run(cmd, shell=True).returncode
      if rc != 0:
        print(RED + f"[Warning] liquid BAR submitTinker exited with code {rc}; some jobs may not have been dispatched" + ENDC)
  return


def result():
  print(YELLOW + " Checking the completeness of the BAR analysis ..." + ENDC)
  proceed, gasperturbsteps, gasenes, liquidperturbsteps, liquidenes, fep_gasperturbsteps, fep_gasenes, fep_liquidperturbsteps, fep_liquidenes = checkbar(phases, orderparams, homedir, ignoregas, verbose, liquidtotalsnapshot, gastotalsnapshot)

  if not proceed:
    return

  FEgas, Errgas = [], []
  FEliquid, Errliquid = [], []
  fep_FEgas, fep_Errgas = [], []
  fep_FEliquid, fep_Errliquid = [], []

  for enefile, fe_list, err_list in [
      (gasenes, FEgas, Errgas),
      (fep_gasenes, fep_FEgas, fep_Errgas),
      (liquidenes, FEliquid, Errliquid),
      (fep_liquidenes, fep_FEliquid, fep_Errliquid),
  ]:
    for ene in enefile:
      result_pair = _read_free_energy(ene)
      if result_pair is not None:
        fe_list.append(result_pair[0])
        err_list.append(result_pair[1])

  FEall = FEgas + FEliquid
  Errall = Errgas + Errliquid
  totFE = np.array(FEall).sum()
  totErr = np.sqrt(np.square(np.array(Errall)).sum())

  with open(os.path.join(homedir, "result.txt"), "w") as fo:
    header = "%20s%20s%30s%20s" % ("StateA", "StateB", "FreeEnergy(kcal/mol)", "Error(kcal/mol)")
    print(GREEN + header + ENDC)
    fo.write(header + "\n")
    for i in range(len(FEgas) - 1, -1, -1):
      state0, state1 = gasperturbsteps[i]
      line = "%20s%20s%25.4f%20.4f" % (state0, state1, FEgas[i], Errgas[i])
      print(line)
      fo.write(line + "\n")
    for i in range(len(FEliquid)):
      state0, state1 = liquidperturbsteps[i]
      line = "%20s%20s%25.4f%20.4f" % (state0, state1, FEliquid[i], Errliquid[i])
      print(line)
      fo.write(line + "\n")
    sumline = "%40s%25.4f%20.4f" % ("SUM OF THE TOTAL FREE ENERGY (FE0)", totFE, totErr)
    print(GREEN + sumline + ENDC)
    fo.write(sumline + "\n")

    fep_header = "    %20s%20s%27s%20s" % ("GAS", "LIQUID", "GAS+LIQUID", "GAS+LIQUID+FE0")
    print(YELLOW + fep_header + ENDC)
    fo.write(fep_header + "\n")
    for i in range(len(fep_FEliquid)):
      feg = fep_FEgas[i] if ignoregas == 0 else 0.0
      fel = fep_FEliquid[i]
      line = f"    FEP_{i+1:03d}{feg:14.4f}{fel:19.4f}{feg+fel:25.4f}{totFE+feg+fel:20.4f}"
      print(line)
      fo.write(line + "\n")
  return


def _load_settings():
  """Load and validate settings from settings.yaml. Returns config dict."""
  rootdir = os.path.join(os.path.split(__file__)[0])
  homedir = os.getcwd()

  settings_path = os.path.join(homedir, "settings.yaml")
  if not os.path.isfile(settings_path):
    sys.exit(RED + "Please provide 'settings.yaml' file; " + ENDC +
             GREEN + f"An example is here for you {os.path.join(rootdir, 'dat', 'settings.yaml')}" + ENDC)

  # Use the object-style API so this works on ruamel.yaml >=0.18,
  # which removed the module-level yaml.load(..., Loader=...) shim.
  yaml_parser = yaml.YAML(typ='safe', pure=True)
  with open(settings_path) as f:
    settings = yaml_parser.load(f)
  return rootdir, homedir, settings


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('act', help="Actions to take.", choices=['setup', 'dynamic', 'bar', 'result', 'auto'], type=str.lower)

  # global settings
  global inputaction
  inputaction = vars(parser.parse_args())['act']
  global rootdir, homedir
  rootdir, homedir, FEsimsettings = _load_settings()

  global verbose
  verbose = int(FEsimsettings.get('verbose', 1))

  global prm, checkingtime
  lig = FEsimsettings['gas_xyz']
  box = FEsimsettings['box_xyz']
  prm = FEsimsettings["parameters"]
  checkingtime = FEsimsettings["checking_time"]
  global natom
  with open(lig) as f:
    natom = int(f.readline().split()[0])
  global nodes
  nodes = FEsimsettings.get("node_list")
  if nodes is None:
    if getpass.getuser() != "liuchw":
      sys.exit(RED + "[Error] node_list must be provided in settings.yaml" + ENDC)
    else:
      nodes = []

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
  if 'lambda_window_file' in FEsimsettings:
    orderprmfile = FEsimsettings['lambda_window_file']
    print(YELLOW + f" [Warning] You are responsible for your {orderprmfile}" + ENDC)
  else:
    print(GREEN + f" [GOOD] You are using lambda window settings from {orderprmfile}" + ENDC)

  global manualelescale
  manualelescale = FEsimsettings.get('manual_ele_scale', False)

  with open(orderprmfile) as f:
    for line in f:
      if "#" not in line:
        d = line.split()
        if len(d) >= 2:
          orderparams.append([float(d[0]), float(d[1])])

  global copyarcforperturb
  copyarcforperturb = False
  if "copy_arc_for_perturb" in FEsimsettings:
    copyarcforperturb = bool(FEsimsettings["copy_arc_for_perturb"])
    print(YELLOW + f" [Warning] Copy ARC file from e100-v100: {copyarcforperturb}" + ENDC)

  for i in range(1, 100):  # allow 99 files at a time to do FEP
    if os.path.isfile(os.path.join(homedir, f"{prm}_{i:02d}")):
      fakelambda = round(1.0 + i * 0.1, 1)
      orderparams.append([fakelambda, fakelambda])
      print(GREEN + f" [GOOD] You are doing one-step FEP for parameter file {prm}_{i:02d}" + ENDC)

  # liquid phase specific settings
  global phases, liquidkeylines
  phases = ['liquid']
  _lk = FEsimsettings.get('liquid_key')
  if _lk and os.path.isfile(_lk):
    liquidkey = _lk
    print(YELLOW + f"[Warning] You are responsible for your {liquidkey}" + ENDC)
  else:
    liquidkey = os.path.join(rootdir, "dat", "liquid.key")
  with open(liquidkey) as f:
    liquidkeylines = f.readlines()

  global liquidmdexe, liquidbarexe
  liquidmdexe = '$DYNAMIC9'
  liquidminexe = '$MINIMIZE9'
  liquidbarexe = '$BAR9'
  global phase_xyz, phase_key, phase_dynamic, phase_minimize
  phase_xyz = {'liquid': box}
  phase_key = {'liquid': liquidkeylines}
  global liquidtotaltime, liquidtimestep, liquidwriteout, liquidtotalstep, liquidtotalsnapshot
  global liquidT, liquidP, liquidensemble
  liquidtotaltime = FEsimsettings["liquid_md_total_time"]
  liquidtimestep = FEsimsettings["liquid_md_time_step"]
  liquidwriteout = FEsimsettings["liquid_md_write_freq"]
  liquidT = FEsimsettings["liquid_md_temperature"]
  liquidP = FEsimsettings["liquid_md_pressure"]
  liquidensemble = FEsimsettings["liquid_md_ensemble"].upper()
  liquidtotalstep = int((1000000.0 * liquidtotaltime) / liquidtimestep)
  liquidtotalsnapshot = int(1000 * liquidtotaltime / liquidwriteout)
  # BAR skips the first 20 % as equilibration (startsnapshot = total/5 + 1);
  # need at least 6 snapshots so the sampled window is non-empty.
  if liquidtotalsnapshot < 6:
    sys.exit(RED + f"[Error] liquid_md_total_time ({liquidtotaltime} ns) / "
             f"liquid_md_write_freq ({liquidwriteout} ps) yields only "
             f"{liquidtotalsnapshot} snapshots; need at least 6 for BAR" + ENDC)
  os.makedirs('liquid', exist_ok=True)

  # gas phase specific settings
  global ignoregas
  global gaskeylines
  if 'gas_key' in FEsimsettings:
    gaskey = FEsimsettings['gas_key']
    print(YELLOW + f"[Warning] You are responsible for your {gaskey}" + ENDC)
  else:
    gaskey = os.path.join(rootdir, "dat", "gas.key")

  with open(gaskey) as f:
    gaskeylines = f.readlines()

  global gasmdexe, gasbarexe
  gasmdexe = '$DYNAMIC8'
  gasminexe = '$MINIMIZE8'
  gasbarexe = '$BAR8'
  global gastotaltime, gastimestep, gaswriteout, gastemperature, gastotalstep, gastotalsnapshot
  gastotaltime = FEsimsettings["gas_md_total_time"]
  gastimestep = FEsimsettings["gas_md_time_step"]
  gastemperature = FEsimsettings["gas_md_temperature"]
  gaswriteout = FEsimsettings["gas_md_write_freq"]
  gastotalstep = int((1000000.0 * gastotaltime) / gastimestep)
  gastotalsnapshot = int(1000 * gastotaltime / gaswriteout)
  # Gas phase is optional (gastotaltime == 0 disables it), but if enabled
  # it must also yield enough snapshots for BAR to have a non-empty window.
  if gastotaltime > 0 and gastotalsnapshot < 6:
    sys.exit(RED + f"[Error] gas_md_total_time ({gastotaltime} ns) / "
             f"gas_md_write_freq ({gaswriteout} ps) yields only "
             f"{gastotalsnapshot} snapshots; need at least 6 for BAR" + ENDC)
  phase_xyz = {'liquid': box, 'gas': lig}
  phase_key = {'liquid': liquidkeylines, 'gas': gaskeylines}
  phase_dynamic = {'liquid': liquidmdexe, 'gas': gasmdexe}
  phase_minimize = {'liquid': liquidminexe, 'gas': gasminexe}

  # check xyz files
  with open(box) as f:
    boxlines = f.readlines()
  natomliquid = int(boxlines[0].split()[0])
  if natomliquid != len(boxlines) - 2:
    print(YELLOW + f"[Warning] No box info in {box}. Can be in liquid.key instead." + ENDC)
  else:
    box_tokens = boxlines[1].split()
    if len(box_tokens) < 3:
      sys.exit(RED + f"[Error] Box file {box} line 2 must contain at least 3 values for box dimensions (a, b, c)" + ENDC)
    a, b, c = float(box_tokens[0]), float(box_tokens[1]), float(box_tokens[2])
    if min(a, b, c) < 30.0:
      sys.exit(RED + "[Error] Please provide a bigger box (>30*30*30)" + ENDC)

  with open(lig) as f:
    liglines = f.readlines()
  natomgas = int(liglines[0].split()[0])
  if natomgas < 5:
    gastotaltime = 0.0
    print(YELLOW + f" [Warning] I set the simulation time to 0 since it only contains {natomgas} atoms" + ENDC)
  ignoregas = 1 if gastotaltime == 0.0 else 0
  if ignoregas == 0:
    phases.insert(0, 'gas')
    os.makedirs('gas', exist_ok=True)

  actions = {'setup': setup, 'dynamic': dynamic, 'bar': bar, 'result': result}
  if inputaction in actions:
    actions[inputaction]()
  if inputaction == 'auto':
    actions['setup']()
    actions['dynamic']()

    checkorderparams = orderparams[:-1] if copyarcforperturb else orderparams

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
      bar_good, *_ = checkbar(phases, orderparams, homedir, ignoregas, verbose, liquidtotalsnapshot, gastotalsnapshot)

    actions['result']()
