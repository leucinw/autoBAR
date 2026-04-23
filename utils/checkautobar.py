
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import subprocess
from collections import deque

RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\33[93m'

def format_lambda_name(phase, elb, vlb):
    """Format a lambda window name from phase and lambda values."""
    return f"{phase}-e{round(elb * 100):03d}-v{round(vlb * 100):03d}"

def _read_last_lines(filepath, n=5):
    """Read the last n lines of a file efficiently using a fixed-size deque."""
    with open(filepath) as f:
        return deque(f, maxlen=n)

def _append_unique_ene(enedir, enefile, ene_list, perturb_list, perturb_pair):
    """Append an .ene path and its perturbation pair if not already present."""
    path = os.path.join(enedir, enefile)
    if path not in ene_list:
        ene_list.append(path)
        perturb_list.append(perturb_pair)

def _bar_sh_steps_match(shpath, expected_start, expected_total):
    """Return True if the BAR step-2 line in shpath uses the expected snapshot range.

    The step-2 line written by autoBAR looks like:
      $BAR? 2 barfile <start> <total> 1 <start> <total> 1 > enefile
    so parts[1]=='2', parts[3]==start, parts[4]==total.
    Returns False if the file is absent or cannot be parsed.
    """
    if not os.path.isfile(shpath):
        return False
    with open(shpath) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 5 and parts[1] == '2':
                try:
                    return int(parts[3]) == expected_start and int(parts[4]) == expected_total
                except (ValueError, IndexError):
                    pass
    return False

def count_arc_snapshots(file_path):
    if not os.path.exists(file_path):
        return "File not found."

    with open(file_path, 'rb') as f:
        first_line = f.readline().decode().split()
        if not first_line:
            return 0
        n_atoms = int(first_line[0])

        second_line = f.readline().decode().split()
        if not second_line:
            return 1
        stride = (n_atoms + 1) if second_line[0] == "1" else (n_atoms + 2)

    # Use wc -l which is implemented in C and vastly faster than Python loops
    result = subprocess.run(['wc', '-l', file_path],
                            capture_output=True, text=True)
    total_lines = int(result.stdout.split()[0])

    # Handle missing trailing newline
    with open(file_path, 'rb') as f:
        f.seek(-1, os.SEEK_END)
        if f.read(1) != b'\n':
            total_lines += 1

    return total_lines // stride

def checkdynamic(liquidtotalsnapshot, gastotalsnapshot, phases, orderparams, homedir, verbose):
  statuslist = []
  phase_simsnapshot = {'liquid': liquidtotalsnapshot, 'gas': gastotalsnapshot}
  for phase in phases:
    for elb, vlb in orderparams:
      fname = format_lambda_name(phase, elb, vlb)
      arcfile = fname + ".arc"
      arcpath = os.path.join(homedir, phase, arcfile)
      if os.path.isfile(arcpath):
        simsnapshot = count_arc_snapshots(arcpath)
        per = int(simsnapshot / phase_simsnapshot[phase] * 100)
        if simsnapshot >= phase_simsnapshot[phase]:
          if verbose > 0:
            print(GREEN + f"{fname:>20s}: " + u'\u2584' * (int(per / 2)) + f" [{per:>3d}%]" + ENDC)
          statuslist.append(True)
        else:
          if verbose > 0:
            print(YELLOW + f"{fname:>20s}: " + u'\u2584' * (int(per / 2)) + f" [{per:>3d}%]" + ENDC)
          statuslist.append(False)
      else:
        # Only treat a missing arc as "OK" for the gas phase when gas is
        # disabled (gas_md_total_time = 0). A missing liquid arc must
        # always fail, regardless of the gas setting.
        if phase == 'gas' and gastotalsnapshot == 0:
          statuslist.append(True)
        else:
          if verbose > 0:
            print(RED + f"{arcfile} does not exist" + ENDC)
          statuslist.append(False)
  completed = all(statuslist)
  return completed, phase_simsnapshot

def checkbar(phases, orderparams, homedir, ignoregas, verbose, liquidtotalsnapshot=0, gastotalsnapshot=0):
  statuslist = []
  gasenes = []
  liquidenes = []
  gasperturbsteps = []
  liquidperturbsteps = []
  fep_gasenes = []
  fep_liquidenes = []
  fep_gasperturbsteps = []
  fep_liquidperturbsteps = []
  ene_to_sh = {}  # maps ene_path → expected sh_path for the current pair
  checkphases = phases
  if ignoregas == 1:
    checkphases = ['liquid']

  for phase in checkphases:
    for (elb0, vlb0), (elb1, vlb1) in zip(orderparams, orderparams[1:]):
      if elb1 > 1.0:
        elb0, vlb0 = 1.0, 1.0

      fname0 = format_lambda_name(phase, elb0, vlb0)
      fname1 = format_lambda_name(phase, elb1, vlb1)

      elb1_int = round(elb1 * 100)
      vlb1_int = round(vlb1 * 100)
      is_fep = elb1_int > 100 and vlb1_int > 100

      enedir = os.path.join(homedir, phase)
      if is_fep:
        idx = int(elb1_int / 10 - 10)
        enedir = os.path.join(homedir, phase, f'FEP_{idx:02d}')

      e0 = f"{round(elb0 * 100):03d}"
      e1 = f"{round(elb1 * 100):03d}"
      v0 = f"{round(vlb0 * 100):03d}"
      v1 = f"{round(vlb1 * 100):03d}"

      if phase == 'gas':
        enefile = fname1 + ".ene"
        enepath = os.path.join(enedir, enefile)
        _gas_start = int(gastotalsnapshot / 5.0) + 1
        ene_to_sh[enepath] = (os.path.join(enedir, f"bar_e{e1}-v{v1}_e{e0}-v{v0}.sh"), _gas_start, gastotalsnapshot)
        if is_fep:
          _append_unique_ene(enedir, enefile, fep_gasenes, fep_gasperturbsteps, [fname1, fname0])
        else:
          _append_unique_ene(enedir, enefile, gasenes, gasperturbsteps, [fname1, fname0])

      if phase == 'liquid':
        enefile = fname0 + ".ene"
        enepath = os.path.join(enedir, enefile)
        _liq_start = int(liquidtotalsnapshot / 5.0) + 1
        ene_to_sh[enepath] = (os.path.join(enedir, f"bar_e{e0}-v{v0}_e{e1}-v{v1}.sh"), _liq_start, liquidtotalsnapshot)
        if is_fep:
          _append_unique_ene(enedir, enefile, fep_liquidenes, fep_liquidperturbsteps, [fname0, fname1])
        else:
          _append_unique_ene(enedir, enefile, liquidenes, liquidperturbsteps, [fname0, fname1])

  for enefile in gasenes + liquidenes + fep_gasenes + fep_liquidenes:
    if not os.path.isfile(enefile):
      if verbose > 0:
        print(RED + f" {os.path.basename(enefile)}: free energy file (.ene) not found!" + ENDC)
      statuslist.append(False)
    elif enefile in ene_to_sh:
      sh_path, start_snap, total_snap = ene_to_sh[enefile]
      if not os.path.isfile(sh_path):
        if verbose > 0:
          print(YELLOW + f" {os.path.basename(enefile)}: stale .ene from different lambda settings, rerun needed" + ENDC)
        statuslist.append(False)
      elif total_snap > 0 and not _bar_sh_steps_match(sh_path, start_snap, total_snap):
        if verbose > 0:
          print(YELLOW + f" {os.path.basename(enefile)}: stale .ene from different simulation time, rerun needed" + ENDC)
        statuslist.append(False)
      else:
        last_lines = _read_last_lines(enefile)
        barfinished = any("BAR Estimate of -T*dS" in line for line in last_lines)
        statuslist.append(barfinished)
    else:
      last_lines = _read_last_lines(enefile)
      barfinished = any("BAR Estimate of -T*dS" in line for line in last_lines)
      statuslist.append(barfinished)

  finished = statuslist.count(True)
  targeted = len(statuslist)
  if finished != targeted:
    if verbose > 0:
      print(YELLOW + f" Finished {finished} out of {targeted} bar analysis steps" + ENDC)
  else:
    if verbose > 0:
      print(GREEN + f" Finished {finished} out of {targeted} bar analysis steps" + ENDC)
  completed = all(statuslist)
  return completed, gasperturbsteps, gasenes, liquidperturbsteps, liquidenes, fep_gasperturbsteps, fep_gasenes, fep_liquidperturbsteps, fep_liquidenes
