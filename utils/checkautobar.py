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


def lambda_fname(phase: str, elb: float, vlb: float) -> str:
    """Return the canonical filename stem for a lambda window."""
    return f"{phase}-e{round(elb * 100):03d}-v{round(vlb * 100):03d}"


def count_arc_snapshots(file_path):
    """Return the number of snapshots in a Tinker .arc trajectory file.

    Returns 0 if the file does not exist or is empty.
    """
    if not os.path.exists(file_path):
        return 0

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
    """Check whether all MD trajectories have reached the target snapshot count.

    Returns (completed, phase_simsnapshot) where completed is True when all
    windows are done and phase_simsnapshot maps phase names to required counts.
    """
    statuslist = []
    phase_simsnapshot = {'liquid': liquidtotalsnapshot, 'gas': gastotalsnapshot}
    for phase in phases:
        for i in range(len(orderparams)):
            elb, vlb = orderparams[i]
            fname = lambda_fname(phase, elb, vlb)
            arcfile = fname + ".arc"
            if os.path.isfile(os.path.join(homedir, phase, arcfile)):
                simsnapshot = count_arc_snapshots(os.path.join(homedir, phase, arcfile))
                if simsnapshot == phase_simsnapshot[phase]:
                    per = int(simsnapshot/phase_simsnapshot[phase]*100)
                    if verbose > 0:
                        print(GREEN + f"{fname:>20s}: " + u'\u2584'*(int(per/2)) + f" [{per:>3d}%]" + ENDC)
                    statuslist.append(True)
                else:
                    per = int(simsnapshot/phase_simsnapshot[phase]*100)
                    if verbose > 0:
                        print(YELLOW + f"{fname:>20s}: " + u'\u2584'*(int(per/2)) + f" [{per:>3d}%]" + ENDC)
                    statuslist.append(False)
            else:
                if gastotalsnapshot == 0:
                    statuslist.append(True)
                else:
                    if verbose > 0:
                        print(RED + f"{arcfile} does not exist" + ENDC)
                    statuslist.append(False)
    completed = all(statuslist)
    return completed, phase_simsnapshot


def checkbar(phases, orderparams, homedir, ignoregas, verbose):
    """Check whether all BAR analysis steps have completed.

    Returns a tuple of (completed, gasperturbsteps, gasenes, liquidperturbsteps,
    liquidenes, fep_gasperturbsteps, fep_gasenes, fep_liquidperturbsteps,
    fep_liquidenes).
    """
    statuslist = []
    gasenes = []
    liquidenes = []
    gasperturbsteps = []
    liquidperturbsteps = []
    fep_gasenes = []
    fep_liquidenes = []
    fep_gasperturbsteps = []
    fep_liquidperturbsteps = []
    checkphases = phases
    if ignoregas == 1:
        checkphases = ['liquid']

    for phase in checkphases:
        for i in range(len(orderparams)-1):
            elb0, vlb0 = orderparams[i]
            elb1, vlb1 = orderparams[i+1]
            if (elb1 > 1.0):
                elb0, vlb0 = 1.0, 1.0

            elb0s = f"{round(elb0*100):03d}"
            elb1s = f"{round(elb1*100):03d}"
            vlb0s = f"{round(vlb0*100):03d}"
            vlb1s = f"{round(vlb1*100):03d}"
            fname0 = f"{phase}-e{elb0s}-v{vlb0s}"
            fname1 = f"{phase}-e{elb1s}-v{vlb1s}"

            enedir = os.path.join(homedir, phase)

            if (int(elb1s) > 100) and (int(vlb1s) > 100):
                idx = int(int(elb1s)/10 - 10)
                enedir = os.path.join(homedir, phase, f'FEP_{idx:02d}')

            if phase == 'gas':
                enefile = fname1 + ".ene"
                fname = fname1
                if (int(elb1s) > 100) and (int(vlb1s) > 100):
                    if os.path.join(enedir, enefile) not in fep_gasenes:
                        fep_gasenes.append(os.path.join(enedir, enefile))
                        fep_gasperturbsteps.append([fname1, fname0])
                else:
                    if os.path.join(enedir, enefile) not in gasenes:
                        gasenes.append(os.path.join(enedir, enefile))
                        gasperturbsteps.append([fname1, fname0])

            if phase == 'liquid':
                enefile = fname0 + ".ene"
                fname = fname0
                if (int(elb1s) > 100) and (int(vlb1s) > 100):
                    if os.path.join(enedir, enefile) not in fep_liquidenes:
                        fep_liquidenes.append(os.path.join(enedir, enefile))
                        fep_liquidperturbsteps.append([fname0, fname1])
                else:
                    if os.path.join(enedir, enefile) not in liquidenes:
                        liquidenes.append(os.path.join(enedir, enefile))
                        liquidperturbsteps.append([fname0, fname1])

    for enefile in gasenes + liquidenes + fep_gasenes + fep_liquidenes:
        if not os.path.isfile(enefile):
            if verbose > 0:
                print(RED + " " + fname + f": free energy file (.ene) not found!" + ENDC)
            statuslist.append(False)
        else:
            with open(enefile) as fh:
                tail = deque(fh, maxlen=5)
            barfinished = any("BAR Estimate of -T*dS" in line for line in tail)
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
