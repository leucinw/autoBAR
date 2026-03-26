#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import numpy as np


def scaledownele(xyz, prm, elb):
    """Scale down electrostatic parameters in a Tinker .prm file by factor *elb*.

    Reads atom types from *xyz* and returns a list of modified parameter lines
    (as strings) for AMOEBA/AMOEBA+ electrostatic terms (multipole, polarize,
    chgtrn, chgpen, bndcflux, angcflux).  Only atoms present in the xyz file
    are modified.
    """
    prmstrs = []
    atomtypes = set(np.loadtxt(xyz, usecols=(5,), skiprows=1, dtype=str).flat)
    with open(prm) as fh:
        prmlines = fh.readlines()
    for i, line in enumerate(prmlines):
        # for AMOEBA/AMOEBA+
        # multipole
        if 'multipole' in line.lower():
            ss = line.split()
            if (ss[0].lower() == 'multipole') and (ss[1].isdigit()) and (ss[1] in atomtypes):
                if elb != 1.0:
                    prmstrs.append('  '.join(ss[:-1]) + f"{float(ss[-1])*elb:10.5f}")
                    # dipole
                    ss = prmlines[i+1].split()
                    prmstrs.append(f"        {float(ss[0])*elb:10.5f}{float(ss[1])*elb:10.5f}{float(ss[2])*elb:10.5f}")
                    # quadrupole
                    ss = prmlines[i+2].split()
                    prmstrs.append(f"        {float(ss[0])*elb:10.5f}")
                    ss = prmlines[i+3].split()
                    prmstrs.append(f"        {float(ss[0])*elb:10.5f}{float(ss[1])*elb:10.5f}")
                    ss = prmlines[i+4].split()
                    prmstrs.append(f"        {float(ss[0])*elb:10.5f}{float(ss[1])*elb:10.5f}{float(ss[2])*elb:10.5f}")

        # polarize
        if 'polarize' in line.lower():
            ss = line.split()
            if (ss[0].lower() == 'polarize') and (ss[1].isdigit()) and (ss[1] in atomtypes):
                if elb != 1.0:
                    prmstrs.append('  '.join(ss[0:2]) + f"{float(ss[2])*elb:10.5f}   " + '  '.join(ss[3:]))

        # for AMOEBA+
        # charge transfer
        if 'chgtrn' in line.lower():
            ss = line.split()
            if (ss[0].lower() == 'chgtrn') and (ss[1].isdigit()) and (ss[1] in atomtypes):
                if elb != 1.0:
                    prmstrs.append('  '.join(ss[0:2]) + f"{float(ss[2])*elb:10.5f}   " + '  '.join(ss[3:]))
        # charge penetration
        if 'chgpen' in line.lower():
            ss = line.split()
            if (ss[0].lower() == 'chgpen') and (ss[1].isdigit()) and (ss[1] in atomtypes):
                if elb != 1.0:
                    prmstrs.append('  '.join(ss[0:2]) + f"{float(ss[2])*elb:10.5f}   " + '  '.join(ss[3:]))
        # bndcflux
        if 'bndcflux' in line.lower():
            ss = line.split()
            if (ss[0].lower() == 'bndcflux') and (ss[1].isdigit()) and (ss[1] in atomtypes) and (ss[2] in atomtypes):
                if elb != 1.0:
                    prmstrs.append('  '.join(ss[0:3]) + f"{float(ss[3])*elb:10.5f}")
        # angcflux
        if 'angcflux' in line.lower():
            ss = line.split()
            if (ss[0].lower() == 'angcflux') and (ss[1].isdigit()) and (ss[1] in atomtypes) and (ss[2] in atomtypes) and (ss[3] in atomtypes):
                if elb != 1.0:
                    prmstrs.append('  '.join(ss[0:4]) + f"{float(ss[4])*elb:10.5f}{float(ss[5])*elb:10.5f}{float(ss[6])*elb:10.5f}{float(ss[7])*elb:10.5f}")
    return prmstrs
