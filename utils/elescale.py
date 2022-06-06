
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import numpy as np

def scaledownele(xyz, prm, elb):
  prmstrs = []
  atomtypes = np.loadtxt(xyz, usecols=(5), skiprows=1, unpack=True, dtype='str')
  prmlines = open(prm).readlines()
  for i in range(len(prmlines)):
    # for AMOEBA/AMOEBA+
    # multipole
    if 'multipole' in prmlines[i].lower():
      ss = prmlines[i].split()
      if (ss[0].lower() == 'multipole') and (ss[1].isdigit()) and (ss[1] in atomtypes):
        if elb != 1.0:
          prmstrs.append('  '.join(ss[:-1]) +  "%10.5f"%(float(ss[-1])*elb))
          # dipole
          ss = prmlines[i+1].split()
          prmstrs.append('        %10.5f%10.5f%10.5f'%(float(ss[0])*elb, float(ss[1])*elb,float(ss[2])*elb))
          # quadrupole 
          ss = prmlines[i+2].split()
          prmstrs.append('        %10.5f'%(float(ss[0])*elb))
          ss = prmlines[i+3].split()
          prmstrs.append('        %10.5f%10.5f'%(float(ss[0])*elb, float(ss[1])*elb))
          ss = prmlines[i+4].split()
          prmstrs.append('        %10.5f%10.5f%10.5f'%(float(ss[0])*elb, float(ss[1])*elb,float(ss[2])*elb))

    # polarize
    if 'polarize' in prmlines[i].lower():
      ss = prmlines[i].split()
      if (ss[0].lower() == 'polarize') and (ss[1].isdigit()) and (ss[1] in atomtypes):
        if elb != 1.0:
          prmstrs.append('  '.join(ss[0:2]) +  "%10.5f   "%(float(ss[2])*elb) + '  '.join(ss[3:]))

    # for AMOEBA+
    # charge transfer
    if 'chgtrn' in prmlines[i].lower():
      ss = prmlines[i].split()
      if (ss[0].lower() == 'chgtrn') and (ss[1].isdigit()) and (ss[1] in atomtypes):
        if elb != 1.0:
          prmstrs.append('  '.join(ss[0:2]) +  "%10.5f   "%(float(ss[2])*elb) + '  '.join(ss[3:]))
    # charge penetration
    if 'chgpen' in prmlines[i].lower():
      ss = prmlines[i].split()
      if (ss[0].lower() == 'chgpen') and (ss[1].isdigit()) and (ss[1] in atomtypes):
        if elb != 1.0:
          prmstrs.append('  '.join(ss[0:2]) +  "%10.5f   "%(float(ss[2])*elb) + '  '.join(ss[3:]))
    # bndcflux
    if 'bndcflux' in prmlines[i].lower():
      ss = prmlines[i].split()
      if (ss[0].lower() == 'bndcflux') and (ss[1].isdigit()) and (ss[1] in atomtypes) and (ss[2] in atomtypes):
        if elb != 1.0:
          prmstrs.append('  '.join(ss[0:3]) +  "%10.5f"%(float(ss[3])*elb))
    # angcflux
    if 'angcflux' in prmlines[i].lower():
      ss = prmlines[i].split()
      if (ss[0].lower() == 'angcflux') and (ss[1].isdigit()) and (ss[1] in atomtypes) and (ss[2] in atomtypes) and (ss[3] in atomtypes):
        if elb != 1.0:
          prmstrs.append('  '.join(ss[0:4]) +  "%10.5f%10.5f%10.5f%10.5f"%(float(ss[4])*elb,float(ss[5])*elb,float(ss[6])*elb,float(ss[7])*elb))
  return prmstrs
