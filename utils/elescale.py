
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import numpy as np


def _matches_keyword(tokens, keyword, atomtypes, num_type_fields=1):
  """Check if tokens represent a parameter line matching the keyword and atom types."""
  if tokens[0].lower() != keyword:
    return False
  type_fields = tokens[1:1 + num_type_fields]
  return all(t.isdigit() and t in atomtypes for t in type_fields)


def _scale_single_value_param(tokens, elb, prefix_count, trailing=True):
  """Scale the value at tokens[prefix_count] and optionally append trailing tokens."""
  scaled = "%10.5f" % (float(tokens[prefix_count]) * elb)
  prefix = '  '.join(tokens[:prefix_count])
  if trailing and len(tokens) > prefix_count + 1:
    return prefix + scaled + "   " + '  '.join(tokens[prefix_count + 1:])
  return prefix + scaled


def _scale_values(tokens, elb):
  """Scale all tokens as floats by elb and format them."""
  return ''.join("%10.5f" % (float(v) * elb) for v in tokens)


def scaledownele(xyz, prm, elb):
  if elb == 1.0:
    return []

  prmstrs = []
  atomtypes = set(np.loadtxt(xyz, usecols=(5), skiprows=1, unpack=True, dtype='str').flat)

  with open(prm) as f:
    prmlines = f.readlines()

  for i, line in enumerate(prmlines):
    tokens = line.split()
    if not tokens:
      continue
    keyword = tokens[0].lower()

    # for AMOEBA/AMOEBA+
    # multipole (spans 5 lines: monopole, dipole, quadrupole)
    if keyword == 'multipole' and _matches_keyword(tokens, 'multipole', atomtypes):
      prmstrs.append('  '.join(tokens[:-1]) + "%10.5f" % (float(tokens[-1]) * elb))
      # dipole
      dp = prmlines[i + 1].split()
      prmstrs.append('        ' + _scale_values(dp[:3], elb))
      # quadrupole
      q1 = prmlines[i + 2].split()
      prmstrs.append('        ' + _scale_values(q1[:1], elb))
      q2 = prmlines[i + 3].split()
      prmstrs.append('        ' + _scale_values(q2[:2], elb))
      q3 = prmlines[i + 4].split()
      prmstrs.append('        ' + _scale_values(q3[:3], elb))

    # polarize, charge transfer, charge penetration (same pattern: scale 3rd field)
    elif keyword in ('polarize', 'chgtrn', 'chgpen'):
      if _matches_keyword(tokens, keyword, atomtypes):
        prmstrs.append(_scale_single_value_param(tokens, elb, prefix_count=2))

    # bndcflux (2 atom-type fields, scale 4th field)
    elif keyword == 'bndcflux':
      if _matches_keyword(tokens, 'bndcflux', atomtypes, num_type_fields=2):
        prmstrs.append(_scale_single_value_param(tokens, elb, prefix_count=3, trailing=False))

    # angcflux (3 atom-type fields, scale fields 4-7)
    elif keyword == 'angcflux':
      if _matches_keyword(tokens, 'angcflux', atomtypes, num_type_fields=3):
        prefix = '  '.join(tokens[:4])
        prmstrs.append(prefix + _scale_values(tokens[4:8], elb))

  return prmstrs
