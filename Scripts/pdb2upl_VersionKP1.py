#
# This is an example script to simulate upl file from a PDB file.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module (CTRL+B)
#

import __main__
s = __main__.main_session

import numpy as np
from math import sqrt
from myseq import AAA_dict

# pdb
pdb_file = s.open_filedialog('Select a PDB file.',
                             'PDB file (*.pdb)', '')

if pdb_file == '' or not pdb_file.endswith('.pdb'):
  print("You didn't select a PDB file")
  raise SystemExit

model_number = s.show_inputdialog('Model number',
                                  'Select the model number in PDB',
                                  '1')

# distance cutoff
ans = \
    s.show_inputdialog('Cutoff', 'Distance cutoff (e.g. H_5.5, C_7.7 or N_7.7)', 
                      'H_5.5') 
try:
  cutoffAtom = ans.split('_')[0]
  cutoffDist = float(ans.split('_')[1])
except:
  raise SystemExit
#amino acid sequence
AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
def distance3D(atom1, atom2):
  """ takes two coordinates. ex: ((26.266, 25.413, 2.842),
                                  (26.913, 26.639, -3.51))
      returns the distance
  """
  return sqrt((atom1[0] - atom2[0]) ** 2 +
              (atom1[1] - atom2[1]) ** 2 +
              (atom1[2] - atom2[2]) ** 2)

def readPDB(pdb, modelnumber):
  """
  Parameters
  ----------
  pdb : address to a .pdb file
  modelnumber : the model number in the pdb (starting from 1)
  Returns
  -------
  pdb_list: a list of each [atom's coordinates, aa, atom]
  """
  pdbLines, modelList = [], []
  tempLines = open(pdb, 'r').readlines()
  
  # clean
  for line in tempLines:
    if line[0:4] in ['MODE', 'ATOM', 'ENDM']:
      pdbLines.append(line)

  # fill modelList
  for line in pdbLines:
    if line[0:5] == 'MODEL':
      modelList.append([])
    if line[0:4] != 'ATOM':
      continue
    if line[12:16].strip()[0] not in ['C', 'N', 'H', 'O']:
        continue
    aaa = line[17:20].strip()
    atm = line[12:16].strip()
    nSeq = int(line[23:26].strip()) # <residue sequence>
    x = float(line[30:38].strip())
    y = float(line[38:46].strip())
    z = float(line[46:54].strip())
    c = line[21:22].strip() # <chain_ID>

    # in case MODEL not in PDB.
    if len(modelList) == 0:
      modelList.append([])
    try:
      modelList[-1].append( [nSeq, x, y, z, aaa, AAA_dict[aaa], atm, c] )
    except:
      continue

  return modelList[modelnumber - 1]

def createDistanceMatrix(pdb_list):
  distMat = np.zeros( (len(pdb_list), len(pdb_list), 4, 4 ) )
  keyList = []
  c_list = 'ABCD'
  for i in range(len(pdb_list)):
    nSeq, x, y, z, aaa, a, atm, c = pdb_list[i]
    keyList.append(c+'_'+a+str(nSeq)+atm)
    for j in range(i+1, len(pdb_list)):
      nSeq2, x2, y2, z2, aaa2, a2, atm2, c2 = pdb_list[j]
      dist = distance3D( (x, y, z), (x2, y2, z2))
      distMat[i, j, c_list.index(c), c_list.index(c2)] = \
        distMat[j, i, c_list.index(c2), c_list.index(c)] = dist
  return distMat, keyList

# read PDB
pdb_list = readPDB(pdb_file, int(model_number))
print(len(pdb_list))
distMat, keyList = createDistanceMatrix(pdb_list)
peak_list = []

#methyls
s12s32 = (
 ('G', 'HA1',  'HA3'),
 ('C', 'HB1',  'HB3'),
 ('D', 'HB1',  'HB3'),
 ('E', 'HB1',  'HB3'),
 ('E', 'HG1',  'HG3'),
 ('F', 'HB1',  'HB3'),
 ('H', 'HB1',  'HB3'),
 ('I', 'HG11', 'HG13'),
 ('K', 'HB1',  'HB3'),
 ('K', 'HD1',  'HD3'),
 ('K', 'HG1',  'HG3'),
 ('K', 'HE1',  'HE3'),
 ('L', 'HB1',  'HB3'),
 ('M', 'HB1',  'HB3'),
 ('M', 'HG1',  'HG3'),
 ('M', 'H',    'H1' ),
 ('M', 'H',    'H2' ),
 ('M', 'H',    'H3' ),
 ('N', 'HB1',  'HB3'),
 ('P', 'HB1',  'HB3'),
 ('P', 'HD1',  'HD3'),
 ('P', 'HG1',  'HG3'),
 ('Q', 'HB1',  'HB3'),
 ('Q', 'HG1',  'HG3'),
 ('R', 'HB1',  'HB3'),
 ('R', 'HD1',  'HD3'),
 ('R', 'HG1',  'HG3'),
 ('S', 'HB1',  'HB3'),
 ('W', 'HB1',  'HB3'),
 ('Y', 'HB1',  'HB3'),
)

# changing from 2/3 nomeclature to 1/2

def s32tos12(szA, szHX):
  for grp in s12s32:
    if (grp[0] == szA) and (grp[2] == szHX):
      return grp[1]

  return szHX

content = ''
for i in range(len(keyList)):
  key = keyList[i]
  nSeq, x, y, z, aaa, a, atm, c = pdb_list[i]
  if atm[0] not in ['H', 'N', 'O']:
    continue

  for j in range(i+1, len(keyList)):
    nSeq2, x2, y2, z2, aaa2, a2, atm2, c2 = pdb_list[j]
    if atm2[0] not in ['H', 'N', 'O']:
      continue

    if atm2[0] not in ['H', 'N', 'O']:
      continue

    if 'H' not in [atm[0], atm2[0]]:
      continue

    if distMat[i, j, 0, 0] > float(cutoffDist):
      continue

    if 'H' in [atm[0], atm2[0]] and 'O' in [atm[0], atm2[0]]:
      if distMat[i, j, 0, 0] > 2.2:
        continue

    if 'N' in [atm[0], atm2[0]] and 'O' in [atm[0], atm2[0]]:
      if distMat[i, j, 0 , 0] > 3.3:
        continue

    if distMat[i, j, 0, 0] > 0.01 and distMat[i, j, 0, 0] < float(cutoffDist):
      catm = s32tos12(a, atm)
      catm2 = s32tos12(a2,atm2)
      dist = distMat[i, j, 0,0] * 1.4

      line = '%3d %3s %4s %3d %3s %4s %.3f' % \
       (nSeq, aaa, catm, nSeq2, aaa2, catm2, dist)
      print(line)
      content += f'{line}\n'



out_file = s.save_filedialog('Save as.',
                             'UPL file (*.upl)', '')
if out_file != '':
  f = open(out_file, 'w')
  f.write(content)
  f.close()