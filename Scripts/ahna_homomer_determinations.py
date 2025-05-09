
#@title POKY Builder Job ID and Oligomer Information
import os
#make directory to unzip Job ID zip files
os.system (f' cd /Users/Karen/AHNA/2N74')

# JobID from POKY Builder and PDB ID
job_id = 230615_175442_137
PDB_id = '2N74'
# number of oligomers
num_oligomer = 2 

# add protons to AlphaFold PDBs by following python script from POKY

# editing and establishing AlphaFold-Multimer predictions directory 
# AF_path = os.listdir(f'/Users/Karen/AHNA/2N74/2N74_AF/')
# for predictions in AF_path:  
#   if predictions.startswith("ranked"):
#       print (predictions)

#Install Dependencies for NOESY Simulations
# pip -qq install nmrglue
# pip -qq install pynmrstar
import numpy as np
import nmrglue as ng
import pynmrstar
import json
import re
from math import sqrt

#amino acid sequence
AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

protein_attached_heavy_atoms = {
  'A': {'H':'N', 'HA':'CA', 'HB1':'CB', 'HB2':'CB', 'HB1':'CB', 'MB':'CB', 'HB':'CB'},
  'C': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB', 'HG':'S'},
  'D': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB', 'HD2':'OD2'},
  'E': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HE2':'OE2', 'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG'},
  'F': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD1':'CD1', 'HD2':'CD2', 'QD':'CQD', 'HD':'CD',
        'HE1':'CE1', 'HE2':'CE2', 'QE':'CQE', 'HE':'CE', 'HZ':'CZ'},
  'G': {'H':'N', 'HA2':'CA', 'HA1':'CA', 'QA':'CA', 'HA':'CA'},
  'H': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD1':'ND1', 'HD2':'CD2', 'HE1':'CE1', 'HE2':'NE2'},
  'I': {'H':'N', 'HA':'CA', 'HB':'CB',
        'HD11':'CD1', 'HD12':'CD1', 'HD13':'CD1', 'MD1':'CD1', 'HD1':'CD1',
        'HG12':'CG1', 'HG11':'CG1', 'QG1':'CG1', 'HG1':'CG1',
        'HG21':'CG2', 'HG22':'CG2', 'HG23':'CG2', 'MG2':'CG2', 'HG2':'CG2'},
  'K': {'H':'N', 'HA':'CA',
        'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB', 'HD2':'CD', 'HD1':'CD', 'QD':'CD', 'HD':'CD',
        'QE':'CE','HE':'CE','HE2':'CE', 'HE1':'CE', 'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG',
        'HZ1':'NZ', 'HZ2':'NZ', 'HZ1':'NZ'},
  'L': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD11':'CD1', 'HD12':'CD1', 'HD13':'CD1', 'MD1':'CD1', 'HD1':'CD1',
        'HD21':'CD2', 'HD22':'CD2', 'HD23':'CD2', 'MD2':'CD2', 'HD2':'CD2', 'QMD':'CQD', 'HD':'CD',
        'HG':'CG'},
  'M': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HE1':'CE', 'HE2':'CE', 'HE3':'CE', 'ME':'CE', 'HE':'CE',
        'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG'},
  'N': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD21':'ND2', 'HD22':'ND2', 'QD2':'ND2', 'HD2':'ND2'},
  'P': {'H2':'N', 'H3':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD2':'CD', 'HD1':'CD', 'QD':'CD', 'HD':'CD', 'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG'},
  'Q': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HE21':'NE2', 'HE22':'NE2', 'QE2':'NE2', 'HE2':'NE2',
        'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG'},
  'R': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD2':'CD', 'HD1':'CD', 'QD':'CD', 'HD':'CD',
        'HE':'NE', 'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG',
        'HH11':'NH1', 'HH12':'NH1', 'QH1':'NH1', 'HH1':'NH1',
        'HH21':'NH2', 'HH22':'NH2', 'QH2':'NH2', 'HH2':'NH2', 'QQH':'NQH', 'HH':'NH'},
  'S': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB', 'HG':'OG'},
  'T': {'H':'N', 'HA':'CA', 'HB':'CB', 'HG1':'OG1',
         'HG21':'CG2', 'HG22':'CG2', 'HG23':'CG2', 'MG2':'CG2', 'HG2':'CG2'},
  'V': {'H':'N', 'HA':'CA', 'HB':'CB',
         'HG11':'CG1', 'HG12':'CG1', 'HG13':'CG1', 'MG1':'CG1', 'HG1':'CG1',
         'HG21':'CG2', 'HG22':'CG2', 'HG23':'CG2', 'MG2':'CG2', 'HG2':'CG2', 'QMG':'CQG', 'HG':'CG'},
  'W': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
         'HD1':'CD1', 'HE1':'NE1', 'HE3':'CE3',
         'HH2':'CH2', 'HZ2':'CZ2', 'HZ3':'CZ3'},
  'Y': {'H':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD1':'CD1', 'HD2':'CD2', 'QD':'CQD', 'HD':'CD',
        'HE1':'CE1', 'HE2':'CE2', 'QE':'CQE', 'HE':'CE', 'HH':'OH'},
}

# change of nomenclature
s12s32 = (
 ('A', 'HB1',  'HB3'),
 ('G', 'HA1',  'HA3'),
 ('C', 'HB1',  'HB3'),
 ('D', 'HB1',  'HB3'),
 ('D', 'H' ,   'HN' ),
 ('E', 'HB1',  'HB3'),
 ('E', 'HG1',  'HG3'),
 ('F', 'HB1',  'HB3'),
 ('H', 'HB1',  'HB3'),
 ('I', 'HG11', 'HG13'),
 ('K', 'HB1',  'HB3'),
 ('K', 'HD1',  'HD3'),
 ('K', 'HG1',  'HG3'),
 ('K', 'HE1',  'HE3'),
 ('K', 'HZ1',  'HZ3'),
 ('L', 'HB1',  'HB3'),
 ('M', 'H' ,   'H1'),
 ('M', 'H' ,   'H2'),
 ('M', 'H' ,   'H3'),
 ('M', 'HB1',  'HB3'),
 ('M', 'HG1',  'HG3'),
 ('M', 'H'  ,  'HN' ),
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
 ('Y', 'HB1',  'HB3')
)

# changing from 2/3 nomeclature to 1/2

def s32tos12(szA, szHX):
  for grp in s12s32:
    if (grp[0] == szA) and (grp[2] == szHX):
      return grp[1]

  return szHX

#NOESY Simulation: Setting Up Parsing BMRB and PDB
############################## PARSING BMRB ####################################

AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
  'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
  'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
  'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
  'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

def read_bmrb(fname=None, ename=None):
  try:
    if fname != None:
      strentry = pynmrstar.Entry.from_file(fname)
    elif ename != None:
      strentry = pynmrstar.Entry.from_database(ename)
    else:
      return {}
  except:
    return {}

  # Title
  info_saveframes = strentry.get_saveframes_by_category('entry_information')
  title = info_saveframes[0]['Title'][0].strip()

  # extract sequences
  sequence_loops = strentry.get_loops_by_category('_Entity_comp_index')
  sequence_dict = json.loads(sequence_loops[0].get_json())
  sequence_list = list(map(lambda x: [x[1], x[2], AAA_dict[x[2]]],
                           sequence_dict['data']))
  fasta = ''.join(list(map(lambda x: x[2], sequence_list)))

  # extract chemical shifts
  atom_chem_shift_loops = strentry.get_loops_by_category('atom_chem_shift')

  # First, try from bond and chemical shift information
  # shift_list includes a set of lists:
  #     [nseq, 3-letter-seq, 1-letter-seq, chemical shift]
  atom_chem_shift_dict = json.loads(atom_chem_shift_loops[0].get_json())
  shift_list = list(map(lambda x: [x[4], x[6], AAA_dict[x[6]], x[7], x[10]],
                              atom_chem_shift_dict['data']))

  return {'Title': title, 'Fasta': fasta,
          'Sequence List': sequence_list,
          'Chemical Shifts': shift_list}
results = read_bmrb(ename=25793)
print (results)
raise SystemExit
############################## PARSING PDB #####################################

def distance3D(atom1, atom2):
  """ takes two coordinates. ex: ((26.266, 25.413, 2.842),
                                  (26.913, 26.639, -3.51))
      returns the distance
  """
  return sqrt((atom1[0] - atom2[0]) ** 2 +
              (atom1[1] - atom2[1]) ** 2 +
              (atom1[2] - atom2[2]) ** 2)
def s32tos12(szA, szHX):
  for grp in s12s32:
    if (grp[0] == szA) and (grp[2] == szHX):
      return grp[1]

  return szHX

def readPDB(pdb_file, modelnumber):
  """
  Parameters
  ----------
  uploaded_pdb : address to a .pdb file
  modelnumber : the model number in the pdb (starting from 1)
  Returns
  -------
  pdb_list: a list of each [Nseq, atom's coordinates, aa, atom]
  """
  pdbLines, modelList = [], []
  tempLines = open(pdb_file, 'r').readlines()

  # clean
  for line in tempLines:
    #print(line)
    if line[0:4] in ['MODE', 'ATOM', 'ENDM']:
      pdbLines.append(line)

  # fill modelList
  for line in pdbLines:
    if line[0:5] == 'MODEL':
      modelList.append([])
    if line[0:4] != 'ATOM':
      continue
    if line[12:16].strip()[0] not in ['H']:
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

    modelList[-1].append( [nSeq, x, y, z, aaa, AAA_dict[aaa], atm, c] )


  return modelList[modelnumber-1]

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

def get_shift(seqidx, atomname):
  filterAtoms = list(filter(lambda x :x[0] == str(seqidx) and x[3] ==atomname,
                            chemShifts))
  if len(filterAtoms) == 0:
    return -9999
  return float(filterAtoms[0][4])

def distance2height(HHdist, offset = 0.0):
    """ distance scaling using r^-6 approximation for H-H """

    min_hh = 1.7  # closest distance between H atoms (1.70 A)
    max_hh = 5.5  # farthest distance observed between H atoms

    max_I = 10**6 # arbitrary value
    min_I = 10**4 # arbitrary value

    approx = -6.0

    dist = min(HHdist, max_hh-offset)
    dist = max(dist, min_hh-offset)

    # A * r**approx + B = I
    # A * min_hh**approx + B = max_I    arbituary max
    # A * max_hh**approx + B = min_I    arbituary min

    # A * (min_hh**approx - max_hh**approx) = max_I - min_I
    #

    A = (max_I - min_I) / ((min_hh-offset)**approx - (max_hh-offset)**approx)       # 75227372.04997641
    B = min_I - A * (max_hh-offset)**approx                     # -19444.257193289242

    #print(A, B, dist, approx, offset)
    return A * dist**approx + B

#@title NOESY Simulation: Simulate NOESY Peak List
regex = r"\d+"
c_list = 'ABCD'

def parse_key(key):
    c = key.split('_')[0]
    a = key.split('_')[1][0]
    nSeq = int(re.search(regex, key.split('_')[1]).group())
    sep = c+'_'+a+str(nSeq)
    atm = key.split(sep)[1]
    return c, a, nSeq, atm
noesyType = 'nnoe' # ["nnoe", "cnoe"]
def create_noesy_peak_list(keyList,noesy_type={noesyType}):
  # key: c+'_'+a+str(nSeq)+atm
  cnt = 0
  lines = ''
  for key in keyList:
    # break down keyList
    c, a, nSeq, atm = parse_key(key)
    atm = s32tos12(a,atm)
    # getting heavy atom 'N' or 'C'
    ncatm = protein_attached_heavy_atoms[a][atm]
    if noesy_type == 'nnoe':
      if ncatm[0] != 'N':
        continue
    elif noesy_type == 'cnoe':
      if ncatm[0] != 'C':
        continue

    # get heavy atom shift
    nc_shift = get_shift(nSeq, ncatm)

    # get proton atom shift
    h_shift = get_shift(nSeq, atm)

    if nc_shift < -1000 or h_shift < -1000:
      continue

    grp1 = f'{a}{nSeq}{ncatm}'
    grp2 = f'{atm}'

    for key2 in keyList:
      c2, a2, nSeq2, atm2 = parse_key(key2)
      dist = distMat[nSeq, nSeq2, c_list.index(c), c_list.index(c2)]
      if key != key2 and dist == 0:
        continue
      if dist > 5:
        continue
      h2_shift = get_shift(nSeq2, atm2)

      if h2_shift < -1000:
        continue

      height = distance2height(dist)
      if nSeq != nSeq2:
        grp3 = f'{a2}{nSeq2}{atm2}'
      else:
        grp3 = f'{atm2}'

      asgn = f'{grp1}-{grp2}-{grp3}'
      line = '  %20s %8.3f %8.3f %8.3f %15d' % (asgn, nc_shift, h_shift,
                                                h2_shift, int(height))
      #print(line)
      lines += line + '\n'
    return lines
############################# Write Peak List ##################################
experimental_pdb_file = '/AHNA/2N74' + str(PDB_ID) + '.pdb'
Alt_pdb_file = '/AHNA/2N74/2N74_AF'

for PDB in os.listdir(pdb_dir_alt):
  interfacePath = os.path.join(pdb_dir_alt, PDB)
  interfacePDB = interfacePath
  pdb_list_alt = readPDB(interfacePDB,1)
  pdb_list = readPDB(experimental_pdb_file, 0)

#createDistanceMatrix(pdb_list)
distMat, keyList = createDistanceMatrix(pdb_list)
chemShifts = results['Chemical Shifts']

#createDistanceMatrix(pdb_list_alt)
distMat_Alt, keyList_alt = createDistanceMatrix(pdb_list_alt)


line = create_noesy_peak_list(keyList)
line_alt = create_noesy_peak_list(keyList_alt)

PDB_text =open (f'NOESY_{PDB_ID}.txt', 'w')
PDB_text.write(line)
PDB_text.close()

Alt_list =open (f'{interfacePDB[:-4]}.txt', 'w')
Alt_list.write(line_alt)
Alt_list.close()
if os.path.exists(Alt_pdb_file):
  for files in os.listdir(Alt_pdb_file):
    pdb_dir_alt = os.path.join(Alt_pdb_file, files)
    print (pdb_dir_alt)
    for PDB in os.listdir(pdb_dir_alt):
      if PDB.endswith('.pdb'):
        interfacePath = os.path.join(pdb_dir_alt, PDB)
        interfacePDB = interfacePath
        print (interfacePDB)

        pdb_list_alt = readPDB(interfacePDB, 1)
        pdb_list = readPDB(experimental_pdb_file, 1)

        #createDistanceMatrix(pdb_list)
        distMat, keyList = createDistanceMatrix(pdb_list)
        chemShifts = results['Chemical Shifts']

        #createDistanceMatrix(pdb_list_alt)
        distMat_Alt, keyList_alt = createDistanceMatrix(pdb_list_alt)

        line = create_noesy_peak_list(keyList)
        line_alt = create_noesy_peak_list(keyList_alt)

        PDB_text =open (f'NOESY_{PDB_ID}.txt', 'w')
        PDB_text.write(line)
        PDB_text.close()

        Alt_list =open (f'{interfacePDB[:-4]}.txt', 'w')
        Alt_list.write(line_alt)
        Alt_list.close()

#@title Compare NOESY Peak Lists

def compareNOESY (PDB, ALT):

  mat1 = np.zeros((500,160,360))  # bucket size -> 0.1 for N (90-140)
  mat2 = np.zeros((500,160,360))  # bucket size -> 0.05 for HN (6 - 14)
                                  # bucket size -> 0.05 for H (-4 - 14)

  for d in PDB:
    try:
      i1 = int((d[0]-90)/0.1 +0.5) # subtracting by minimal values
      i2 = int((d[1]-6)/0.05 +0.5)
      i3 = int((d[2]+4)/0.05 +0.5)

      mat1[i1,i2,i3] = 1

    except:
      pass

  for d in ALT:
    try:
      i1 = int((d[0]-90)/0.1 +0.5)
      i2 = int((d[1]-6)/0.05 +0.5)
      i3 = int((d[2]+4)/0.05 +0.5)

      mat2[i1,i2,i3] = 1

    except:
      pass

  compared_mat = np.logical_and(mat1, mat2)
  return np.sum(compared_mat)

PDB_text = '/content/NOESY_2N74.txt'
Alt_text = '/content/interfacePDB/content'
for text in os.listdir(Alt_text):
    txt_dir_alt = os.path.join(Alt_text, text)
    if txt_dir_alt.endswith('.txt'):

      with open(PDB_text) as f1, open(txt_dir_alt) as f2:
  # Read data directly into NumPy arrays for efficiency
        PDB = np.genfromtxt(f1, delimiter='', usecols = (1,2,3)) #peak list NOESY
        ALT = np.genfromtxt(f2, delimiter='', usecols=(1,2,3)) # peak list Alt NOESY
        similarities = compareNOESY(PDB, ALT)
        print (txt_dir_alt,similarities)

#@title Creating Interface .tbl File from Alphafold-Multimer

import numpy as np
from math import sqrt

#amino acid sequence
AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
# change of nomenclature
s12s32 = (
 ('A', 'HN'  ,  'H' ),
 ('N', 'HN'  ,  'H' ),
 ('G', 'HN'  ,  'H' ),
 ('G', 'HA1',  'HA3'),
 ('C', 'HB1',  'HB3'),
 ('D', 'HN'  ,  'H' ),
 ('D', 'HB1',  'HB3'),
 ('E', 'HB1',  'HB3'),
 ('E', 'HG1',  'HG3'),
 ('E', 'HN',   'H'),
 ('F', 'HN'  ,  'H' ),
 ('F', 'HB1',  'HB3'),
 ('H', 'HB1',  'HB3'),
 ('I', 'HN'  ,  'H' ),
 ('I', 'HG11', 'HG13'),
 ('K', 'HB1',  'HB3'),
 ('K', 'HD1',  'HD3'),
 ('K', 'HG1',  'HG3'),
 ('K', 'HE1',  'HE3'),
 ('L', 'HN'  ,  'H' ),
 ('L', 'HB1',  'HB3'),
 ('M', 'HB1',  'HB3'),
 ('M', 'HG1',  'HG3'),
 ('N', 'HB1',  'HB3'),
 ('P', 'HB1',  'HB3'),
 ('P', 'HD1',  'HD3'),
 ('P', 'HG1',  'HG3'),
 ('Q', 'HB1',  'HB3'),
 ('Q', 'HG1',  'HG3'),
 ('R', 'HN'  ,  'H' ),
 ('R', 'HB1',  'HB3'),
 ('R', 'HD1',  'HD3'),
 ('R', 'HG1',  'HG3'),
 ('T', 'HN'  ,  'H' ),
 ('V', 'HN'  ,  'H' ),
 ('S', 'HN'  ,  'H' ),
 ('S', 'HB1',  'HB3'),
 ('W', 'HN'  ,  'H' ),
 ('W', 'HB1',  'HB3'),
 ('Y', 'HB1',  'HB3'),
)
# distance cutoff
cutoffDist = "H_5.5" # @param ["H_5.5", "C_7.7", "N_7.7"] {allow-input: true}
cutoffDist= cutoffDist.split('_')[1]
model_number = 1 #@param {type:"slider", min: 1, max: 4}

def distance3D(atom1, atom2):
  """ takes two coordinates. ex: ((26.266, 25.413, 2.842),
                                  (26.913, 26.639, -3.51))
      returns the distance
  """
  return sqrt((atom1[0] - atom2[0]) ** 2 +
              (atom1[1] - atom2[1]) ** 2 +
              (atom1[2] - atom2[2]) ** 2)

def readPDB(pdb_file, modelnumber):
  """
  Parameters
  ----------
  uploaded_pdb : address to a .pdb file
  modelnumber : the model number in the pdb (starting from 1)
  Returns
  -------
  pdb_list: a list of each [atom's coordinates, aa, atom]
  """
  pdbLines, modelList = [], []
  tempLines = open(pdb_file, 'r').readlines()

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
    if nSeq < 100 and c == 'A':
      modelList[-1].append( [nSeq, x, y, z, aaa, AAA_dict[aaa], atm, c] )
    if nSeq > 100 and c !='A':
      nSeqB = nSeq - 100
      modelList[-1].append( [nSeqB, x, y, z, aaa, AAA_dict[aaa], atm, c] )

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

# changing from 2/3 nomeclature to 1/2

def s32tos12(szA, szHX):
  for grp in s12s32:
    if (grp[0] == szA) and (grp[2] == szHX):
      return grp[1]

  return szHX

############################ Write .tbl File ###################################
altPath = '/content/interfacePDB/content/230615175442137_e8e3c_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_000_H.pdb' #@param
#@markdown Path of Alphafold-Multimer with the most similarities.
pdb_list= readPDB(altPath, int(model_number))
distMat, keyList = createDistanceMatrix(pdb_list)
line = ''
newDist = ''
for i in range(len(pdb_list)):
  nSeq, x, y, z, aaa, a, atm, c = pdb_list[i]
  if atm[0] not in ['H', 'N', 'O']:
    continue

  for j in range(i+1, len(pdb_list)):
    nSeq2, x2, y2, z2, aaa2, a2, atm2, c2 = pdb_list[j]
    dist = distance3D( (x, y, z), (x2, y2, z2))
    if atm2[0] not in ['H', 'N', 'O']:
      continue

    if 'H' not in [atm[0], atm2[0]]:
      continue

    if distMat[i, j, 0, 1] > float(cutoffDist):
      continue

    if 'H' in [atm[0], atm2[0]] and 'O' in [atm[0], atm2[0]]:
      if distMat[i, j, 0, 1] > 2.2:
        continue

    if 'N' in [atm[0], atm2[0]] and 'O' in [atm[0], atm2[0]]:
      if distMat[i, j, 0 , 1] > 3.3:
        continue

    if distMat[i, j, 0, 1] > 0.01 and distMat[i, j, 0, 1] < float(cutoffDist):
      catm = s32tos12(a, atm)
      catm2 = s32tos12(a2,atm2)
      dist = distMat[i, j, 0,1] # giving 'rubber band' 40% more resistance
      #print (dist)
      distRound = round(dist, 2)
      dist2 = round((distRound - 1.7),2) # 1.7A is from VDW

      line += f'''assign (segid {c} and resid {nSeq} and name {catm}) \
                  (segid {c2} and resid {nSeq2} and name {catm2}) \
                  {distRound} {dist2} 0.0  \n'''
      print(line)
      interface_path_tbl = '/content/AHNA/BestEvaluated/alt.tbl'
      f = open(interface_path_tbl, 'w')
      f.write(line)
      f.close()


