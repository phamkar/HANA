import os
import numpy as np
import json
import re
import math
from math import sqrt
from itertools import product
import proton_nomenclature 

############################################################################################################################################################
#                                                                    PACSY Abundance                                                                       #                                         
############################################################################################################################################################
m_pAllHX = ["H", "HN", "HA","HA2","HA3","HB","HB1","HB2","HB3",
            "HD", "HD1","HD11","HD12","HD13",
            "HD2","HD21","HD22","HD23", "HD3","HE","HE1","HE2","HE3","HE21",
            "HE22","HG","HG1","HG11", "HG12","HG13","HG2","HG21","HG22",
            "HG23","HG3","HH11","HH12","HH2","HH21","HH22","HZ","HZ2","HZ3"]
m_pHeavyNC = ["N","NE","NE1","NE2","ND1","ND2","NH1",
              "NH2","NZ","CA","CB","CD","CD1","CD2",
              "CE","CE1","CE2","CE3","CG","CG1","CG2",
              "CH2","CZ","CZ2","CZ3"]
m_pA = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W",
        "Y","V"]
m_pAAA = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE", "LEU",
          "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]

AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
     'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

def AA2Index(a):
  try:
    return m_pA.index(a)
  except:
    try:
      return m_pAAA.index(a)
    except:
      pass
  return -1

def GetAttachedProtons(szAAA, szHeavyAtom):
    attached_protons = {
        "N": {"PRO": [], "default": ["H"]},
        "CA": {"GLY": ["HA2", "HA3"], "default": ["HA"]},
        "CB": {
            "THR": ["HB2", "HB3"],
            "ILE": ["HB2", "HB3"],
            "VAL": ["HB2", "HB3"],
            "default": ["HB"],
        },
        "ALA": {"default": []},
        "ARG": {
            "CG": ["HG2", "HG3", "HG"],
            "CD": ["HD2", "HD3", "HD"],
            "NE": ["HE"],
            "NH1": ["HH11", "HH12", "HH1"],
            "NH2": ["HH21", "HH22", "HH2"],
            "NH": ["HH1", "HH2", "HH"],
            "default": [],
        },
        "ASN": {"ND2": ["HD21", "HD22", "HD2"], "default": []},
        "ASP": {"default": []},
        "CYS": {"default": []},
        "GLN": {
            "NE2": ["HE21", "HE22", "HE2"],
            "CG": ["HG2", "HG3", "HG"],
            "default": [],
        },
        "GLU": {"CG": ["HG2", "HG3", "HG"], "default": []},
        "GLY": {"default": []},
        "HIS": {
            "ND1": ["HD1"],
            "NE1": ["HE1"],
            "CD2": ["HD2"],
            "CE1": ["HE1"],
            "default": [],
        },
        "ILE": {
            "CG1": ["HG12", "HG13", "HG1"],
            "CG2": ["HG2"],
            "CD1": ["HD1"],
            "default": [],
        },
        "LEU": {"CG": ["HG"], "CD1": ["HD1"], "CD2": ["HD2"], "default": []},
        "LYS": {
            "NZ": ["HZ"],
            "CG": ["HG2", "HG3", "HG"],
            "CD": ["HD2", "HD3", "HD"],
            "CE": ["HE2", "HE3", "HE"],
            "default": [],
        },
        "MET": {"CE": ["HE"], "CG": ["HG2", "HG3", "HG"], "default": []},
        "PHE": {
            "CD1": ["HD1"],
            "CD2": ["HD2"],
            "CD": ["HD1", "HD2", "HD"],
            "CE1": ["HE1"],
            "CE2": ["HE2"],
            "CE": ["HE1", "HE2", "HE"],
            "CZ": ["HZ"],
            "default": [],
        },
        "PRO": {"CG": ["HG2", "HG3", "HG"], "CD": ["HD2", "HD3", "HD"]},
        "SER": {"default": []},
        "THR": {"CG2": ["HG2"], "default": []},
        "TRP": {
            "NE1": ["HE1"],
            "CD1": ["HD1"],
            "CD2": ["HD2"],
            "CD": ["HD1", "HD2", "HD"],
            "CE3": ["HE3"],
            "CZ2": ["HZ2"],
            "CZ3": ["HZ3"],
            "CZ": ["HZ2", "HZ3", "HZ"],
            "CH2": ["HH2"],
            "default": [],
        },
        "TYR": {
            "CD1": ["HD1"],
            "CD2": ["HD2"],
            "CD": ["HD1", "HD2", "HD"],
            "CE1": ["HE1"],
            "CE2": ["HE2"],
            "CE": ["HE1", "HE2", "HE"],
            "default": [],
        },
        "VAL": {
            "CG1": ["HG1"],
            "CG2": ["HG2"],
            "CG": ["HG1", "HG2", "HG"],
            "default": [],
        },
    }

    return attached_protons.get(szHeavyAtom, {}).get(szAAA,
                attached_protons[szHeavyAtom].get("default", []))


def GetAttachedHeavyAtom(szAAA, szHAtom):
  if szHAtom in ["HN", "H"]:
      return "N"
  if len(szHAtom) < 2:
    return ""

  szPostfix = szHAtom[1:]
  szPostfix2 = szHAtom[1:-1]

  for i, szHeavyAtom in enumerate(["C" + szPostfix, "N" + szPostfix,
                                   "C" + szPostfix2, "N" + szPostfix2]):
    attachList = GetAttachedProtons(szAAA, szHeavyAtom)
    if i > 1 and len(szHAtom < 3):
       return ''
    if szHAtom in attachList:
      return szHeavyAtom

  return ''

def GetIndexForAllProton(szAtom):
  try:
    return m_pAllHX.index(szAtom)
  except:
    return -1

def GetPacsyAbundance(szA, szAtom, szA2, szAtom2):
  iA = AA2Index(szA)
  iA2 = AA2Index(szA2)
  iH = GetIndexForAllProton(szAtom)
  iH2 = GetIndexForAllProton(szAtom2)

  if -1 in [iA, iA2, iH, iH2]:
    return 0

  f = open(f'/Users/Karen/AHNA/Scripts/PacsyCloseQRef/pacsymat/pacsyoccurrence_{szA}_all.txt', 'r')
  lines = f.readlines()
  f.close()
  line = lines[iA2*43+iH]
  return int(line.strip().split(',')[iH2])

############################################################################################################################################################
#                                                                 Parsing Chemical Shifts                                                                  # 
############################################################################################################################################################

def read_pokybuilder_prot(prot_path):
# This function read the *.prot file, which gives the atoms and their chemical shifts existing within the compound. 
# The function returns a list and dictionary of the sequence number, atom name, and chemical shifts.
  f = open(prot_path, 'r')
  lines = f.readlines()
  f.close() 

  # start parsing
  cs_list = []
  cs_dict ={}
  for line in lines:
    sp_list = line.strip().split() #  1 64.2 0.200 CA 1
    try:
      if sp_list[0][0] == '#':
        continue
      cs, atm, nseq = float(sp_list[1]), sp_list[3], int(sp_list[4])
      cs_list.append([nseq, atm, cs])
      cs_dict[nseq,atm] = [cs]

    except:
      continue

  return cs_dict, cs_list

############################################################################################################################################################
#                                                                 Proton Nomenclature                                                                      # 
############################################################################################################################################################
def protonNomenclatureChange (aa, atom):
  prot = read_pokybuilder_prot(prot_path)
  try:
    prot_aa = proton_nomenclature.protonNomenclature[aa] # 'V' -> {'HG':['HG11','HG12','HG13', 'HG21','HG22','HG23']
    nomenclature = prot_aa[atom] # 'HG' -> ['HG11','HG12','HG13', 'HG21','HG22','HG23']
  except:
    return [atom,]
  return nomenclature

############################################################################################################################################################
#                                                              NOESY Simulation: Parsing PDB                                                               #                                                   #
############################################################################################################################################################

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


############################################################################################################################################################
#                                                              NOESY Simulation: Simulate NOESY Peak List                                                  #                                                               
############################################################################################################################################################

def distance3D(atom1, atom2):
  """ takes two coordinates. ex: ((26.266, 25.413, 2.842),
                                  (26.913, 26.639, -3.51))
      returns the distance
  """
  return sqrt((atom1[0] - atom2[0]) ** 2 +
              (atom1[1] - atom2[1]) ** 2 +
              (atom1[2] - atom2[2]) ** 2)

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

def get_shift(seqidx, atomname, chem_shifts):
  filterAtoms = list(filter(lambda x :x[0] == seqidx and x[1] ==atomname,
                            chem_shifts))
  if len(filterAtoms) == 0:
    return -9999
  return float(filterAtoms[0][2])

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

def parse_key(key):
    c = key.split('_')[0]
    a = key.split('_')[1][0]
    nSeq = int(re.search(regex, key.split('_')[1]).group())
    sep = c+'_'+a+str(nSeq)
    atm = key.split(sep)[1]
    return c, a, nSeq, atm

def create_noesy_peak_list(keyList,noesy_type,chem_shifts):
  # key: c+'_'+a+str(nSeq)+atm
  lines = ''
  for key in keyList:
    # break down keyList
    c, a, nSeq, atm = parse_key(key)
    atm = proton_nomenclature.s32tos12(a,atm)
    # getting heavy atom 'N' or 'C'
    ncatm = proton_nomenclature.protein_attached_heavy_atoms[a][atm]
    if noesy_type == 'nnoe':
      if ncatm[0] != 'N':
        continue
    elif noesy_type == 'cnoe':
      if ncatm[0] != 'C':
        continue

    # get heavy atom shift
    nc_shift = get_shift(nSeq, ncatm, chem_shifts)

    # get proton atom shift
    h_shift = get_shift(nSeq, atm, chem_shifts)

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
      h2_shift = get_shift(nSeq2, atm2, chem_shifts)

      if h2_shift < -1000:
        continue

      height = distance2height(dist)

      if atm == atm2:
        continue
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

############################################################################################################################################################
#                                                              Compare and Parse NOESY Peak Lists                                                          #                                                               
############################################################################################################################################################
def reorganize_peaks(user_peaks):
  # clean peaks
  user_peaks = list(filter(lambda x: not np.isnan(x[0]) and not np.isnan(x[1])
                                      and not np.isnan(x[2]), user_peaks))

  w1_shifts = list(map(lambda x: x[0], user_peaks))
  w2_shifts = list(map(lambda x: x[1], user_peaks))
  w3_shifts = list(map(lambda x: x[2], user_peaks))

  w1_avg, w1_std = np.average(w1_shifts), np.std(w1_shifts)
  w2_avg, w2_std = np.average(w2_shifts), np.std(w2_shifts)
  w3_avg, w3_std = np.average(w3_shifts), np.std(w3_shifts)

  # heavy atoms - highest average
  ha_idx = [w1_avg, w2_avg, w3_avg].index(np.max([w1_avg, w2_avg, w3_avg]))

  # through-space hydrogens - higher std. dev. among two remainders
  std_list = [w1_std, w2_std, w3_std]
  del std_list[ha_idx] # without heavy atom
  thru_hidx = std_list.index(np.max(std_list))
  if ha_idx <= thru_hidx:
    thru_hidx += 1

  # one left will be attached hydrogen
  attached_hidx = [0, 1, 2]
  attached_hidx.remove(ha_idx)
  attached_hidx.remove(thru_hidx)
  attached_hidx = attached_hidx[0]

  peaks = list(map(lambda x: [x[ha_idx], x[attached_hidx], x[thru_hidx]],
                   user_peaks))
  return peaks
 
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

def parse_NOESY_list(NOESY_list):
  f = open(NOESY_list, 'r')
  lines = f.readlines()
  f.close()
  nc_shift_list, H_shift_list, h_shift_list = [],[],[]
  peak_height_dict = {}
  for line in lines:
    line = line.strip().split()
    if len(line) < 4:
      continue
    if line[0] == 'Assignment':
      continue
    if line[0] != '?-?-?':
      line[0] = '?-?-?'
    nc_shift = line[1]
    H_shift = line[2]
    h_shift = line[3]
    peak_height= line[4]
    peak_height_dict[float(H_shift), float(nc_shift), float(h_shift)] = peak_height
    nc_shift_list.append(nc_shift)
    H_shift_list.append(H_shift)
    h_shift_list.append(h_shift)
  return nc_shift_list, H_shift_list, h_shift_list, peak_height_dict

############################################################################################################################################################
#                                                              Finding NOESY Candidates                                                                    #                                                               
############################################################################################################################################################

def read_pokybuilder_seq(seq_file):
  f = open(seq_file, 'r')
  lines = f.readlines()
  f.close()

  # start parsing
  seq_dict = {}
  for line in lines:
    sp_list = line.strip().split() #  MET 1
    try:
      aaa, nseq = sp_list
      a = AAA_dict[aaa]
      nseq = int(nseq)
      seq_dict[nseq] = [aaa, a]
    except:
      continue

  nstart = min(list(seq_dict.keys()))
  #print (nstart)
  nend = max(list(seq_dict.keys()))

  return seq_dict, nstart, nend # ({1: ['GLY', 'G'],...},1,114)

def matching_seq_atom(cs, cs_list, tol):
  filtered_cs_list = list(filter(lambda x: abs(x[2] - cs) < tol, cs_list))
  seq_atom_list = []
  for nseq, atm, shift in filtered_cs_list:
    seq_atom_list.append([abs(shift-cs), nseq, atm])
  seq_atom_list.sort()
  return seq_atom_list

def matching_bond(cs, cs_list, tol, cs2, cs_list2, tol2):
  filtered_cs_list = list(filter(lambda x: abs(x[2] - cs) < tol, cs_list))
  filtered_cs_list2 = list(filter(lambda x: abs(x[2] - cs2) < tol2, cs_list2))
  seq_atom_list = []

  for nseq, atm, shift in filtered_cs_list:
    for nseq2, atm2, shift2 in filtered_cs_list2:
      if nseq != nseq2:
        continue
      if atm == 'N' or atm == 'C' and atm2 == 'H':
        seq_atom_list.append([sqrt((abs(shift-cs)/tol)**2 + (abs(shift2-cs2)/tol2)**2),
                              nseq, atm, nseq2, atm2])
        continue
      if len(atm) < 2 or len(atm2) < 2:
        continue
      # HB3 -> B3, CB -> B
      if atm2[1:].find(atm[1:]) == 0:
        seq_atom_list.append([sqrt((abs(shift-cs)/tol)**2 + (abs(shift2-cs2)/tol2)**2),
                      nseq, atm, nseq2, atm2])

  seq_atom_list.sort()
  return seq_atom_list
  
############################################################################################################################################################
#                                                              Calculating Distance Probabilities                                                          #                                                               
############################################################################################################################################################

def dUPL (CaliConstant, peak_intensity, CaliRatio):
  # dUpl = pow( dCaliConst / dPeakIntensity, 1 / 6.0) * dCaliRatio * 1.5
  distance_candidates = pow((CaliConstant/peak_intensity),1/6.0)*float(CaliRatio)*1.5
  return distance_candidates

def distance_probability (pdbList):
  #obtain x,y,z coordinates from alphafold file
  readPDB_ = readPDB(pdbList,1)
  distMat, keyList = createDistanceMatrix(readPDB_) # wil give distMat and keyList
  for i in range(distMat.shape[0]):
    for j in range(distMat.shape[1]):
      for k in range(distMat.shape[2]):
        for l in range(distMat.shape[3]):
          if k == l:
            distMat[i, j, k, l] = 0
            continue
          if distMat[i, j, k, l] == 0:
            continue
          distMat[i, j, k, l] = \
            math.exp((1-max(distMat[i, j, k, l] - 6, 0)**1.5))
  return distMat, keyList

############################################################################################################################################################
#                                                              NOESY Simulation: Write Peak List                                                           #                                                               
############################################################################################################################################################
job_id = '240926_194514_267'
PDB_id = '2N74'
num_oligomer = 2 

os.system(f'mkdir /Users/Karen/AHNA/2N74/2N74_DeepMind/Dimer')

work_directory = f'/Users/Karen/AHNA/2N74/' # set parent directory

# use link to access colab notebook to insert H into AlphaFold PDB
link = f'https://colab.research.google.com/drive/1e8FbwUuS1sdFqk5uj5KVkRAOroIqUOQA#scrollTo=j1b-8LoRjF74'

#AlphaFold PDB Files
original_pdb_filepath = f'{work_directory}/monomer/2N74.pdb' # PDB file from PDB website
AF_pdb_filepath = f'{work_directory}/2N74_DeepMind/content/interface' # AlphaFold with added H

# get chemical shifts
prot_path =f'{work_directory}{job_id}/pokybuilder.prot'
cs_dict, cs_list = read_pokybuilder_prot(prot_path)

# get sequence information
seq_file = f'{work_directory}{job_id}/pokybuilder.seq'
seq_dict, nstart, nend = read_pokybuilder_seq(seq_file)

# NOESY types  
NoesyType = 'cnoe' # or 'nnoe'

regex = r"\d+"
c_list = 'ABCD'

# open original PDB for protein
f= open(original_pdb_filepath, 'r')
lines = f.readlines
f.close()
original_readPDB = readPDB(original_pdb_filepath, 1)

distMat, keyList = createDistanceMatrix(original_readPDB)
original_line = create_noesy_peak_list(keyList, NoesyType, cs_list)

# write out NOESY peak list
PDB_text =open (f'/Users/Karen/AHNA/2N74/2N74_DeepMind/Dimer/original_NOESY_{PDB_id}.txt', 'w')
PDB_text.write(original_line)
PDB_text.close()

# open AF PDB of protein
for PDB in os.listdir(AF_pdb_filepath):
  AF_path = os.path.join(AF_pdb_filepath, PDB)
  if AF_path.endswith('.pdb'):
    AF_readPDB = readPDB(AF_path,1)

    AF_distMat, AF_keyList = createDistanceMatrix(AF_readPDB)
    AF_line = create_noesy_peak_list(AF_keyList, NoesyType, cs_list)
    
    # add chain B to AF structures
    f = open(AF_path, 'r')
    lines = f.readlines()
    f.close()
    # write out NOESY peak list for each AF
    AF_PDB_text = open (f'/Users/Karen/AHNA/2N74/2N74_DeepMind/Dimer/NOESY_{PDB[:-4]}.txt', 'w')
    AF_PDB_text.write(AF_line)
    AF_PDB_text.close()
############################################################################################################################################################
#                                                                   NOESY Comparison                                                                       #                                                               
############################################################################################################################################################
# compare NOESY peak lists
PDB_NOESY_list = f'/Users/Karen/AHNA/2N74/2N74_DeepMind/Dimer/original_NOESY_{PDB_id}.txt'
AF_PDB_list_path = f'/Users/Karen/AHNA/2N74/2N74_DeepMind/Dimer'

for txt in os.listdir(AF_PDB_list_path):
  AF_PDB_list = os.path.join(AF_PDB_list_path, txt) 
  with open(PDB_NOESY_list) as f1, open(AF_PDB_list) as f2:
      # Read data directly into NumPy arrays for efficiency
      PDB = np.genfromtxt(f1, delimiter='', usecols = (1,2,3)) #peak list NOESY
      ALT = np.genfromtxt(f2, delimiter='', usecols=(1,2,3)) # peak list Alt NOESY987
      similarities = compareNOESY(PDB, ALT)
      #print (AF_PDB_list, similarities)

pdbList = f'/Users/Karen/AHNA/2N74/2N74_DeepMind/content/interface/fold_2n74_deepmind_model_0_H.pdb'

# distance probability of AF structure
distMat_probability, keyList_probability = distance_probability(pdbList)
# parse NOESY file
Noesy_file_path = f'/Users/Karen/AHNA/2N74/2N74_DeepMind/Dimer/NOESY_fold_2n74_deepmind_model_0_H.txt'
nc_shift, H_shift, h_shift, peak_height_dict = parse_NOESY_list(Noesy_file_path)
# get chemical shifts for N,C,H; all from *.prot file
N_cs_list = list(filter(lambda x: x[1][0] == 'N', cs_list))
C_cs_list = list(filter(lambda x: x[1][0] == 'C', cs_list))
H_cs_list = list(filter(lambda x: x[1][0] == 'H', cs_list))

# NOE tolerances for N,C,H
tol_dict = {'N': 0.35, 'H': 0.03, 'C': 0.4}

tolNC = tol_dict[NoesyType[0].upper()]
tolH = tol_dict['H']

if NoesyType == 'cnoe':
  NC_cs_list = C_cs_list
else:
  NC_cs_list = N_cs_list

for i in range(len(nc_shift)): #nc_shift is ['124.040', '115.257',...]
  # noesy peak
    NC_cs, H_cs, h_cs = nc_shift[i], H_shift[i], h_shift[i]
    bond_matching = matching_bond(float(NC_cs), NC_cs_list, tolNC, float(H_cs), H_cs_list, tolH)
    h_matching_seq = matching_seq_atom(float(h_cs), H_cs_list, tolH)
    combination = product (bond_matching, h_matching_seq)
    combination = list(combination)

    for j, [NH_cand, h_cand] in enumerate(combination):

      # NH_cand: [distance, nseq, 'N', nseq, 'H']
      # h_cand: [distance, nseq, 'H']
      NH_dist, nseq, NCatom, _, Hatom = NH_cand # distance b/w NH
      h_dist, nseq2, hatom = h_cand # dist b/w NH and h (thru-space)
      A = seq_dict[nseq][1]
      A2 = seq_dict[nseq2][1]
      frq = GetPacsyAbundance(A, Hatom, A2, hatom)
      combination[j] = combination[j] + (frq,)#([1.1644394678885182, 55, 'N', 55, 'H'],
                                              #[0.0009999999999994458, 114, 'H'], 1669)
      Hatom_new = protonNomenclatureChange(A, Hatom)
      hatom_new = protonNomenclatureChange(A2, hatom)

      proton_dist_probability = 0
      for Hatm in Hatom_new:
        #print(Hatm)
        key = f'A_{A}{nseq}{Hatm}'
        k = keyList_probability.index(key)
        for hatm in hatom_new:
          key2 = f'B_{A2}{nseq2}{hatm}'
          l = keyList_probability.index(key2)
          proton_dist_probability = max(proton_dist_probability,
                                        distMat_probability[k,l,0,1])

      combination[j] = combination[j] + (proton_dist_probability,)
    freq_list = list(map(lambda x: x[-2], combination))
    combined_freq = sum(freq_list)
    combination = tuple(map(lambda x: x + (x[-2]/combined_freq,), combination)) #0.4422535211267606

    #Bayes' Theorem?
    combination = tuple(map(lambda x: x + (x[-1]*x[-2],), combination))

    #marginal probability and normalize
    prob_list = list(map(lambda x: x[-1], combination)) #[0.4422535211267606, 0.0, 0.0, 0.0]
    combined_prob = sum(prob_list)
    if combined_prob == 0:
      continue
    combination_ = tuple(map(lambda x:(x[-1]/combined_prob,), combination))
    enduranceScore = combination_*20
    print (enduranceScore) 




# pdb_list= readPDB(pdbList, int(1))

# distMat, keyList = createDistanceMatrix(pdb_list)
# line = ''
# for i in range(len(pdb_list)):
#   nSeq, x, y, z, aaa, a, atm, c = pdb_list[i]
#   #print ('b')
#   if atm[0] not in ['H', 'N', 'O']:
#     #print ('c')
#     continue

#   for j in range(i+1, len(pdb_list)):
#     nSeq2, x2, y2, z2, aaa2, a2, atm2, c2 = pdb_list[j]
#     dist = distance3D( (x, y, z), (x2, y2, z2))
#     #print ('d')
#     if atm2[0] not in ['H', 'N', 'O']:
#       #print ('e')
#       continue

#     if 'H' not in [atm[0], atm2[0]]:
#       #print ('f')
#       continue

#     if distMat[i, j, 0, 1] > float(5.5):
#       #print ("g")
#       continue

#     if 'H' in [atm[0], atm2[0]] and 'O' in [atm[0], atm2[0]]:
#       if distMat[i, j, 0, 1] > 2.2:
#         continue

#     if 'N' in [atm[0], atm2[0]] and 'O' in [atm[0], atm2[0]]:
#       if distMat[i, j, 0 , 1] > 3.3:
#         continue

#     if distMat[i, j, 0, 1] > 0.01 and distMat[i, j, 0, 1] < float(5.5):
#       catm = s32tos12(a, atm)
#       catm2 = s32tos12(a2,atm2)
#       dist = distMat[i, j, 0,1] # giving 'rubber band' 40% more resistance
#       #print (dist)
#       distRound = round(dist, 2)
#       dist2 = round((distRound - 1.7),2) # 1.7A is from VDW

#       line += f'''assign (segid {c} and resid {nSeq} and name {catm}) \
#                   (segid {c2} and resid {nSeq2} and name {catm2}) \
#                   {distRound} {dist2} 0.0 ENDURANCE {enduranceScore}  \n'''
#       #print(line)



#       f = open(f'{work_directory}/alt.tbl', 'w')
#       f.write(line)
#       f.close()
