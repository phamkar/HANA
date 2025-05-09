''' 
___________________________________________________________________________________________________________________________________________________________
Benchmark 2.0
This script is for benchmarking homodimer proteins. The inputs are chemical shifts, chemical sequence, and AF protein PDB files. From the PDB files will be 
used to generate a NOESY peak list, where all assignments are assumed to be unknown and candidates are determined. From there NOESY peak list can be used
to create restraints and be cycled through PONDEROSA. 
___________________________________________________________________________________________________________________________________________________________
'''

import os
import numpy as np
import re
import math
from math import sqrt
from itertools import product
import proton_nomenclature 
import PacsyCloseQRef as PACSY
import json
import pynmrstar

############################################################################################################################################################
#                                                                 Parsing Chemical Shifts                                                                  # 
############################################################################################################################################################

# for benchmarking, the chemical shifts will be from BMRB database instead of the .prot file from POKY Builder. 
# BMRB contains pseudo atoms 

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
  # sequence list includes a set of lists
    # [nseq, 3-letter-seq, 1-letter-seq]
  sequence_loops = strentry.get_loops_by_category('_Entity_comp_index')
  sequence_dict = json.loads(sequence_loops[0].get_json())
  sequence_list = list(map(lambda x: [x[1], x[2], proton_nomenclature.AAA_dict[x[2]]],
                           sequence_dict['data']))
  fasta = ''.join(list(map(lambda x: x[2], sequence_list)))

  # extract chemical shifts
  atom_chem_shift_loops = strentry.get_loops_by_category('atom_chem_shift')

  # First, try from bond and chemical shift information
  # shift_list includes a set of lists:
  #     [nseq, 3-letter-seq, 1-letter-seq,atom, chemical shift]
  atom_chem_shift_dict = json.loads(atom_chem_shift_loops[0].get_json())
  atom_chem_shift_data = atom_chem_shift_dict['data']
  
  new_atom_chem_shift_data = []
  for item_ls in atom_chem_shift_data:
    new_item_ls =[]
    for i, item in enumerate(item_ls):
      if item != '.':
         new_item_ls.append(item)
    new_atom_chem_shift_data.append(new_item_ls)

  shift_list_a = []
  shift_list_b = []
  for item_ls in new_atom_chem_shift_data:
    line_ls_a = [int(item_ls[4]), item_ls[6], float(item_ls[9])]
    # line_ls_b = [int(item_ls[4])+10000, item_ls[6], item_ls[9]]
    shift_list_a.append(line_ls_a)
    # shift_list_b.append(line_ls_b)
      #if re.fullmatch(r"[a-zA-Z]+", item):
      # print (i,item)
  # bmrb_nseq_a = item_ls[i-1]
  # bmrb_nseq_b = int(item_ls[i-1])+10000
  # bmrb_aaa = item
  # bmrb_atom = item_ls[i+1]
  # bmrb_cs = item_ls[i+3]
  # shift_list_a.append([bmrb_nseq_a, bmrb_aaa, bmrb_atom, bmrb_cs])
  # shift_list_b.append([bmrb_nseq_b, bmrb_aaa, bmrb_atom, bmrb_cs])
  #shift_list_a = list(map(lambda x: [x[4], x[7], x[10]],
  #                             atom_chem_shift_dict['data']))

  #shift_list_b = list(map(lambda x: [str(int(x[4])+10000), x[7], x[10]],atom_chem_shift_dict['data']))
  # shift_list = shift_list_a+shift_list_b
  shift_list = shift_list_a
  return shift_list # [nseq, atom, chemical shift]

def filter_shift_list(cs_list):
    """ 
    ____________________________________________________________________
    Finds identical shifts and modifies the atm name   
    input: chemical shift list from BMRB
    output: new chemical shift list with pseudo atoms identified with *
    ____________________________________________________________________
    """
    for i in range(len(cs_list)):
        nseq, atm, shift = cs_list[i]  # Correct way to access elements
        if atm == 'H':
           cs_list[i][1] = atm.replace('H', 'HN') # change nomenclature and update list
        for j in range(i + 1, len(cs_list)):  # Start from i+1 to avoid comparing with itself
            nseq2, atm2, shift2 = cs_list[j]
            if atm2 == 'H':
               cs_list[j][1] = atm2.replace('H', 'HN')
            if shift == shift2:
                # Modify atm based on the comparison
                match = re.search(r"(\D+)\d+$", atm)
                if match:
                    cs_list[i][1] = re.sub(r'\d+$', '*', atm)

                match2 = re.search(r"(\D+)\d+$", atm2)
                if match2:
                    cs_list[j][1] = re.sub(r'\d+$', '*', atm2)
    seen = set()
    filtered_cs_list = []
    for sublist in cs_list:
        tuple_sublist = tuple(sublist)  # Convert to tuple
        if tuple_sublist not in seen:
            seen.add(tuple_sublist)
            filtered_cs_list.append(sublist)
    return filtered_cs_list

############################################################################################################################################################
#                                                                    Parsing PDB                                                                           #                                                   #
############################################################################################################################################################

def readPDB(pdb_file, modelnumber, nuclei=['H',]):
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
    if line[0:4] in ['MODE', 'ATOM', 'ENDM']:
      pdbLines.append(line)


  # fill modelList
  for line in pdbLines:
    if len(line[12:16]) == 0:
       continue
    if line[12:16].strip()[0] not in nuclei:
      continue
    if line[0:5] == 'MODEL':
      modelList.append([])
    if line[0:4] != 'ATOM':
      continue

    c = line[21:22].strip() # <chain_ID>
    aaa = line[17:20].strip()
    atm = line[12:16].strip()
    nSeq = int(line[23:26].strip()) # <residue sequence>
    x = float(line[30:38].strip())
    y = float(line[38:46].strip())
    z = float(line[46:54].strip())
     
    
    # define and add psuedo atoms
    if atm == 'H':
        atm = atm.replace('H','HN')
    if len(atm)> 3:
        match = re.search(r"(\D+)\d+$", atm)
        if match: 
          atm = re.sub(r'\d+$', '*', atm)
    # # in case MODEL not in PDB.
    
    if len(modelList) == 0:
      modelList.append ([])
    
    # if c == 'B':
    #   nSeq = nSeq - 200
    modelList[-1].append( [nSeq, x, y, z, aaa, proton_nomenclature.AAA_dict[aaa], atm, c] )
    if len(modelList) == 0:
       continue
     
  return modelList[-1]
  # print (modelList[modelnumber-1])
  # return modelList[0]
############################################################################################################################################################
#                                                                 Proton Nomenclature                                                                      # 
############################################################################################################################################################
def protonNomenclatureChange (aa, atom):
  try:
    prot_aa = proton_nomenclature.protonNomenclature[aa] # 'V' -> {'HG':['HG11','HG12','HG13', 'HG21','HG22','HG23']
    nomenclature = prot_aa[atom] # 'HG' -> ['HG11','HG12','HG13', 'HG21','HG22','HG23']
  except:
    return [atom,]
  return nomenclature

############################################################################################################################################################
#                                                                 Read Sequence File                                                                       # 
############################################################################################################################################################

def read_seq_file(seq_file):
    ''' 
    ______________________________________________________________________
    input: protein sequence file path
    output: sequence dictionary, first residue number, last residue number
    key = residue number/ sequence id
    value = 3 letter code amino acid, one letter code amino acid 
    ({1: ['GLY', 'G'],...},1,114)
    ______________________________________________________________________
    '''
    f = open(seq_file, 'r')
    lines = f.readlines()
    f.close()
    # start parsing
    seq_dict = {}
    for line in lines:
        sp_list = line.strip().split() #  MET 1
        try:
            aaa, nseq = sp_list
            a = proton_nomenclature.AAA_dict[aaa]
            nseq = int(nseq)
            seq_dict[nseq] = [aaa, a]
            seq_dict[nseq+10000] = [aaa,a]
        except:
            continue

    nstart = min(list(seq_dict.keys()))
    nend = max(list(seq_dict.keys()))

    return seq_dict, nstart, nend # ({1: ['GLY', 'G'],...},1,114)

############################################################################################################################################################
#                                                     NOESY Simulation: Setting Simulated NOESY Peak List                                                  #                                                               
############################################################################################################################################################

def distance3D(atom1, atom2):
    """ 
    _________________________________________________________
    takes two coordinates. ex: ((26.266, 25.413, 2.842),
                             (26.913, 26.639, -3.51))
    returns the distance

    Euclidean calculation
    _________________________________________________________
    """
    return sqrt((atom1[0] - atom2[0]) ** 2 +
              (atom1[1] - atom2[1]) ** 2 +
              (atom1[2] - atom2[2]) ** 2)

def createDistanceMatrix(pdb_list):
    '''
    _________________________________________________________________________
    input: parsed PDB list
    output: distance matrix and key list

    calculating the distance between each and every atom

    setting the 'keylist': chainID + '_'+ amino acid + sequence number + atom
    _________________________________________________________________________
    '''
    distMat = np.zeros( (len(pdb_list), len(pdb_list), 4, 4 ) )
    keyList = []
    for i in range(len(pdb_list)):
      nSeq, x, y, z, aaa, a, atm, c = pdb_list[i]
      keyList.append(c+'_'+a+str(nSeq)+atm)
      for j in range(i+1, len(pdb_list)):
        nSeq2, x2, y2, z2, aaa2, a2, atm2, c2 = pdb_list[j]
        dist = distance3D( (x, y, z), (x2, y2, z2))
        distMat[i, j, 0, 1] = \
        distMat[j, i, 1, 0] = dist
    return distMat, keyList

def get_shift(seqidx, atomname, cs_list):
    '''
    ____________________________________________________________
    input: sequence id/residue number, atom, and chemical shifts
    output: chemical associated to the sequence id and atom
    ____________________________________________________________
    '''

    filterAtoms = list(filter(lambda x :x[0] == seqidx and x[1] ==atomname,
                            cs_list))
    
    if len(filterAtoms) == 0:
        return -9999
    
    return float(filterAtoms[0][2])

def distance2height(HHdist, offset = 0.0):
    """ 
    __________________________________________________
    input: distance between protons 
    output: NOESY peak height

    distance scaling using r^-6 approximation for H-H 
    __________________________________________________
    """

    min_hh = 1.7  # closest distance between H atoms (1.70 A)
    max_hh = 5.5  # farthest distance observed between H atoms

    max_I = 10**6 # arbitrary value
    min_I = 10**4 # arbitrary value

    approx = -6.0

    dist = min(HHdist, max_hh-offset)
    dist = max(dist, min_hh-offset)

    A = (max_I - min_I) / ((min_hh-offset)**approx - (max_hh-offset)**approx)        # 75227372.04997641
    B = min_I - A * (max_hh-offset)**approx                                          # -19444.257193289242

    return A * dist**approx + B

def parse_key(key):
    '''
    ____________________________________________________________________
    input: keys in KeyList
    output: breakdown of keys: chainID, amino acid, residue number, atom
    ____________________________________________________________________
    '''
    c = key.split('_')[0]
    a = key.split('_')[1][0]
    nSeq = int(re.search(regex, key.split('_')[1]).group())
    sep = c+'_'+a+str(nSeq)
    atm = key.split(sep)[1]
    return c, a, nSeq, atm
def create_noesy_peak_list(pdb_list, noesy_type, cs_list, distMat):
  c_list = 'ABCD'
  lines = ''
  header = '  %20s    %5s   %5s   %5s' % ('Assignments','w1','w2','w3')
  lines += header + '\n'
  lines += '\n'
  
  for i in range(len(pdb_list)):
    nseq, x, y, z, aaa, a, atm, c = pdb_list[i]
    atm = proton_nomenclature.s32tos12(a,atm) # alter nomenclature from 32 system to 12 system for XPLOR-NIH 

    try: # get correlated heavy atom
      ncatm = proton_nomenclature.protein_attached_heavy_atoms[a][atm] 
    except:
      try:
        ncatm = proton_nomenclature.protein_attached_heavy_atoms[a][atm[:2]]
      except:
         continue
    
    if noesy_type == 'nnoe':
      if ncatm[0] != 'N':
        continue
    elif noesy_type == 'cnoe':
      if ncatm[0] != 'C':
        continue

    try: # get shifts from shifts list
      nc_shift = get_shift(nseq, ncatm, cs_list)
    except:
      continue
    try:
      h_shift = get_shift(nseq, atm, cs_list)
    except:
      continue
    
    if nc_shift < -1000 or h_shift < -1000: # if shifts are not present then remove
      continue
    
    grp1 = f'{a}{nseq}{ncatm}' # heavy atom and reside/sequence number
    grp2 = f'{atm}' # connected H

    for j in range(i+1, len(pdb_list)):
      nseq2, x2, y2, z2, aaa2, a2, atm2, c2 = pdb_list[j]
      dist = distMat[i, j, c_list.index(c), c_list.index(c2)]
      if c == c2:
         continue
      if dist == 0:
         continue
      if dist > 6:
         continue
    
      try: # get thru space H shifts from chain B
        h2_shift = get_shift(nseq2, atm2, cs_list)
      except:
        continue
      if h2_shift < -1000:
        continue   
      
      grp3 = f'{atm2}' 
      asgn = f'{grp1}-{grp2}-{grp3}'
      line = '  %20s %8.3f %8.3f %8.3f' % (asgn, nc_shift, h_shift,
                                                h2_shift)
      lines += line + '\n'
      
  lines_list = lines.split('\n')
  seen = set()
  result = []
  for item in lines_list:
     if item not in seen:
        result.append(item)
        seen.add(item)

  return ('\n').join(result)
############################################################################################################################################################
#                                                             Parse and Filter NOESY Peak Lists                                                            #           
############################################################################################################################################################

def parse_NOESY_list(NOESY_list):
  '''
  ______________________________________________________________
  input: simulated NOESY peak list
  output: isolated chemical shifts
  nc = heavy atoms, H = attached protons, h = thru space protons
  ______________________________________________________________
  '''
  f = open(NOESY_list, 'r')
  lines = f.readlines()
  f.close()
  nc_shift_list, H_shift_list, h_shift_list = [],[],[]
  for line in lines:
    line = line.strip().split()
    if len(line) < 4:
      continue
    if line[0] == 'Assignments':
      continue
    nc_shift = line[1]
    H_shift = line[2]
    h_shift = line[3]
    nc_shift_list.append(nc_shift)
    H_shift_list.append(H_shift)
    h_shift_list.append(h_shift)
  return nc_shift_list, H_shift_list, h_shift_list

############################################################################################################################################################
#                                                              Setting NOESY Candidates                                                                    #                                                               
############################################################################################################################################################

def matching_bond(cs, cs_list, tol, cs2, cs_list2, tol2):
    '''
    _______________________________________________________________________________________________________________________________________________________
    input: chemical shifts (cs) from NOESY simulated peak list,
           chemical shifts (cs_list) from BMRB,
           tolerances (tol) set for each atom type
    output: list made up of normalized differences between heavy atoms and protons, 
            sequence number and atom of heavy atom,
            sequence number and atom of corresponding proton
    This function takes into the consideration of chemical shifts of a heavy atom and 'H'. It matched the heavy atoms to a proton, assuming it's the
    attached proton to the heavy atom.
    _______________________________________________________________________________________________________________________________________________________

    '''
    
    filtered_cs_list = list(filter(lambda x: abs(float(x[2]) - cs) < float(tol), cs_list)) # finding proper chemical shifts of HA
    filtered_cs_list2 = list(filter(lambda x: abs(float(x[2]) - cs2) < float(tol2), cs_list2)) # finding proper chemical shifts of protons/'H'
    seq_atom_list = []
    # print (filtered_cs_list)
    # print (filtered_cs_list2)
    for nseq, atm, shift in filtered_cs_list: # for nseq and atom of the filtered list of HA
        shift = float(shift)
        for nseq2, atm2, shift2 in filtered_cs_list2: # for nseq and atom of the filtered list of 'H'
            shift2 = float(shift2)
            if nseq != nseq2: # avoid sequential data
                continue
            if atm == 'N' or atm == 'C' and atm2 == 'H': # appending only N or C, and H
                seq_atom_list.append([sqrt((abs(shift-cs)/tol)**2 + (abs(shift2-cs2)/tol2)**2),
                                    int(nseq), atm, int(nseq2), atm2])
                continue
            if len(atm) < 2 or len(atm2) < 2: # if length of atom is less than two i.e. "H" or "C" then bypass
                continue

            # HB3 -> B3, CB -> B
            if atm2[1:].find(atm[1:]) == 0: # if it is a match
                seq_atom_list.append([sqrt((abs(shift-cs)/tol)**2 + (abs(shift2-cs2)/tol2)**2),
                            int(nseq), atm, int(nseq2), atm2])
    seq_atom_list.sort()
    return seq_atom_list

def matching_seq_atom(cs, cs_list, tol):
  '''
    _______________________________________________________________________________________________________________________________________________________
    input: chemical shifts (cs) from NOESY simulated peak list,
           chemical shifts (cs_list) from BMRB,
           tolerances (tol) set for each atom type
    output: list made up off difference between cs_list and cs,
            sequence number from cs_list,
            atom type from cs_list
    This function matches the thru space 'H' by comparing the chemical shifts within the NOESY peak list to the ones from BMRB database. If the difference 
    between cs and cs_list shifts are below the set tolerance values then it'll be taken into consideration when assigning possible NOESY peak candidates.
    _______________________________________________________________________________________________________________________________________________________
  '''
  filtered_cs_list = list(filter(lambda x: abs(float(x[2]) - cs) < tol, cs_list)) 
  seq_atom_list = []
  for nseq, atm, shift in filtered_cs_list:
    shift = float(shift)
    seq_atom_list.append([abs(shift-cs), int(nseq), atm])
  seq_atom_list.sort()
  return seq_atom_list


############################################################################################################################################################
#                                                             Filtering through NOESY Candidates                                                           # 
############################################################################################################################################################

def aa_contact_map(pdbList):
    '''
    ________________________________________________________________________________________________________________________________________________

    input: parsed PDB file
    output: a matrix that indexes the wanted heavy atom
    This function creates a contact map for indexing the wanted HA to aid in filtering out the combination list, so only the wanted atom is present.
    I think???
    ________________________________________________________________________________________________________________________________________________
    '''
    # create contact map - indexing
    readPDB_list = readPDB(pdbList, 1, nuclei=['C',])
    aa_matrix = np.zeros((len(readPDB_list),len(readPDB_list),4,4)) # for now let's calculate for CA
    for i in range(len(readPDB_list)):
        nSeq, x, y, z, aaa, a, atm, c = readPDB_list[i] 
        if int(nSeq) > 10000:
           nSeq = int(nSeq) - 10000  
        if atm != 'CA':
            continue
        for j in range(i+1, len(readPDB_list)):
            nSeq2, x2, y2, z2, aaa2, a2, atm2, c2 = readPDB_list[j]
            if int(nSeq2) > 10000:
              nSeq2 = int(nSeq2) - 10000
            if atm2 != 'CA':
                continue
            dist = distance3D( (x, y, z), (x2, y2, z2))
            if dist < 15:
                aa_matrix[nSeq - 1, nSeq2 - 1,0, 1] = \
                    aa_matrix[nSeq2 - 1, nSeq - 1, 1, 0] = 1
    return aa_matrix

def distance_probability (pdbList):
    '''
    __________________________________________________________________________________
    input: parsed PDB file from Alphafold-Multimer
    output: the probability of NOESY signal appear, essentially any atoms less than 6A
    __________________________________________________________________________________
    '''
    #obtain x,y,z coordinates from alphafold file
    readPDB_ = readPDB(pdbList,1,nuclei=['H',])
    distMat, keyList = createDistanceMatrix(readPDB_) # wil give distMat and keyList
    for i in range(distMat.shape[0]):
        for j in range(distMat.shape[1]):
            for k in range(distMat.shape[2]):
                for l in range(distMat.shape[3]):
                    if k == l and abs(i-j) < 2:
                        distMat[i, j, k, l] = 1
                        continue
                    elif k == l:
                        distMat[i, j, k, l] = 0
                        continue
                    else:
                       distMat[i, j, k, l] = \
                        math.exp((1-max(distMat[i, j, k, l]-6, 0)**1.5)) # the smaller the value, the less likely the signal would appear
    return distMat, keyList

############################################################################################################################################################
#                                                              Write .tbl File for XPLOR_NIH                                                               # 
############################################################################################################################################################

def aa_dictionary(pdbList):
  readPDB_list = readPDB(pdbList, 1, nuclei = ['H',])
  aa_dict = {}
  for i in range(len(readPDB_list)):
      nSeq, x, y, z, aaa, a, atm, c = readPDB_list[i]
      aa_dict[(atm,nSeq,c)]=(x, y, z, aaa, a)
      for j in range(i+1, len(readPDB_list)):
        nSeq2, x2, y2, z2, aaa2, a2, atm2, c2 = readPDB_list[j]
        if (atm2, nSeq2, c2) not in aa_dict:
          aa_dict[(atm2, nSeq2, c2)] = (x2, y2, z2, aaa2, a2) 
  return aa_dict 
def write_tbl_file(combination, AF_readPDB, aa_dict,seq_dict):
# protons to protons only
# 1 MET    N   2 ASP    H 5.317
# nSeq aaa atom nSeq2 aaa2 atom2 distance
# combination is a tuple
# [HA-H],[H], frq, proton_dist_probability, normailzed frq, Bayes' Theorem, normalized Bayes', endurance score
  proton_matrix = np.zeros((len(AF_readPDB),len(AF_readPDB),4,4)) # (87,87,4,4) for 2N74
  line = ''
  for k in combination:
    # get get attached 'H'and seq_id from chain A
    catm, cnSeq = k[0][4], k[0][1]
    # get amino acid
    aa = seq_dict[int(cnSeq)][1]
    # proton nomenclature change 32 to 12
    catm = proton_nomenclature.s32tos12(aa, catm)
    # get through space 'H' and seq_id from chain B
    catm2,cnSeq2 = k[1][2], k[1][1]
    # get amino acid
    aa2 = seq_dict[int(cnSeq2)][1]
    catm2 = proton_nomenclature.s32tos12(aa2, catm2)
 
    
    # get endurance score
    endurance_score = round(k[-1],2)
    if endurance_score == 0:
       continue
    try:
      x,y,z = aa_dict[catm, int(cnSeq), 'A'][0:3]
    except:
      continue
    try:  
      x2,y2,z2 = aa_dict[catm2,int(cnSeq2),'B'][0:3]
    except:
      continue
    
    dist = distance3D( (x, y, z), (x2, y2, z2))
    # print (dist)
    if dist < 1.7:
      continue
    if dist < 5.5:
      proton_matrix[int(cnSeq)-1,int(cnSeq2)-1,0, 1] = \
        proton_matrix[int(cnSeq2)-1,int(cnSeq)-1, 1, 0] = dist
      distRound = round(dist,2)
      dist2 = round(distRound - 1.7,2)

      line += f'assign (segid A and resid {cnSeq} and name {catm}) \
                  (segid B and resid {cnSeq2} and name {catm2}) \
                  {distRound} {dist2} 0.0 \
                  !ENDURANCE {endurance_score} \n'  
        
  return line

############################################################################################################################################################
#                                                                       Implementation                                                                     # 
############################################################################################################################################################

PDB_id = '5HUZ' 
POKY_job_id = '250429_145141_133'

work_directory = f'/Users/Karen/AHNA/{PDB_id}' # set parent directory
os.system(f"mkdir {work_directory}/Dimer") # to hold NOESY lists and other generated files

# use link to access colab notebook to insert "H" into AlphaFold PDB
link = f'https://colab.research.google.com/drive/1e8FbwUuS1sdFqk5uj5KVkRAOroIqUOQA#scrollTo=j1b-8LoRjF74'

#AlphaFold PDB files with protons added
AF_pdb_filepath = f'{work_directory}/5HUZ_Boltz1' # PDB with protons added

# get chemical shifts
cs_list = read_bmrb(ename=30005)  # [[nseq, atom, chemical shift], ] 
cs_list = filter_shift_list(cs_list) # assign pseudo atoms with *
# get sequence information for homodimers; fasta file will be provided by user
proton_nomenclature.write_seq_file(work_directory, PDB_id)
seq_path = f'{work_directory}/monomer/{PDB_id}.seq'
seq_dict, nstart, nend = read_seq_file(seq_path)
NoesyType = ['cnoe', 'nnoe']

regex = r"\d+"
c_list = 'ABCD'
content = ''

for types in NoesyType:
# open AF PDB of protein with the protons added
  for PDB in os.listdir(AF_pdb_filepath):
    AF_path = os.path.join(AF_pdb_filepath, PDB)
    if AF_path.endswith('H.pdb'):
      AF_readPDB = readPDB(AF_path,1,nuclei = ['H',])
      AF_distMat, AF_keyList = createDistanceMatrix(AF_readPDB)
      AF_line = create_noesy_peak_list(AF_readPDB, types, cs_list,AF_distMat)
      # write out NOESY peak list for each AF structure
      AF_PDB_text = open (f'{work_directory}/Dimer/{types[0].upper()}_NOESY_{PDB_id}.list', 'w')
      AF_PDB_text.write(AF_line)
      AF_PDB_text.close()
  # chosen PDB for analysis
  pdbList = f'{AF_pdb_filepath}/{PDB_id}_H.pdb'
  # print (readPDB(pdbList,1))
  Noesy_file_path = f'{work_directory}/Dimer/{types[0].upper()}_NOESY_{PDB_id}.list'
  # TODO: change this to user input in colab
  nc_shift, H_shift, h_shift = parse_NOESY_list(Noesy_file_path)

  # get chemical shifts for N,C,H; from BMRB
  N_cs_list = list(filter(lambda x: x[1][0] == 'N', cs_list))
  C_cs_list = list(filter(lambda x: x[1][0] == 'C', cs_list))
  H_cs_list = list(filter(lambda x: x[1][0] == 'H', cs_list))
  
  if types == 'cnoe': # or 'nnoe'
    NC_cs_list = C_cs_list
  else:
    NC_cs_list = N_cs_list
  # NOE tolerances for N,C,H
  tol_dict = {'N': 0.35, 'H': 0.03, 'C': 0.4}
  tolNC = tol_dict[types[0].upper()]
  tolH = tol_dict['H']

  ca_contact_matrix = aa_contact_map(pdbList)
  distMat_probability, keyList_probability = distance_probability(pdbList)
  # print (keyList_probability)
  aa_dict = aa_dictionary(pdbList)
  # print (aa_dict)
  for i in range(len(nc_shift)): #nc_shift is ['124.040', '115.257',...]
      NC_cs, H_cs, h_cs = nc_shift[i], H_shift[i], h_shift[i]
      bond_matching = matching_bond(float(NC_cs), NC_cs_list, tolNC, float(H_cs), H_cs_list, tolH) # matching HA to H
      h_matching_seq = matching_seq_atom(float(h_cs), H_cs_list, tolH) # finding thru space h
      combination = product (bond_matching, h_matching_seq) # combination of NH and h
      combination = list(combination)
      for j, [NH_cand, h_cand] in enumerate(combination):
      # NH_cand: [distance, nseq, 'N', nseq, 'H']
      # h_cand: [distance, nseq, 'H']
        NH_dist, nseq, NCatom, _, Hatom = NH_cand # distance b/w NH
        h_dist, nseq2, hatom = h_cand # dist b/w NH and h (thru-space)
        if ca_contact_matrix[int(nseq)-1][int(nseq2)-1][0][1] != 1:
            continue
        A = seq_dict[int(nseq)][1] # get 1 ltr amino acid seq
        A2 = seq_dict[int(nseq2)][1] # get 1 ltr amino acid seq
        frq = PACSY.GetPacsyAbundance(A, Hatom, A2, hatom)  # calculate proton abundance 
        combination[j] = combination[j] + (frq,)    #([1.1644394678885182, 55, 'N', 55, 'H'],[0.0009999999999994458, 114, 'H'], 1669)
        
        try:
            Hatom_new = protonNomenclatureChange(A, Hatom)
            if Hatom.endswith('*'):
                Hatom_new = protonNomenclatureChange(A, Hatom[:2]) # assign protons to psuedo atoms HB* --> ['HB2', 'HB3']
        except:
            continue
        try:    
            hatom_new = protonNomenclatureChange(A2, hatom)
            if hatom.endswith('*'):
                hatom_new = protonNomenclatureChange(A, hatom[:2])
        except:
            continue
        proton_dist_probability = 0
        for Hatm in Hatom_new:
          key = f'A_{A}{nseq}{Hatm}'
          try:
            k = keyList_probability.index(key) # keyList = c+'_'+a+str(nSeq)+atm
          except:
            continue
          for hatm in hatom_new:
            key2 = f'B_{A2}{nseq2}{hatm}'
            try:
              l = keyList_probability.index(key2)
            except:
              continue
            proton_dist_probability = max(proton_dist_probability,
                                        distMat_probability[k,l,0,1])         
        combination[j] = combination[j] + (proton_dist_probability,) # update combination list
          
      # filter combination to have all four elements
      combination = list(filter(lambda x: len(x) == 4, combination))

      # get frq values
      freq_list = list(map(lambda x: x[-2], combination))
      # combine frq values
      combined_freq = sum(freq_list)
      if combined_freq == 0:
        continue
      # normalize frq values to 1
      combination = tuple(map(lambda x: x + (x[-2]/combined_freq,), combination)) #0.4422535211267606
      
      #Bayes' Theorem? Total probability? 
      combination = tuple(map(lambda x: x + (x[-1]*x[-2],), combination))
      
      #marginal probability and normalize
      prob_list = list(map(lambda x: x[-1], combination)) #[0.4422535211267606, 0.0, 0.0, 0.0]
      combined_prob = sum(prob_list)

      if combined_prob == 0:
        continue
      combination = tuple(map(lambda x:x + (x[-1]/combined_prob,), combination))
      # set endurance score to weight of 20: 10.3 + 9.7 = 20
      combination = tuple(map(lambda x: x + (x[-1]*20.0,),combination))
      tbl_line= write_tbl_file(combination, AF_readPDB, aa_dict,seq_dict)
      content += tbl_line
  # print (types, content)
  content_list = content.split('\n') # take content lines and turn into list
  filter_content = list(set(content_list)) # find duplicate lists and remove them
  new_content = '\n'.join(filter_content) # take filtered list and join to strings

  print(new_content)

  f = open(f'{work_directory}/{POKY_job_id}/BestEvaluated/alt.tbl', 'w')
  f.write(new_content)
  f.close()  
