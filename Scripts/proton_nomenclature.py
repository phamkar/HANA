############################################################################################################################################################
#                                                                 Proton Nomenclature                                                                      # 
############################################################################################################################################################
#amino acid sequence
AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

protein_attached_heavy_atoms = {
  'A': {'HN':'N', 'HA':'CA', 'HB1':'CB', 'HB2':'CB', 'HB1':'CB', 'MB':'CB', 'HB':'CB'},
  'C': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB', 'HG':'S'},
  'D': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB', 'HD2':'OD2'},
  'E': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HE2':'OE2', 'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG'},
  'F': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD1':'CD1', 'HD2':'CD2', 'QD':'CQD', 'HD':'CD',
        'HE1':'CE1', 'HE2':'CE2', 'QE':'CQE', 'HE':'CE', 'HZ':'CZ'},
  'G': {'HN':'N', 'HA2':'CA', 'HA1':'CA', 'QA':'CA', 'HA':'CA'},
  'H': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD1':'ND1', 'HD2':'CD2', 'HE1':'CE1', 'HE2':'NE2'},
  'I': {'HN':'N', 'HA':'CA', 'HB':'CB',
        'HD11':'CD1', 'HD12':'CD1', 'HD13':'CD1', 'MD1':'CD1', 'HD1':'CD1', 'HD':'CD1',
        'HG12':'CG1', 'HG11':'CG1', 'QG1':'CG1', 'HG1':'CG1',
        'HG21':'CG2', 'HG22':'CG2', 'HG23':'CG2', 'MG2':'CG2', 'HG2':'CG2', 'HG':'CG2'},
  'K': {'HN':'N', 'HA':'CA',
        'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB', 'HD2':'CD', 'HD1':'CD', 'QD':'CD', 'HD':'CD',
        'QE':'CE','HE':'CE','HE2':'CE', 'HE1':'CE', 'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG',
        'HZ1':'NZ', 'HZ2':'NZ', 'HZ1':'NZ'},
  'L': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD11':'CD1', 'HD12':'CD1', 'HD13':'CD1', 'MD1':'CD1', 'HD1':'CD1',
        'HD21':'CD2', 'HD22':'CD2', 'HD23':'CD2', 'MD2':'CD2', 'HD2':'CD2', 'QMD':'CQD', 'HD':'CD',
        'HG':'CG'},
  'M': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HE1':'CE', 'HE2':'CE', 'HE3':'CE', 'ME':'CE', 'HE':'CE',
        'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG'},
  'N': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD21':'ND2', 'HD22':'ND2', 'QD2':'ND2', 'HD2':'ND2', 'HD':'ND2'},
  'P': {'HT2':'N', 'HT1':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD2':'CD', 'HD1':'CD', 'QD':'CD', 'HD':'CD', 'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG'},
  'Q': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HE21':'NE2', 'HE22':'NE2', 'QE2':'NE2', 'HE2':'NE2', 'HE':'NE2',
        'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG'},
  'R': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD2':'CD', 'HD1':'CD', 'QD':'CD', 'HD':'CD',
        'HE':'NE', 'HG2':'CG', 'HG1':'CG', 'QG':'CG', 'HG':'CG',
        'HH11':'NH1', 'HH12':'NH1', 'QH1':'NH1', 'HH1':'NH1',
        'HH21':'NH2', 'HH22':'NH2', 'QH2':'NH2', 'HH2':'NH2', 'QQH':'NQH', 'HH':'NH'},
  'S': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB', 'HG':'OG'},
  'T': {'HN':'N', 'HA':'CA', 'HB':'CB', 'HG1':'OG1',
         'HG21':'CG2', 'HG22':'CG2', 'HG23':'CG2', 'MG2':'CG2', 'HG2':'CG2', 'HG':'CG2'},
  'V': {'HN':'N', 'HA':'CA', 'HB':'CB',
         'HG11':'CG1', 'HG12':'CG1', 'HG13':'CG1', 'MG1':'CG1', 'HG1':'CG1',
         'HG21':'CG2', 'HG22':'CG2', 'HG23':'CG2', 'MG2':'CG2', 'HG2':'CG2', 'QMG':'CQG', 'HG':'CG'},
  'W': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
         'HD1':'CD1', 'HE1':'NE1', 'HE3':'CE3',
         'HH2':'CH2', 'HZ2':'CZ2', 'HZ3':'CZ3'},
  'Y': {'HN':'N', 'HA':'CA', 'HB2':'CB', 'HB1':'CB', 'QB':'CB', 'HB':'CB',
        'HD1':'CD1', 'HD2':'CD2', 'QD':'CQD', 'HD':'CD',
        'HE1':'CE1', 'HE2':'CE2', 'QE':'CQE', 'HE':'CE', 'HH':'OH'},
}

# change of nomenclature
s12s32 = (
 ('A', 'HN'  ,  'H' ),    
 ('A', 'HB1',  'HB3'),
 ('G', 'HN'  ,  'H' ),
 ('G', 'HA1',  'HA3'),
 ('C', 'HN'  ,  'H' ),
 ('C', 'HB1',  'HB3'),
 ('D', 'HB1',  'HB3'),
 ('D', 'HN' ,   'H' ),
 ('E', 'HN'  ,  'H' ), 
 ('E', 'HB1',  'HB3'),
 ('E', 'HG1',  'HG3'),
 ('F', 'HN'  ,  'H' ),
 ('F', 'HB1',  'HB3'),
 ('H', 'HN'  ,  'H' ), 
 ('H', 'HB1',  'HB3'),
 ('I', 'HN'  ,  'H' ), 
 ('I', 'HG11', 'HG13'),
 ('K', 'HN'  ,  'H' ), 
 ('K', 'HB1',  'HB3'),
 ('K', 'HD1',  'HD3'),
 ('K', 'HG1',  'HG3'),
 ('K', 'HE1',  'HE3'),
 ('K', 'HZ1',  'HZ3'),
 ('L', 'HN'  ,  'H' ), 
 ('L', 'HB1',  'HB3'),
 ('M', 'HN' ,   'H1'), # is this terminal protons for N- terminus
 ('M', 'HN' ,   'H2'), # is this terminal protons for N- terminus
 ('M', 'HN' ,   'H3'), # is this terminal protons for N- terminus
 ('M', 'HB1',  'HB3'),
 ('M', 'HG1',  'HG3'),
 ('M', 'HN'  ,  'H' ),
 ('N', 'HN'  ,  'H' ),
 ('N', 'HB1',  'HB3'),
 ('P', 'HN'  ,  'H' ), 
 ('P', 'HB1',  'HB3'),
 ('P', 'HD1',  'HD3'),
 ('P', 'HG1',  'HG3'),
 ('Q', 'HN'  ,  'H' ),
 ('Q', 'HB1',  'HB3'),
 ('Q', 'HG1',  'HG3'),
 ('R', 'HN'  ,  'H' ), 
 ('R', 'HB1',  'HB3'),
 ('R', 'HD1',  'HD3'),
 ('R', 'HG1',  'HG3'),
 ('S', 'HB1',  'HB3'),
 ('S', 'HN' ,  'H'  ),
 ('T', 'HN'  ,  'H' ),
 ('V', 'HN'  ,  'H' ),
 ('W', 'HN'  ,  'H' ),
 ('W', 'HB1',  'HB3'),
 ('Y', 'HB1',  'HB3')
)

protonNomenclature = {
'A' : {'HB':['HB2', 'HB3'], 'HB2':['HB2',], 'HB3': ['HB3',]},
'C' : {'HB':['HB2', 'HB3'], 'HB2':['HB2',], 'HB3': ['HB3',]},
'D' : {'HB':['HB2', 'HB3'], 'HB2':['HB2',], 'HB3': ['HB3',], 'HD2':['HD2,']},
'E' : {'HB':['HB2','HB3'],'HG':['HG2','HG3'],'HB2':['HB2',],'HB3':['HB3',],
       'HG2':['HG2',],'HG3':['HG3',],'HE2':['HE2',]},
'F' : {'HB':['HB2', 'HB3'],'HD':['HD1','HD2'],'HE':['HE1', 'HE2'],'HB2':['HB2',],
       'HB3': ['HB3',], 'HD1':['HD1',],'HD2':['HD2'],'HE1':['HE1',],
       'HE2':['HE2',]},
'G' : {'HA':['HA2', 'HA3'], 'HA2':['HA2',], 'HA3': ['HA3',]},
'H' : {'HB':['HB2', 'HB3'],'HD':['HD1','HD2'],'HE':['HE1', 'HE2'],'HB2':['HB2',],
       'HB3': ['HB3',], 'HD1':['HD1',],'HD2':['HD2'],'HE1':['HE1',],
       'HE2':['HE2',]},
'I' : {'HG':['HG12','HG12','HG13','HG21','HG22','HG23'],'HG1':['HG12','HG12','HG13'],
       'HG2':['HG21','HG22','HG23'],'HD1':['HD11'], 'HD':['HD11','HD12','HD13'],
       'HG12':['HG12',],'HG13':['HG13',],'HG13':['HG23',],'HG21':['HG21',],
       'HG22':['HG22',],'HG23':['HG23',]},  
'K' : {'HB':['HB2','HB3'],'HG':['HG2','HG3'],'HD':['HD2','HD3'],'HE':['HE2','HE3'],
       'HZ':['HZ1','HZ2','HZ3'],'HB2':['HB2',],'HB3':['HB3',],'HG2':['HG2',],
       'HG3':['HG3',],'HD2':['HD2',],'HD3':['HD3',],'HE2':['HE2',],'HE3':['HE3',],
       'HZ1':['HZ1',],'HZ2':['HZ2',],'HZ3':['HZ3'],},
'L' : {'HB':['HB2', 'HB3'],'HD':['HD11','HD12','HD13','HD21','HD22','HD23'],'HD2':['HD21','HD22','HD23'],
       'HD1':['HD11','HD12','HD13'],'HB2':['HB2',], 'HB3': ['HB3',],'HD11':['HD11',],'HD12':['HD12',],
       'HD13':['HD13',],'HD21':['HD21',],'HD22':['HD22',],'HD23':['HD23',]},
'M' : {'HB':['HB2', 'HB3'],'HG':['HG2', 'HG3'],'HE':['HE1','HE2','HE3'],
       'HB2':['HB2',], 'HB3': ['HB3',],'HG2':['HG2',], 'HG3': ['HG3',],'HE1':['HE1',],
       'HE2':['HE2',],'HE3':['HE3',]},
'N' : {'HB':['HB2', 'HB3'],'HD2':['HD21','HD22'], 'HB2':['HB2',], 'HB3': ['HB3',],
       'HD21':['HD21',],'HD22':['HD22',]},
'P' : {'HB':['HB2', 'HB3'],'HG':['HG2', 'HG3'],'HD':['HD2','HD3'],
       'HB2':['HB2',], 'HB3': ['HB3',],'HG2':['HG2',], 'HG3': ['HG3',],'HD2':['HD2',],
       'HD3':['HD3',]},
'Q' : {'HB':['HB2', 'HB3'],'HG':['HG2', 'HG3'],'HE2':['HE21','HE22'],
       'HB2':['HB2',], 'HB3': ['HB3',],'HG2':['HG2',], 'HG3': ['HG3',],'HE21':['HE21',],
       'HE22':['HE22',]},
'R' : {'HB':['HB2', 'HB3'],'HG':['HG2', 'HG3'],'HD':['HD2','HD3'],'HH':['HH11','HH12','HH21','HH22'],
       'HH1':['HH11','HH12'],'HH2':['HH21','HH22'],'HB2':['HB2',], 'HB3': ['HB3',],
       'HG2':['HG2',], 'HG3': ['HG3',],'HD2':['HD2',],'HD3':['HD3',],'HH11':['HH11'],
       'HH12':['HH12'],'HH21':['HH21'],'HH22':['HH22']},
'S' : {'HB':['HB2', 'HB3'], 'HB2':['HB2',], 'HB3': ['HB3',]},
'T' : {'HG2':['HG21','HG22','HG23'],'HG':['HG1', 'HG21','HG22','HG23'], 'HG1':['HG1',],'HG21': ['HG21',],
       'HG22': ['HG22',],'HG23': ['HG23',]},
'V' : {'HG1':['HG11','HG12','HG13'],'HG':['HG11','HG12','HG13', 'HG21','HG22','HG23'],
       'HG2':['HG21','HG22','HG23'], 'HG11':['HG11',],'HG12':['HG12',],
       'HG13':['HG13',],'HG21':['HG21',],'HG22': ['HG22',],'HG23': ['HG23',]},
'W' : {'HB':['HB2', 'HB3'], 'HB2':['HB2',], 'HB3': ['HB3',]},
'Y' : {'HB':['HB2', 'HB3'],'HD':['HD1', 'HD2'],'HE':['HE1', 'HE2'],'HB2':['HB2',],
       'HB3': ['HB3',], 'HD1':['HD1',],'HD2':['HD2',],'HE1':['HE1',],'HE2':['HE2',]}
}
# create pseudo atoms from PDB to fit XPLOR formating
# [HG11, HG12, HG12] --> [HG*]
# pseudo_atoms = {
# 'A' : {'HB1':'HB*', 'HB2':'HB*', 'HB3':'HB*'},
# 'R' : {'HB1':'HB*', 'HB2':'HB*'}    
# }

a2aaa = (
('A', 'ALA'),
('R', 'ARG'),
('N', 'ASN'),
('D', 'ASP'),
('C', 'CYS'),
('E', 'GLU'),
('Q', 'GLN'),
('G', 'GLY'),
('H', 'HIS'),
('I', 'ILE'),
('L', 'LEU'),
('K', 'LYS'),
('M', 'MET'),
('F', 'PHE'),
('P', 'PRO'),
('S', 'SER'),
('T', 'THR'),
('W', 'TRP'),
('Y', 'TYR'),
('V', 'VAL')
)
import os 
# getting 3 letter aa code
def aaafroma (amino_acid):
       for acid in a2aaa:
              if amino_acid == acid[0]:
                     return acid[1]

def write_seq_file (work_directory,PDB_id):
       # work directory = /Users/Karen/AHNA/{PDB_id}
       seq_content = ''
       for files in os.listdir(f'{work_directory}/monomer'):
            new_path = os.path.join(f'{work_directory}/monomer/{files}')
            if new_path.endswith('.fasta'):          
              f = open(f'{new_path}', 'r')
              lines = f.readlines()
              f.close()
              lines = [row.rstrip('\n') for row in lines][1]
              seq = ''.join(lines[0:])
              for i in range(len(seq)):
                     seq_id = seq[i]
                     # if seq_id.startswith('>'):
                     #        continue
                     aaa = aaafroma(seq_id)
                     seq_content += f'{aaa}  {i+1} \n'
       seq_file = open(f'{work_directory}/monomer/{PDB_id}.seq', 'w')
       seq_file.write(seq_content)
       seq_file.close()
       return seq_file

# changing from 2/3 nomeclature to 1/2
def s32tos12(szA, szHX):
  for grp in s12s32:
    if (grp[0] == szA) and (grp[2] == szHX):
      return grp[1]
  return szHX



