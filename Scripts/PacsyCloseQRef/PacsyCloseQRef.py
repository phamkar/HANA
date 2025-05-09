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
    return 0;

  f = open(f'pacsymat/pacsyoccurrence_{szA}_all.txt', 'r')
  lines = f.readlines()
  f.close()
  line = lines[iA2*43+iH]
  return int(line.strip().split(',')[iH2])

import sys
ret = GetPacsyAbundance(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
print(ret)