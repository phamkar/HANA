pdb_seqres_url = 'https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt'
pdb_entrytype_url = 'https://files.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt'

import os
if not os.path.exists(pdb_seqres_url):
  os.system('wget ' + pdb_seqres_url)

if not os.path.exists(pdb_entrytype_url):
  os.system('wget ' + pdb_entrytype_url)

f = open('pdb_seqres.txt', 'r')
lines = f.readlines()
f.close()

f = open('pdb_entry_type.txt', 'r')
elines = f.readlines()[1:]
f.close()

edict = {}
for line in elines:
  sp_list = line.strip().split()
  edict[sp_list[0]] = sp_list[1:]
  #print (edict[sp_list[0]])

homodimers = []
for i in range(len(lines)-5):
  line = lines[i]
  if line[0] != '>':
    continue

  # only proteins
  if line.find('mol:protein') == -1:
    continue

  # CHAIN A
  if line[6] != 'A':
    continue

  # CHAIN B
  if lines[i+2][6] != 'B':
    continue

  # HOMODIMER?
  if lines[i+1] != lines[i+3]:
    continue

  if lines[i+4][6] != 'A':
    continue

  sp_list = line.split()
  nseq = int(sp_list[2].replace('length:', ''))

  # protein size between 40 and 400
  if nseq > 400 or nseq < 40:
    continue

  #NMR?
  pdb_id = line[1:5]
  try:
    if edict[pdb_id][1] != 'NMR':
      continue
  except:
    print(pdb_id + ' not exist')
    continue
  # this is proper homodimer by NMR.
  print(line.strip())
  print(lines[i+1].strip())
  homodimers.append([pdb_id, line.strip(), lines[i+1].strip()])
print(str(len(homodimers)) + ' homodimers determined by NMR between 40 and 400 found.')

#@title Creating fasta file

import re

def homodimer_verification(PDB_ID):

#Determine if PDB is Homodimer
  for information in range(len(homodimers)):
    information_PDB = homodimers[information]
    ID = information_PDB[0]
    prot_info = information_PDB[1]
    fastaSeq = information_PDB[2]

    fasta_match = re.search (PDB_ID, prot_info)
    if fasta_match:
      _match_result = (f'>sequence_1\n{fastaSeq}\n>sequence2\n{fastaSeq}')

      fasta_file = open(f'{PDB_ID}.fasta', 'w')
      fasta_file.write(_match_result)
      fasta_file.close()
      return fasta_file

# Establish PDB ID
PDB_ID = '1a03' #@param
homodimer_verification(PDB_ID)

