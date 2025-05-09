import os, cmd

pdb_ref = "/Users/Karen/ahna/run/BestEvaluated/2n74_demo/test/2n74_ref.pdb"
ahna = '/Users/Karen/ahna/run/BestEvaluated/2n74_demo/test/' 

def automated_files (ahna):
    files_list = os.listdir(ahna)
    all_files = list()
    for subdir in files_list:
        f1 = os.path.join(subdir, 'ahna_dimer_0.pdb')
        f2 = os.path.join(subdir, 'ahna_dimer_1.pdb')
        all_files.append(f1)
        all_files.append(f2)
    return all_files
#print(automated_files(ahna))  # <Dir Path>  ==> your directory path

cmd.load(pdb_ref, 'ref')
all_files = automated_files(ahna)
cmd.load(automated_files(ahna), 'ahna')

ahna_report = ''
for file in all_files:
    cmd.load(file, 'ahna')
    ret = cmd.align('ref & n. N+CA+C', 'ahna', cutoff=0, mobile_state=1, target_state=1, quiet=0)
    RMSD = ret[0]
    line = f'{file}: {RMSD}'
    ahna_report.append(line)
    #print(line)
    #print(RMSD)



        



