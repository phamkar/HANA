import os, sys
import itertools
import shutil, glob
# print("Current working directory:", os.getcwd())
parent_dir = '/Users/Karen/AHNA/ScaleFactor'
data_dir = '/Users/Karen/AHNA/230615_175442_137/BestEvaluated'
os.system(f'mkdir {parent_dir}')
def start_value ():
    initial_list = []
    increment = 0.5
    number = 0.00
    while number <= 2.0:
        number += increment
        x = round(number,2)     # only two decimal points
        initial_list.append(x)
    return (initial_list)
#print (start_value())

def stop_value ():
    final_list =[]
    increment = 1.0
    number = 0.0
    
    while number <= 20.0:
        number += increment
        y = round(number,2) # only two decimal points
        final_list.append(y)
    return (final_list)
#print (stop_value())

combined_list = [start_value(), stop_value()]
combination = [p for p in itertools.product(*combined_list)]
#print(type(combination)) # combination is list



for i in combination:   # i is a tuple
    #print (i)
    string_i = str(i)
    new_string_i = string_i.replace(')', ')\n')
    #print (i[0])
    #print(len(new_string_i))
    f = open (f'/Users/Karen/AHNA/230615_175442_137/BestEvaluated/scales_factor.txt', 'w')
    f.write(new_string_i)
    f.close()
    run_dir = os.path.join(parent_dir, f'{i[0]}_{i[1]}')# run_dir = '/somewhere/0.01_0.10'
    
    
    
    #print (run_dir)
    #os.mkdir(run_dir)     
    os.system(f'mkdir {run_dir}')
    os.system(f'cp {data_dir}/*.* {run_dir}/')
    scale_file = os.path.join(run_dir, 'scale.txt')
    f = open(scale_file, 'w')
    f.write(str(i[0])+'\n')
    f.write(str(i[1])+'\n')
    f.close()
    

        
    #for pdb_files in os.listdir(run_dir):
     #   pdb_file_name = ('ahna_dimer_0.pdb')
      #  pdb_file_dir = os.path.join(run_dir, pdb_file_name)
        #print (pdb_file_dir)
   
    '''for files in os.listdir(data_dir):
        source = os.path.join(data_dir, files) 
        destination = os.path.join(run_dir, files) 
        try:
            shutil.copy(source, destination)
            print ('Successful for copying ', files)
        except Exception as e:
            print(f'An error occurred: {e}')
    '''
    
    
    
    # f = open(f'{run_dir}/ahna_dimer.py', 'r')
    # content = f.readlines()
    # f.close()
    # for j in content:
    #     #print("j: ", j)
    #     if j == "start_value":
    #         content = content.replace(j, str(i[0]))
    #     if j == "stop_value":        
    #         content = content.replace (j, str(i[1]))
    
    # f = open(f'{run_dir}/ahna_dimer.py', 'w')
    # for line in content:
    #     f.write(line)
    # f.close()
    
    os.system(f'cd {run_dir}; xplor -py ./ahna_dimer.py')
    