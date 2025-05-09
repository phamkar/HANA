#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 13:40:37 2023

@author: Karen
"""

import os
import pandas as pd

main_dir = '/Users/Karen/AHNA/2N74/240926_194514_267/BestEvaluated/XplorCalc/'
sub_dir = os.listdir(main_dir)

stat_files = []

stat_path = []
#print (sub_dir)

for folders in sub_dir:
    if folders.startswith('0'):
        stat_path.append(main_dir + '/' + folders + '/' )

for i in stat_path:
    stat_path_list = os.listdir(i)
    for j in stat_path_list:
       if j.endswith('.pdb.stats'):
            final_file = pd.read_csv (i+j, sep=',', header = [0])
            print (final_file)

