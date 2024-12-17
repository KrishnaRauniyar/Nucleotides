
#Program to search specfic keys consist of thrianges CA, CB and CG atoms.
#Author: Tarikul Islam Milon
#Created on: 06/30/2023

import os
import csv #Changed


# folder path
dir_path = '/ddnB/work/wxx6941/TSR/code/code/Krishna_Code/nt_tsr_key_11428/key_generator/proteins'

OutputFile=open('P_C6_C3_nt_key_extraction_719922263.txt','w')

File=open(f"{dir_path}/sample_details_update.csv",'r')
reader=csv.reader(File)
next(reader)

atom_list=['P','C6','C3\'']
for row in reader:
    FileName = row[0]
    f = open(f'{dir_path}/{FileName}.keys_theta29_dist18', 'r')
    next(f)
    for row in f:
        data = row.split()
        Atom1 = data[1]
        Atom2 = data[3]
        Atom3 = data[5]
        if Atom1 in atom_list and Atom2 in atom_list and Atom3 in atom_list:
            OutputFile.write(f'{FileName}\t{row}')

print('Completed Successfully')
