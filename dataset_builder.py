import sys
import os
import shutil

#
# This script is used to automate the process of extracting
# PDB (Protein Data Bank) files from a list of specified PDB IDs 
# and chain IDs in an input file.
#
# ex. "dataset_builder.py PDB_list_fileName.txt"
#

def remove(string):
    # Funzione per rimuovere gli spazi da una stringa
    return string.replace(".pdb","")

# Authorization command to execute batch_download.sh script 
auth_command = 'chmod +rwx batch_download.sh'
os.system(auth_command)

pdbs_list_filename = sys.argv[1]

# create dataset directory 
name = pdbs_list_filename.split(".")[0]
dataset_dir_name = name + '_dataset'
os.mkdir(dataset_dir_name)

pdb_list = []
with open(pdbs_list_filename, "r") as f:
    pdb_list = f.readlines()

# Executes extractor_pdb.py script and moves extracted chain pdb file
# into dataset directory
for line in pdb_list:
    split_list = line.split()
    pdb_id = remove(split_list[0])
    chain_id = split_list[1]
    command = 'python3 extractor_pdb.py ' + pdb_id + ' ' + chain_id
    os.system(command)
    file_to_move = pdb_id + '_chain_' + chain_id + '.pdb'
    shutil.move(file_to_move, dataset_dir_name + "/" + file_to_move)
