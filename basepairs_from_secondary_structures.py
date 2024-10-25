# Formats the PDB's 2D structures for all PDBs in the given dataset.
# Generates .pairs (2D structure file format required by ARTS tool) files for all PDBs in the dataset.
# 
# USAGE: python3 basepairs_from_secondary_structures.py <dataset_name> 

import sys
import os
import subprocess
import shutil

dataset_name = sys.argv[1]

dataset_pdbs = os.listdir("./Datasets/" + dataset_name)

for pdb in dataset_pdbs:
    if(not pdb.endswith(".pdb_new")):
        pdb_id = pdb[0:4]
        pdb_name = pdb[0:len(pdb) - 4]
        secondary_structure_path = "./2D_structures/" + dataset_name + "/" + pdb_name + "/" + pdb + ".out"
        subprocess.run(["python3", "secondary_structure_formatter.py", secondary_structure_path])
        basepairs_file_path = "./2D_structures/" + dataset_name + "/" + pdb_name + "/" + pdb_name + ".pairs" 
        shutil.move(basepairs_file_path, "./Datasets/" + dataset_name + "/")
