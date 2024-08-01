import os
import subprocess
import shutil

dataset_pdbs = os.listdir("Dataset")

for pdb in dataset_pdbs:
    if(not pdb.endswith(".pdb_new")):
        pdb_id = pdb[0:4]
        pdb_name = pdb[0:len(pdb) - 4]
        secondary_structure_path = "./PDBs_secondary_structures/" + pdb_name + "/" + pdb + ".out"
        subprocess.run(["python3", "secondary_structure_formatter.py", secondary_structure_path])
        basepairs_file_path = "./PDBs_secondary_structures/" + pdb_name + "/" + pdb_name + ".pairs" 
        shutil.move(basepairs_file_path, "./Dataset/")
