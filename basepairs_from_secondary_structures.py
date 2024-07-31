import os
import subprocess
import shutil

dataset_pdbs = os.listdir("dataset")

for pdb in dataset_pdbs:
    pdb_id = pdb[0:4]
    secondary_structure_path = "./PDBs_secondary_structures/" + pdb_id + "/" + pdb + ".out"
    subprocess.run(["python3", "secondary_structure_formatter.py", secondary_structure_path])
    basepairs_file_path = "./PDBs_secondary_structures/" + pdb_id + "/" + pdb[0:len(pdb) - 4] + ".pairs" 
    shutil.move(basepairs_file_path, "./dataset/")
