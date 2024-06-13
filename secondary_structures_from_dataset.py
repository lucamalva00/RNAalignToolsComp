import os
import subprocess

dataset_pdbs = os.listdir("tRNA_dataset")

for pdb in dataset_pdbs:
    pdb_id = pdb[0:4]
    subprocess.run(["python3", "secondary_structure_generator.py", pdb_id])

