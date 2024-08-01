import os
import subprocess

dataset_pdbs = os.listdir("Dataset")

for pdb in dataset_pdbs:
    subprocess.run(["python3", "secondary_structure_generator.py", pdb])

