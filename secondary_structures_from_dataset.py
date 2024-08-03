# Script used for deriving 2D structures from PDBs in a given dataset.
# The script uses "secondary_structure_generator.py" for deriving the 2D structure for each PDB.

# USAGE: python3 secondary_structures_from_dataset.py <dataset_name> .

# Note: the dataset on wich the script runs must be located under the path "./Datasets".

import sys
import os
import subprocess

dataset_name = sys.argv[1]

dataset_pdbs = os.listdir("./Datasets/" + dataset_name)

for pdb in dataset_pdbs:
    subprocess.run(["python3", "secondary_structure_generator.py", pdb, dataset_name])

