# This script use RNAVIEW to obtain the secondary structure from a pdb file

# Usage: 'python3 secondary_structure_generator.py <pdb_name> <dataset_name>'.
# Note: the dataset from wich generate the 2D structures must be located in the "./Dataset" directory.

import sys
import os
import subprocess
import shutil

pdb_name = sys.argv[1]
dataset_name = sys.argv[2]

dataset_path = "./Datasets/" + dataset_name
pdb_list = os.listdir(dataset_path)
path_file = ""
rnaview_path = "RNAVIEW/bin/rnaview"
auth_command = 'chmod +rwx RNAVIEW/bin/rnaview'
auth_command2 = 'chmod +rx RNAVIEW/BASEPARS/misc_rna.par'
os.system(auth_command)
os.system(auth_command2)

for pdbname in pdb_list:
    if (pdbname.startswith(pdb_name) and  not pdbname.endswith("new")):
        path_file += dataset_path + "/" + pdbname
 

subprocess.run([rnaview_path, "-p", path_file])

os.makedirs("./2D_structures/" + dataset_name + "/" + pdb_name[0:len(pdb_name) - 4])

datasetlist = os.listdir(dataset_path)

for filename in datasetlist:
    if (filename.endswith(".out") or filename.endswith(".xml") or filename.endswith(".ps")
        or filename.endswith("tmp.pdb")):
        shutil.move(dataset_path + "/" + filename, "./2D_structures/" + dataset_name + "/" + pdb_name[0:len(pdb_name)-4] + "/" )

shutil.move("base_pair_statistics.out", "./2D_structures/" + dataset_name + "/" + pdb_name[0:len(pdb_name)-4] + "/")