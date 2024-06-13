# This script use RNAVIEW to obtain the secondary structure from a pdb file
# Usage: 'python3 secondary_structure_generator.py <pdb_id>'

import sys
import os
import subprocess
import shutil

pdb_id = sys.argv[1]
pdb_list = os.listdir("tRNA_dataset")
path_file = "tRNA_dataset/"
rnaview_path = "RNAVIEW/bin/rnaview"
auth_command = 'chmod +rwx RNAVIEW/bin/rnaview'
auth_command2 = 'chmod +rx RNAVIEW/BASEPARS/misc_rna.par'
os.system(auth_command)
os.system(auth_command2)

for pdbname in pdb_list:
    if (pdbname.startswith(pdb_id) and  not pdbname.endswith("new")):
        path_file = path_file + pdbname
 

subprocess.run([rnaview_path, "-p", path_file])

os.makedirs("./PDBs_secondary_structures/" + pdb_id)

datasetlist = os.listdir("tRNA_dataset")

for filename in datasetlist:
    if (filename.endswith(".out") or filename.endswith(".xml") or filename.endswith(".ps")
        or filename.endswith("tmp.pdb")):
        shutil.move("tRNA_dataset/" + filename, "./PDBs_secondary_structures/" + pdb_id + "/" )

shutil.move("base_pair_statistics.out", "./PDBs_secondary_structures/" + pdb_id + "/")