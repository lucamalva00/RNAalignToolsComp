#Thi script use RNAVIEW to obtain the secondary structure from a pdb file
#Usage: 'python3 2D_structure_generator.py <pdb_id>'

import sys
import os
import subprocess
import shutil
pdb_id = sys.argv[1]
pdb_list = os.listdir("tRNA_dataset")
path_file = "tRNA_dataset/"
rnaview_path = "RNAView-master/bin/rnaview"
auth_command = 'chmod +rwx RNAView-master/bin/rnaview'
os.system(auth_command)

for pdbname in pdb_list:
    if (pdbname.startswith(pdb_id)):
        path_file=path_file + pdbname

rnaview_file_path = "../" + path_file 

subprocess.run([rnaview_path, "-p", path_file])

os.mkdir("extracted_chains_2d_structure/"+ pdb_id + "_2d_structure")
datasetlist = os.listdir("tRNA_dataset")

for filename in datasetlist:
    if (filename.endswith(".out") or filename.endswith(".xml") or filename.endswith(".ps")
        or filename.endswith("tmp.pdb")):
        shutil.move("tRNA_dataset/"+filename,"extracted_chains_2d_structure/"+ pdb_id + "_2d_structure" )

shutil.move("base_pair_statistics.out","extracted_chains_2d_structure/"+ pdb_id + "_2d_structure")