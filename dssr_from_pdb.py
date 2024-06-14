import subprocess
import os

pdb_list=os.listdir("./LocalSTAR3D/PDB")
for pdb in pdb_list:
    subprocess.run(["./LocalSTAR3D/STAR3D_struct_info/x3dna-dssr","-i=LocalSTAR3D/PDB/"+pdb,"-o=LocalSTAR3D/STAR3D_struct_info/"+pdb[0:4]+".dssr"])

files=os.listdir("./LocalSTAR3D/STAR3D_struct_info")
for file in files:
    if(not file.endswith(".dssr")):
        os.remove("./LocalSTAR3D/STAR3D_struct_info/"+file)
