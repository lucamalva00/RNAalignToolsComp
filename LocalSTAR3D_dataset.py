import subprocess


pdb_list = []
with open("tRNA.txt","r") as file:
    pdb_list=file.readlines()

for pdb in pdb_list:
    splitted_list=pdb.split(" ")
    pdb_id=splitted_list[0]
    pdb_chain=splitted_list[1]
    subprocess.run(["java","-cp","./LocalSTAR3D/LocalSTAR3D.jar","Preprocess",pdb_id[0:4],pdb_chain])

