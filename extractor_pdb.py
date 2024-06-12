import sys
import subprocess
import os
from Bio.PDB import PDBParser, PDBIO


def remove(string):
    return string.replace(" ","")

def run_bash_script(script_path, options):
    # Function to execute bash script with the specified options
    try:
        subprocess.run([script_path] + options, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running the bash script: {e}")
        return False
    return True

# Obtain pdb id from command line args
pdb_id = sys.argv[1]
# Obtain chain id from command line args
chain_id = sys.argv[2]


# Write the id in a new text file 
with open("pdb_id.txt", "w") as file_with_id:
    file_with_id.write(pdb_id)

# Execute bash script to download pdb file from the id
script_path = "./batch_download.sh"
options = ["-f", "pdb_id.txt", "-p"]
tar_file = pdb_id + ".pdb.gz"
success = run_bash_script(script_path, options)
if success:
    # Unzip downloaded pdb file
    subprocess.run(["gzip" , "-d" , tar_file])
    # Remove temp text file containing the pdb id
    os.remove("pdb_id.txt")
    print("Bash script executed successfully.")
else:
    print("Error executing bash script.")

io = PDBIO()
parser = PDBParser()

# Parse the downloaded pdb file 
source_pdb = parser.get_structure(pdb_id, pdb_id + '.pdb')

# Select correct chain and save the sub-structure extracted in a pdb file
for chain in source_pdb.get_chains():
    if chain.get_id() == chain_id:
        io.set_structure(chain)
        io.save(source_pdb.get_id() + "_chain_" + chain.get_id() + ".pdb")

os.remove(pdb_id + '.pdb')