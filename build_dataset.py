# This scripts constructs the pdb dataset from a list of pdb IDs.
# Takes in input a text file with a comma separated list of pdb IDs and produces a folder containing the associated pdb files.
#
# USAGE: python3 build_dataset.py <pdb_IDs.txt> (where pdb_IDs.txt is a text file containing a comma separated list of pdb IDs).
# Note: This script uses "batch_download.sh" from RCSB to download pdb files given the associated IDs.

import sys
import subprocess
import os

# Function to execute bash script with the specified options
def run_bash_script(script_path, options):
    try:
        subprocess.run([script_path] + options, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running the bash script: {e}")
        return False
    return True

pdb_IDs_file = sys.argv[1]

os.mkdir("Dataset")

# Authorization command to execute batch_download.sh script 
auth_command = 'chmod +rwx batch_download.sh'
os.system(auth_command)

script_path = "./batch_download.sh"
options = ["-f", pdb_IDs_file, "-o", "Dataset", "-p"]

success = run_bash_script(script_path, options)

if success:
    dataset_tar_files = os.listdir("./Dataset")
    for tar_file in dataset_tar_files:
        subprocess.run(["gzip" , "-d" , "./Dataset/" + tar_file])
else:
    print("Error executing bash script")
