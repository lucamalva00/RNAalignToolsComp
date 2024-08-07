import sys
import os
import csv

dataset_name = sys.argv[1]
alignments_filename = sys.argv[2]

dataset_PDBs = os.listdir("./Datasets/" + dataset_name)

PDBs = []

for pdb_file in dataset_PDBs:
    if pdb_file.endswith(".pairs") or pdb_file.endswith(".pdb_new"):
        continue
    else:
        PDBs.append(pdb_file)

# list of the rows of csv alignments results. ([[pdb_name_a, pdb_name_b, distance], ...])
alignments = []
invalid_PDBs = []



with open(alignments_filename, "r", newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        alignments.append(row)


def validate_pdb(pdb_name):
    for row in alignments:
        if pdb_name in row:
            if len(row) == 3:
                return True
            continue
    return False

for pdb in PDBs:
    pdb_name = pdb[0:len(pdb) - 4]

    valid = validate_pdb(pdb_name) 

    if valid:
        continue
    else:
        invalid_PDBs.append(pdb_name)

print(invalid_PDBs)




