import sys
import os
import subprocess

rna_align_path="RMalign/RMalign"

dataset_path = sys.argv[1]
dataset_list=os.listdir(dataset_path)
comp_pairs_list=[]
for i in range(len(dataset_list)-1):
    for j in range(i+1,len(dataset_list)):
        comp_pairs_list.append((dataset_list[i],dataset_list[j]))

rmalign_to_dataset_path="tRNA_dataset/"

for pair in comp_pairs_list:
    first_elem=pair[0]
    second_elem=pair[1]
    first_elem_chain=first_elem[len(first_elem)-5]

    second_elem_chain=second_elem[len(second_elem)-5]
    result = subprocess.run([rna_align_path,"-A",rmalign_to_dataset_path + first_elem,"-Ac",
                             first_elem_chain,"-B",rmalign_to_dataset_path + second_elem,"-Bc",second_elem_chain ],capture_output=True, text=True)

    with open("all_to_all_outputs.txt", "a") as f:
        f.write(result.stdout)
with open("all_to_all_outputs.txt","a") as f:
    f.write("Total number of comparisons: " + str(len(comp_pairs_list)))