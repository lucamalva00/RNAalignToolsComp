import sys
import subprocess
import os
import csv
#import secondary_structure_formatter
#import secondary_structure_generator

tool_name = sys.argv[1]
dataset_path=sys.argv[2]

def exec_RMalign():
    rna_align_path="RMalign/RMalign"
    dataset_list=os.listdir(dataset_path)
    comp_pairs_list=[]
    for i in range(len(dataset_list)-1):
        for j in range(i+1,len(dataset_list)):
            comp_pairs_list.append((dataset_list[i],dataset_list[j]))

    rmalign_to_dataset_path="tRNA_dataset/"
    outputs_list=[]

    for pair in comp_pairs_list:
        alignment_list=[]
        first_elem=pair[0]
        second_elem=pair[1]
        first_elem_chain=first_elem[len(first_elem)-5]

        second_elem_chain=second_elem[len(second_elem)-5]
        result = subprocess.Popen([rna_align_path,"-A",rmalign_to_dataset_path + first_elem,"-Ac",
                                first_elem_chain,"-B",rmalign_to_dataset_path + second_elem,"-Bc",second_elem_chain ],stdout=subprocess.PIPE)
        stdout_list=result.stdout.read().decode().splitlines()
        alignment_list.append(first_elem[0]+first_elem[1]+first_elem[2]+first_elem[3])
        alignment_list.append(second_elem[0]+second_elem[1]+second_elem[2]+second_elem[3])
        
        for line in stdout_list:
            if(line.startswith("Aligned length=")):
                splitted_line=line.split(",")
                splitted_splitted_line=splitted_line[1].split("=")
                alignment_list.append(splitted_splitted_line[1].strip())
        outputs_list.append(alignment_list)        

    with open("RMalign_output.csv","w",newline='') as file:
        writer=csv.writer(file)
        writer.writerows(outputs_list)

def exec_ARTS():
    arts_path="arts/arts"
    dataset_list=os.listdir(dataset_path)
    comp_pairs_list=[]
    for i in range(len(dataset_list)-1):
        for j in range(i+1,len(dataset_list)):
            comp_pairs_list.append((dataset_list[i],dataset_list[j]))

    arts_to_dataset_path="tRNA_dataset/"
    outputs_list=[]
    for pdb in dataset_list:
        pdb_id=pdb[0]+pdb[1]+pdb[2]+pdb[3]
        subprocess.run(["python3","secondary_structure_generator.py",pdb_id])


match tool_name:
    case "RMalign":
        exec_RMalign()
    case "ARTS":
        exec_ARTS()
