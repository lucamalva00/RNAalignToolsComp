import sys
import subprocess
import os
import csv
import time
#import secondary_structure_formatter
#import secondary_structure_generator

# Constructs the list consisting of all pdb pairs from the given dataset and returns it
def construct_pairs(dataset_path):
    pairs = []
    files_list = os.listdir(dataset_path)
    pdb_list = []
    

    for pdb_elem in files_list:
        if(not pdb_elem.endswith(".pairs")):
            pdb_list.append(pdb_elem)
            
    for i in range(len(pdb_list)-1):
        for j in range(i+1,len(pdb_list)):
            pairs.append((pdb_list[i],pdb_list[j]))
    return pairs

tool_name = sys.argv[1]
dataset_path = sys.argv[2]

pairs = construct_pairs(dataset_path)

def exec_RMalign():
    rna_align_path = "RMalign/RMalign"

    rmalign_to_dataset_path="tRNA_dataset/"

    outputs = []

    for pair in pairs:
        alignment_list = []

        first_elem = pair[0]
        second_elem = pair[1]
        first_elem_chain = first_elem[len(first_elem)-5]
        second_elem_chain = second_elem[len(second_elem)-5]

        result = subprocess.Popen([rna_align_path,"-A",rmalign_to_dataset_path + first_elem,"-Ac", first_elem_chain,"-B",rmalign_to_dataset_path + second_elem,"-Bc",second_elem_chain ],stdout=subprocess.PIPE)
        
        stdout_lines = result.stdout.read().decode().splitlines()

        alignment_list.append(first_elem[0]+first_elem[1]+first_elem[2]+first_elem[3])
        alignment_list.append(second_elem[0]+second_elem[1]+second_elem[2]+second_elem[3])
        
        for line in stdout_lines:
            if(line.startswith("Aligned length=")):
                splitted_line = line.split(",")
                splitted_splitted_line = splitted_line[1].split("=")
                alignment_list.append(splitted_splitted_line[1].strip())

        outputs.append(alignment_list)        

    with open("RMalign_output.csv", "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(outputs)

def exec_ARTS():
    arts_path = "./arts/arts"

    outputs = []
    print(len(pairs))
    for pair in pairs:
        print(pairs.index(pair))
        alignment_list = []

        first_elem = pair[0]
        second_elem = pair[1]

        alignment_list.append(first_elem[0:4])
        alignment_list.append(second_elem[0:4])

        result = subprocess.run([arts_path, "./tRNA_dataset/" + first_elem, "./tRNA_dataset/" + second_elem], capture_output=True)
        if(result.returncode==255):
            continue

        arts_out = []

        while not os.path.exists("arts.out"):
            time.sleep(1)
        
        with open("arts.out", "r") as f:
            arts_out = f.readlines()
            print("opened")

        os.remove("arts.out")

        for line in arts_out:
            if (line.startswith("*** RMSD")):
                print(line)
                splitted_line = line.split(" ")
                alignment_list.append(splitted_line[3].strip())
                break

        print(alignment_list)    

        
        outputs.append(alignment_list)
    
    with open("ARTS_output.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(outputs)



    


match tool_name:
    case "RMalign":
        exec_RMalign()
    case "ARTS":
        exec_ARTS()
