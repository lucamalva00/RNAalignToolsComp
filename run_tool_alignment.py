import sys
import subprocess
import os
import csv
import time


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
        alignment_data = []

        first_elem = pair[0]
        second_elem = pair[1]
        first_elem_chain = first_elem[len(first_elem)-5]
        second_elem_chain = second_elem[len(second_elem)-5]

        result = subprocess.Popen([rna_align_path,"-A",rmalign_to_dataset_path + first_elem,"-Ac", first_elem_chain,"-B",rmalign_to_dataset_path + second_elem,"-Bc",second_elem_chain ],stdout=subprocess.PIPE)
        
        stdout_lines = result.stdout.read().decode().splitlines()

        alignment_data.append(first_elem[0]+first_elem[1]+first_elem[2]+first_elem[3])
        alignment_data.append(second_elem[0]+second_elem[1]+second_elem[2]+second_elem[3])
        
        for line in stdout_lines:
            if(line.startswith("Aligned length=")):
                splitted_line = line.split(",")
                splitted_splitted_line = splitted_line[1].split("=")
                alignment_data.append(splitted_splitted_line[1].strip())

        outputs.append(alignment_data)        

    with open("RMalign_output.csv", "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(outputs)

def exec_ARTS():
    arts_path = "./arts/arts"

    outputs = []

    for pair in pairs:
        alignment_data = []

        first_elem = pair[0]
        second_elem = pair[1]

        alignment_data.append(first_elem[0:4])
        alignment_data.append(second_elem[0:4])

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
                alignment_data.append(splitted_line[3].strip())
                break

        print(alignment_data)    

        
        outputs.append(alignment_data)
    
    with open("ARTS_output.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(outputs)

def exec_CLICK():
    data_base_path = "./data/tRNA_dataset/"

    os.system("chmod +rwx ./Click/Parameters.inp")
    

    outputs = []

    for pair in pairs:
        alignment_data = []

        first_elem = pair[0]
        second_elem = pair[1]

        alignment_data.append(first_elem[0:4])
        alignment_data.append(second_elem[0:4])

        process = subprocess.Popen(["./click", data_base_path + first_elem, data_base_path + second_elem, "-s 0"], stdout=subprocess.PIPE, cwd="./Click")
        process.wait()

        result_list = process.stdout.read().decode().splitlines()

        rmsd = result_list[1].split("=")[1]
        
        alignment_data.append(rmsd)
        outputs.append(alignment_data)

        dataset_dir = os.listdir("./Click/data/tRNA_dataset")

        for pdb in dataset_dir:
            if(pdb.endswith(".clique") or ("-" in pdb)):
                os.remove("./Click/data/tRNA_dataset/" + pdb)
    
    with open("CLICK_output.csv", "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(outputs)

def exec_SARA():
    output = []

    for pair in pairs:
        print(pairs.index(pair))
        error = False
        alignment_data = []

        first_chain = pair[0]
        second_chain = pair[1]

        first_chain_id = first_chain[len(first_chain) - 5]
        second_chain_id = second_chain[len(second_chain) - 5]

        alignment_data.append(first_chain[0:4])
        alignment_data.append(second_chain[0:4])

        sara_command = ["python", "runsara.py", "./tRNA_dataset/" + first_chain, first_chain_id, "./tRNA_dataset/" + second_chain, second_chain_id, "-s", "-o output.txt"]

        result = subprocess.Popen(sara_command, cwd="./SARA", stdout=subprocess.PIPE, text=True)
        result_list = result.stdout.readlines()
        
        for elem in result_list:
            if(elem.startswith("Error")):
                error = True

        if error:
            continue

        out_lines = []

        while not os.path.exists("./SARA/output.txt"):
            time.sleep(1)

        with open("./SARA/output.txt", "r") as f:
            out_lines = f.readlines()
        
        for line in out_lines:
            if line.startswith("REMARK SARARMS"):
                line_elems = line.split(" ")
                alignment_data.append(line_elems[len(line_elems) - 1].strip())

        output.append(alignment_data)

        os.remove("./SARA/output.txt")

        print(alignment_data)

    with open("SARA_output.csv", "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(output)
        


        


match tool_name:
    case "RMalign":
        exec_RMalign()
    case "ARTS":
        exec_ARTS()
    case "CLICK":
        exec_CLICK()
    case "SARA":
        exec_SARA()
