import sys
import subprocess
import os
import signal
import csv
import time
import shutil


# Constructs the list consisting of all pdb pairs from the given dataset and returns it
def construct_pairs(dataset_path):
    pairs = []
    files_list = os.listdir(dataset_path)
    pdb_list = []
    

    for pdb_elem in files_list:
        if(not pdb_elem.endswith(".pairs") and not pdb_elem.endswith(".pdb_new")):
            pdb_list.append(pdb_elem)
    
    for i in range(len(pdb_list)-1):
        for j in range(i+1,len(pdb_list)):
            pairs.append((pdb_list[i],pdb_list[j]))
    return pairs


# Prints the current alignment performed
def print_alignment(index, first_elem, second_elem):
    print(str(index) + ": Performing alignment between " + first_elem + " and " + second_elem + "\n")

tool_name = sys.argv[1]
dataset_name = sys.argv[2]

dataset_path = "./Datasets/" + dataset_name

pairs = construct_pairs(dataset_path)

def exec_RMalign():
    rna_align_path = "RMalign/RMalign"

    outputs = []

    for pair in pairs:
        alignment_data = []

        first_elem = pair[0]
        second_elem = pair[1]

        first_elem_chain = first_elem[5]
        second_elem_chain = second_elem[5]

        print_alignment(pairs.index(pair), first_elem, second_elem)

        result = subprocess.Popen([rna_align_path, "-A", dataset_path + "/" + first_elem, "-Ac", first_elem_chain, "-B", dataset_path + "/" + second_elem, "-Bc", second_elem_chain], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        stdout_lines = result.stdout.read().decode().splitlines()

        alignment_data.append(first_elem[0:len(first_elem) - 4])
        alignment_data.append(second_elem[0:len(second_elem) - 4])
        
        for line in stdout_lines:
            if(line.startswith("Aligned length=")):
                splitted_line = line.split(",")
                splitted_splitted_line = splitted_line[1].split("=")
                alignment_data.append(splitted_splitted_line[1].strip())

        outputs.append(alignment_data)        

    with open(dataset_name + "_RMalign_output.csv", "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(outputs)

def exec_ARTS():
    arts_path = "./arts/arts"

    outputs = []

    for pair in pairs:
        alignment_data = []

        first_elem = pair[0]
        second_elem = pair[1]

        alignment_data.append(first_elem[0:len(first_elem) - 4])
        alignment_data.append(second_elem[0:len(second_elem) - 4])

        print_alignment(pairs.index(pair), first_elem, second_elem)

        result = subprocess.run([arts_path, dataset_path + "/" + first_elem, dataset_path + "/" + second_elem], capture_output=True)

        print("Return code: " + str(result.returncode) + "\n")

        if(result.returncode != 0):
            continue

        arts_out = []

        while not os.path.exists("arts.out"):
            time.sleep(1)
        
        with open("arts.out", "r") as f:
            arts_out = f.readlines()
            print("arts.out opened\n")

        os.remove("arts.out")

        for line in arts_out:
            if (line.startswith("*** RMSD")):
                print(line + "\n")
                splitted_line = line.split(" ")
                alignment_data.append(splitted_line[3].strip())
                break

        print(alignment_data)    
        outputs.append(alignment_data)
    
    with open(dataset_name + "_ARTS_alignments.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(outputs)

def exec_CLICK():
    command_dataset_path = "./data/" + dataset_name
    data_base_path = "./Click/data/" + dataset_name 
    if(not os.path.exists(data_base_path)):
        shutil.copytree(dataset_path,data_base_path + "/")
    os.system("chmod +rwx ./Click/Parameters.inp")
    

    outputs = []

    for pair in pairs:
        alignment_data = []

        first_elem = pair[0]
        second_elem = pair[1]

        alignment_data.append(first_elem[0:4])
        alignment_data.append(second_elem[0:4])

        print_alignment(pairs.index(pair), first_elem, second_elem)

        process = subprocess.Popen(["./click", command_dataset_path + "/" + first_elem, command_dataset_path + "/" + second_elem, "-s", "0"], stdout=subprocess.PIPE, cwd="./Click")
        process.wait()

        result_list = process.stdout.read().decode().splitlines()

        rmsd = result_list[1].split("=")[1]
        
        alignment_data.append(rmsd)
        outputs.append(alignment_data)
        dataset_dir = os.listdir("./Click/data/" + dataset_name)

        for pdb in dataset_dir:
            if(pdb.endswith(".clique") or ("-" in pdb)):
                os.remove("./Click/data/" + dataset_name + "/" + pdb)
    
    with open(dataset_name + "_CLICK_alignments.csv", "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(outputs)

def exec_SARA():
    # invalid pdb list
    invalid_pdbs = ["5mre_b.pdb", "5lzs_i.pdb", "5t2c_A.pdb"]
    output = []
    errors_log = []
    command_dataset_path = "./" + dataset_name + "/"
    sara_dataset_path = "./SARA/" + dataset_name

    if(not os.path.exists(sara_dataset_path)):
        shutil.copytree(dataset_path,sara_dataset_path + "/")


    for pair in pairs:
        error = False
        alignment_data = []

        first_elem = pair[0]
        second_elem = pair[1]

        print_alignment(pairs.index(pair), first_elem[0:len(first_elem) - 4], second_elem[0:len(second_elem) - 4])

        if (first_elem in invalid_pdbs) or (second_elem in invalid_pdbs):
            output.append(alignment_data)
            continue

        first_chain_id = first_elem[5]
        second_chain_id = second_elem[5]

        

        alignment_data.append(first_elem[0:len(first_elem) - 4])
        alignment_data.append(second_elem[0:len(second_elem) - 4])

        sara_command = ["python", "runsara.py", command_dataset_path + first_elem, first_chain_id, command_dataset_path + second_elem, second_chain_id, "-a", "C3\'", "-s", "-o", "output.txt"]
        
        process = subprocess.Popen(sara_command, cwd="./SARA", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()

        if err:
            errors_log.append((pairs.index(pair), err.decode()))
            alignment_data.append('')
            output.append(alignment_data)
            continue

        result_list = out.decode().splitlines()
        
        for elem in result_list:
            if(elem.startswith("Error")):
                errors_log.append((pairs.index(pair), elem))
                alignment_data.append('')
                output.append(alignment_data)
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

    with open(dataset_name + "_SARA_alignments.csv", "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(output)
    
    with open("SARA_errors_log.txt", "w") as errors_f:
        for log in errors_log:
            errors_f.write("On alignment n." + str(log[0]) + ": " + log[1] + "\n")
        

def exec_USalign():
    output = []

    for pair in pairs:
        alignment_data = []

        first_elem = pair[0]
        second_elem = pair[1]

        alignment_data.append(first_elem[0:len(first_elem) - 4])
        alignment_data.append(second_elem[0:len(second_elem) - 4])

        print_alignment(pairs.index(pair), first_elem[0:len(first_elem) - 4], second_elem[0:len(second_elem) - 4])


        process = subprocess.Popen(["./USalign/USalign", dataset_path + "/" + first_elem, dataset_path + "/" + second_elem, "-mol", "RNA"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()

        print("Return code: " + str(process.returncode) + "\n")

        out_lines = process.stdout.read().decode().splitlines()

        for line in out_lines:
            if line.startswith("Aligned length"):
                rmsd_string = line.split(",")[1]
                rmsd = rmsd_string.split("=")[1].strip()
                alignment_data.append(rmsd)

        output.append(alignment_data)
    
    with open(dataset_name + "_USalign_alignments.csv", "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(output)

def exec_STAlign():
    output = []

    for pair in pairs:
        first_elem = pair[0]
        second_elem = pair[1]

        alignment_data = [first_elem, second_elem]

        threshold = 5 

        print_alignment(pairs.index(pair), first_elem, second_elem)

        cmd = ["java", "-jar", "./STalign/STAlign.jar", "-d", "-af", dataset_path + "/" + first_elem, dataset_path + "/" + second_elem, "-t", str(threshold)]

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        out, err = process.communicate()

        output_lines = out.decode().splitlines()

        while len(output_lines) < 2:
            threshold += 1
            if threshold > 15:
                break
            print("Using threshold " + str(threshold) + "\n")
            cmd[len(cmd) - 1] = str(threshold)
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            outs, errs = process.communicate()
            output_lines = outs.decode().splitlines()

        if len(output_lines) == 1:
            output.append(alignment_data)
            continue

        distance = output_lines[1].split("=")[1].strip()
        
        alignment_data.append(distance)
        alignment_data.append(str(threshold))

        output.append(alignment_data)

    with open(dataset_name + "_STAlign_alignments.csv", "w", newline='') as file:
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
    case "USalign":
        exec_USalign()
    case "STAlign":
        exec_STAlign()
