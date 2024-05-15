import sys
import subprocess
import os
import gzip
import string

def remove(string):
    # Funzione per rimuovere gli spazi da una stringa
    return string.replace(" ","")

def run_bash_script(script_path, options):
    # Funzione per eseguire uno script bash con le opzioni specificate
    try:
        subprocess.run([script_path] + options, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running the bash script: {e}")
        return False
    return True

# Ottieni gli ID PDB e della catena da riga di comando
pdb_id = sys.argv[1]
chain_id = sys.argv[2]


# Salva l'ID PDB in un file
with open("pdb_id.txt", "w") as file_with_id:
    file_with_id.write(pdb_id)

# Esegui lo script bash per scaricare il file PDB
script_path = "./batch_download.sh"
options = ["-f", "pdb_id.txt", "-p"]
tar_file = pdb_id + ".pdb.gz"
success = run_bash_script(script_path, options)
if success:
    # Decomprimi il file PDB scaricato
    subprocess.run(["gzip" , "-d" , tar_file])
    # Rimuovi il file temporaneo con l'ID PDB
    os.remove("pdb_id.txt")
    print("Bash script executed successfully.")
else:
    print("Error executing bash script.")

# Apri il file PDB completo
complete_pdb_file = open(pdb_id + ".pdb", "r")
# Inserisci ogni riga del file come elemento della lista lines
lines = complete_pdb_file.readlines()
complete_pdb_file.close()
os.remove(pdb_id + ".pdb")

extracted_lines = []

# Estrai le righe relative alla catena specificata
for line in lines:
    splitted_line = line.split()
    #In case that atom type is separated by a space
    if len(splitted_line) > 12:
        composed_elem = ""
        if len(splitted_line[3]) == 1:
            composed_elem = splitted_line[3] + "  " + splitted_line[4]
        else:
            composed_elem = splitted_line[3] + " " + splitted_line[4]
        splitted_line[3] = composed_elem
        splitted_line.remove(splitted_line[4])
    
    if splitted_line[0] == "ATOM": # Verifica che sia la parte interessata del file
        if splitted_line[4] == chain_id:  # Verifica se la riga corrisponde alla catena specificata
            extracted_line = (
                '{}   {:>4}  {:<4}  {} {}  {:>2}     {:>7}  {:<7}  {:>7}  {:>3} {:>4}           {}'.format(
                    splitted_line[0], splitted_line[1], splitted_line[2], splitted_line[3],
                    splitted_line[4], splitted_line[5], splitted_line[6], splitted_line[7],
                    splitted_line[8], splitted_line[9], splitted_line[10], splitted_line[11]
                )
            )  # Formatta la riga estratta
            extracted_lines.append(extracted_line)  # Aggiungi la riga estratta all'elenco

# Scrivi le righe estratte in un nuovo file PDB
with open("extracted_chain_" + pdb_id + ".pdb", "w") as extracted_chain:
    for line in extracted_lines:
        extracted_chain.write(line + "\n")
