import os

id_list=[]

with open("tRNA_idList.txt","r") as f:
    id_list=f.read().split(",")

download_id_list=os.listdir("Dataset")
downloaded=[]
for download_id in download_id_list:
    downloaded.append(download_id[0:4])
valid=[]
invalid=[]
for id in id_list:
    if(id in downloaded):
        valid.append(id)
    else:
        invalid.append(id)


with open("Dataset_info.txt","w") as f:
    for v in valid:
        f.write(v + " - SCARICATO\n")
    for i in invalid:
        f.write(i + " - NON SCARICATO\n")    
    f.write("Numero file scaricati: " + str(len(valid)) + "\n")
    f.write("Numero file non scaricati: " + str(len(invalid)) + "\n")
    f.write("ID totali: " + str(len(id_list)) + "\n")
