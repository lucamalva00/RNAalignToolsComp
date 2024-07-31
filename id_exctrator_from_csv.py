import sys 
import csv

csv_name = sys.argv[1]

id_list=[]
with open(csv_name,"r") as f:
    file_reader=csv.reader(f,delimiter=";")
    for row in file_reader:
        id_list.append(row[0])

with open("tRNA_idList.txt","w") as file:
    for i in range(0,len(id_list)):
        if(i==len(id_list)-1):
            file.write(id_list[i])
            break
        file.write(id_list[i] + ",")