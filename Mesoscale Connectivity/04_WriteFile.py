import numpy as np
import os,csv
with open('../Data/SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Re=[str(x[0]) for x in temp]

# connection=np.load("./Allen_connectivity.npy",allow_pickle=True)
connection=np.load("./SingleCell_connectivity.npy",allow_pickle=True)

connection=connection/np.amax(connection)

centers=np.load("region_center.npy",allow_pickle=True)
tract_lengths=np.load("tract_lengths.npy",allow_pickle=True)

## weights.txt
new_contents=[]
for x in connection:
    tt=""
    for j in x:
        tt=tt+str(j)+" "
    tt=tt[0:-1]+'\n'
    new_contents.append(tt)
write_path=r".\allen_2mm\ConnectivityAllen2mm_single"
name="weights.txt"
f=open(os.path.join(write_path, name),'w+')
f.writelines(new_contents)
f.close()

## centres.txt
new_contents=[]
name=[x+"_R" for x in TVMB_Re ]+[x+"_L" for x in TVMB_Re]
count=0
for x in centers:
    tt=name[count]+" "
    count+=1
    for j in x:
        tt=tt+str(j)+" "
    tt=tt[0:-1]+'\n'
    new_contents.append(tt)
write_path=r".\allen_2mm\ConnectivityAllen2mm_single"
name="centres.txt"
f=open(os.path.join(write_path, name),'w+')
f.writelines(new_contents)
f.close()

## tract_lengths.txt
new_contents=[]
for x in tract_lengths:
    tt=""
    for j in x:
        tt=tt+str(j)+" "
    tt=tt[0:-1]+'\n'
    new_contents.append(tt)
write_path=r".\allen_2mm\ConnectivityAllen2mm_single"
name="tract_lengths.txt"
f=open(os.path.join(write_path, name),'w+')
f.writelines(new_contents)
f.close()
