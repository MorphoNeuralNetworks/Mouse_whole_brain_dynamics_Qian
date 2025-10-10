import numpy as np

par_list=[]
for i in range(0,100):
    for j in range(0,200):
        par_list.append((i*0.001+0.05, j*0.000001+0.000001))

# save to file
with open("par_list.txt","w") as f:
    for (c,n) in par_list:
        f.write(f"{c} {n}\n")
