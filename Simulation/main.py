import subprocess
import time
import shutil
import os
def func(par):
    coupling_v,noise_v=par[0],par[1]
    file_end="w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    writepath="/Simulation/bash/"
    bash_file = open(writepath+f"{file_end}.sh",'w')
    bash_file.write(f"#!/bin/bash\n")
    bash_file.write(f"#SBATCH --job-name={file_end}\n")
    bash_file.write("#SBATCH --nodes=1\n")
    bash_file.write("#SBATCH --ntasks-per-node=1\n")
    bash_file.write("#SBATCH --cpus-per-task=1\n")
    # bash_file.write("#SBATCH --mem-per-cpu=2G\n")
    bash_file.write("#SBATCH --partition=debug\n")
    bash_file.write(f"#SBATCH --output=/Simulation/out/{file_end}.out\n")
    pp="""
# init environment
cd ./tvb_data
export PATH=`pwd`/bin:$PATH
export PYTHONPATH=`pwd`/lib/python3.10:`pwd`/lib/python3.10/site-packages

if [ ${LD_LIBRARY_PATH+1} ]; then
  export LD_LIBRARY_PATH=`pwd`/lib:`pwd`/bin:$LD_LIBRARY_PATH
else
  export LD_LIBRARY_PATH=`pwd`/lib:`pwd`/bin
fi
if [ ${LD_RUN_PATH+1} ]; then
  export LD_RUN_PATH=`pwd`/lib:`pwd`/bin:$LD_RUN_PATH
else
  export LD_RUN_PATH=`pwd`/lib:`pwd`/bin
fi
cd ../bin
"""
    bash_file.write(pp+"\n")
    
    pp=f"""
par1={round(coupling_v,4)}
par2={round(noise_v,6)}

# run python
cd /Simulation
/TVB_Distribution/tvb_data/bin/python /Simulation/RWWmodel.py $par1 $par2
"""
    bash_file.write(pp+"\n")
    bash_file.close()

writepath=r"/Simulation/bash/"
if os.path.exists(writepath):
    shutil.rmtree(writepath)  
os.mkdir(writepath) 

writepath=r"/Simulation/out/"
if os.path.exists(writepath):
    shutil.rmtree(writepath)  
os.mkdir(writepath) 

writepath=r"/Simulation/BOLD/"
if os.path.exists(writepath):
    shutil.rmtree(writepath)  
os.mkdir(writepath) 

# 0.096 [0-0.5]
# TVMB: 0.0051 x^2/2 0.000013 
# PANS: 0.015        0.00011 [0-]

par_list=[]
for i in range(0,100):
    for j in range(0,200):
        par_list.append((i*0.001+0.05,j*0.000001+0.000001))

for par in par_list:
    func(par)
    coupling_v,noise_v=par[0],par[1]
    file_end="w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    subprocess.run(f"sbatch /Simulation/bash/{file_end}.sh", shell=True)  

     