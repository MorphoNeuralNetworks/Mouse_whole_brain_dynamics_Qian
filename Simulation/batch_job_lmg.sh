#!/bin/bash
#SBATCH --job-name=tvb_array
#SBATCH --array=0-19999
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=debug
#SBATCH --output=/Simulation/out/%A_%a.out

# --- init environment ---
cd /TVB_Distribution/tvb_data
export PATH=`pwd`/bin:$PATH
export PYTHONPATH=`pwd`/lib/python3.10:`pwd`/Lib/site-packages
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

# --- pick parameters for this task ---
read par1 par2 < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /Simulation/par_list.txt)

# --- run python simulation ---
cd /Simulation
/TVB_Distribution/tvb_data/bin/python /Simulation/RWWmodel.py $par1 $par2
