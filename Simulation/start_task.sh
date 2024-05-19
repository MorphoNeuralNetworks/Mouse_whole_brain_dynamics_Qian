#!/bin/bash
#SBATCH --job-name=start_task
#SBATCH --partition=debug
#SBATCH --output=./start_task.out

# init environment
cd /TVB_Distribution/tvb_data
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

# run python
cd /Simulation
/TVB_Distribution/tvb_data/bin/python /Simulation/main.py
