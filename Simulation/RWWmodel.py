from tvb.simulator.lab import *
from tvb.simulator.plot.tools import *
import numpy as np
import time,shutil
import os,sys

p = r'./ConnectivityAllen2mm.zip'
conn=connectivity.Connectivity.from_file(p)

coupling_v,noise_v= 0.05, 0.00005#float(sys.argv[1]),float(sys.argv[2])
sim = simulator.Simulator(
    model=models.ReducedWongWang(w=np.array([1.0]), I_o=np.array([0.3])),
    connectivity=conn,
    coupling=coupling.Linear(a=np.array([coupling_v])), # 0.096
    integrator=integrators.EulerStochastic(dt=0.1, noise=noise.Additive(nsig=np.array([noise_v]))), # 0.000013
    monitors=(monitors.Bold(period=1e3),monitors.TemporalAverage(period=1e3)),
    simulation_length=720e3
)
sim.configure()
# Run the simulation
(bold_time, bold_data), _ = sim.run()
os.makedirs("BOLD", exist_ok=True)
writepath=r"./BOLD/"
file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
np.save(writepath+"bold_time"+file_end+".npy",bold_time)
np.save(writepath+"bold_data"+file_end+".npy",bold_data)