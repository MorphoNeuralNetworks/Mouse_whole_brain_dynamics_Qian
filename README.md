# Mouse whole brain dynamics

This is the code for the CCN2024 paper "Single-cell morphological data provide refined simulations of resting-state" and Penghao Qian's Master thesis.  

## Main Idea
We simulated the whole mouse brain's resting state using an extensive dataset of 1876 fully reconstructed neurons, revealing stronger and more varied connections than previous tracer injection-based brain connectivity measurements indicated. After optimizing global coupling and background noise parameters, we tested the simulation's alignment with experimental data, finding that simulations using single-cell connectivity have increased predictive power compared to tracer-based connectomes. Our findings underscore the importance of incorporating detailed single-cell information to accurately model brain dynamics, offering insights into the mouse brain's functional architecture.   
  
Specifically we will include the following goals:  
1. Develop an algorithm to generate mesoscale connectivity based on single-cell morphology data and bouton locations. And compare to Allen mesoscale connectivity.

2. Analyze the simulation difference between the networks constructed using s-type or c-type.

3. Compare the differences in simulation under different perturbation on neuron morphology.

