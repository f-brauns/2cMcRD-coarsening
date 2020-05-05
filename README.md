# 2cMcRD-coarsening

This repository provides supplementary material for the publication http://arxiv.org/abs/2005.01495

 **Coarsening in (nearly) mass-conserving reaction-diffusion systems**

Authors: Fridtjof Brauns, Henrik Weyer, Jacob Halatek, Junghoon Yoon, Erwin Frey


## Mathematica Notebooks

Each notebook starts with a *Setup* section which must be evaulated first. Each setup section also contains a test (minimal usage example) for the provided functions.

*continuation-setup.nb* is a pure setup file that provides an functions and utils for pseudo-arclength continuation. It must remain in the same directory as *mesa-splitting_Brusselator.nb*, which loads it during setup. All other notebooks are independent.

## COMSOL simulations and Python evaluation scripts

The COMSOL files contain no data. Running the simulations and export the data requires COMSOL Version 5.4. 

Data is exported from COMSOL in a tabular format as txt file. For conversion to hdf5, Python scripts are prvided in ./evaluation-scripts/ (chage the file names in the scripts as needed before running). The evaulation scripts used to count the number of peaks and generate the plots laod the hdf5 files.

