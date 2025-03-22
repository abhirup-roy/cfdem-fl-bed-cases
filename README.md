# CFDEM Fluidised Bed Cases

This repository allows users to test the validity of Bond Number of simulations with previous models - using the humble fluidised bed. This is part of [publication pending]. 
Contained in this repo the CFDEM code (written by of Hanqiao Che and adapted by myself) for carrying out a fluidised bed CFDEM simulation, as well as a Python pre/post-processing package for the simulation.
While the paper is pending publication, if you end up using this package - feel free to cite Hanqiao's paper on coarse-grained CFDEM simulations (this paper was crucial in the making of this project)

> 📘 **Evaluation of coarse-grained CFD-DEM models with the validation of PEPT measurements, Che et al. (2023)** <br>
> Hanqiao Che, Dominik Werner, Jonathan Seville, Tzany Kokalova Wheldon, Kit Windows-Yule, Evaluation of coarse-grained CFD-DEM models with the validation of PEPT measurements, Particuology,
Volume 82,
2023,
Pages 48-63,
ISSN 1674-2001,
https://doi.org/10.1016/j.partic.2022.12.018.

## Pre-Post Processing 
The prepost directory has a self-contained library  for pre-post processing. Documentation can be found [here](prepost/README.md)

Included in this repository are 3 Python example files. These can act as examples on how to import and use the library

* `run_templating.py` - This contains the format for intialising the LIGGGHTS scripts for the DEM case
* `plot_fluidn_curbes.py` - This shows how to plot fluidisation curves for the simulation
* `find_bondno.py` - This shows how to calculate the simulation and theoretical Bond number
