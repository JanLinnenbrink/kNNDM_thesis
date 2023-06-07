# Development of a new spatial cross-validation strategy and its application for machine learning based modelling of ecological data

This repository contains the R-scripts used in the master thesis "Development of a new spatial cross-validation strategy and its application for machine learning based modelling of ecological data" by Jan Linnenbrink. The thesis is partly based on "kNNDM: k-fold Nearest Neighbour Distance Matching Cross-Validation for map accuracy estimation" by Jan Linnenbrink, Carles Milà, Marvin Ludwig and Hanna Meyer. The manuscript has been submitted to the journal *Geoscientific Model Development*.

The repository is structured as follows:

* [simulation_virtualSpecies](simulation_virtualSpecies/): This folder contains the code to reproduce the first simulation.
	* [figures.Rmd](simulation_virtualSpecies/figures.Rmd): R-Markdown file to reproduce the figures.
	* [code](simulation_virtualSpecies/code/): contains the R-code to reproduce all analyses.
	* [data](simulation_virtualSpecies/data/): contains the data used in the simulation.
	* [results](simulation_virtualSpecies/results/): contains the generated results.
	* [figures](simulation_virtualSpecies/figures/): contains the figures.

* [simulation_AGB_OCS](simulation_AGB_OCS/): This folder contains the code to reproduce the second simulation.
	* [summary.R](simulation_virtualSpecies/summary.R): The R-Script to generate the figures of simulation
	* [R](simulation_AGB_OCS/R/): contains the R-scripts to run the deBruin simulations.
	* [material](simulation_AGB_OCS/material/): contains the resulting tables and figures.

* [global_map](global_map/): This folder contains the code to reproduce the global map.
	* [code](global_map/code/): contains the R-code to reproduce the mapping study.
	* [results](global_map/reproduced_random/): contains the reproduced results based on random CV.
	* [results](global_map/reproduced_knndmcv/): contains the reproduced results based on kNNDM CV.
	* [figures](global_map/figures/): contains the figures.

