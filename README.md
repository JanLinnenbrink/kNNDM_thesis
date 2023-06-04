# Development of a new spatial cross-validation strategy and its application for machine learning based modelling of ecological data

This repository contains the R-scripts used in the master thesis "Development of a new spatial cross-validation strategy and its application for machine learning based modelling of ecological data" by Jan Linnenbrink. The thesis is partly based on "kNNDM: k-fold Nearest Neighbour Distance Matching Cross-Validation for map accuracy estimation" by Jan Linnenbrink, Carles Milà, Marvin Ludwig and Hanna Meyer. The manuscript has been submitted to the journal *Geoscientific Model Development*.

The repository is structured as follows:

* [simulation_virtualSpecies](simulation_virtualSpecies/): This folder contains the code to reproduce the first simulation.
	* [figures.Rmd](simulation_virtualSpecies/figures.Rmd): R-Markdown file to reproduce the figures. [figures.pdf](figures.pdf) shows the output as one pdf document.
	* [code](simulation_virtualSpecies/code/): contains the R-code to reproduce all analyses.
	* [data](simulation_virtualSpecies/data/): contains the data used in the simulation.
	* [results](simulation_virtualSpecies/results/): contains the generated results.
	* [figures](simulation_virtualSpecies/figures/): contains the figures.

* [simulation_deBruin](simulation_deBruin/): This folder contains the code to reproduce the second simulation.
	* [figures_sim2.qmd](simulation_deBruin/igures_sim2.qmd): The R-Quarto-Script to generate the figures of simulation 2
	* [R](simulation_deBruin/R/): contains the R-scripts to run the deBruin simulations.
	* [material](simulation_deBruin/material/): contains the resulting tables and figures.

