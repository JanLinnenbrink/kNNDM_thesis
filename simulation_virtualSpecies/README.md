# kNNDM: k-fold Nearest Neighbour Distance Matching Cross-Validation for map accuracy estimation

This repository contains the R-scripts used in the part of the master thesis that consists of the paper "kNNDM: k-fold Nearest Neighbour Distance Matching Cross-Validation for map accuracy estimation" by Jan Linnenbrink, Carles Milà, Marvin Ludwig and Hanna Meyer. The manuscript has been submitted to the journal *Geoscientific Model Development*.

The repository is structured as follows:

* [figures.Rmd](figures.Rmd): R-Markdown file to reproduce the figures. [figures.pdf](figures.pdf) shows the output as one pdf document.

* [code](code/): contains the R-code to reproduce all analyses.
	* [fig_utils.R](code/figures_utils.R): helper functions used in [figures.Rmd](figures.Rmd).
	* [sim_analysis.R](code/sim_analysis.R): R-Script to run the simulation.
	* [sim_functions.R](code/sim_functions.R): contains the functions called from [sim_analysis.R](code/sim_analysis.R).
	* [sim_utils.R](code/sim_utils.R): contains additional helper functions used in [sim_analysis.R](code/sim_analysis.R) and [sim_analysis_W.R](code/sim_analysis_W.R).
	* [sim_landscape.R](code/sim_landscape.R): R-Script to create the landscape dataset used in [sim_analysis.R](code/sim_analysis.R) and [sim_analysis_W.R](code/sim_analysis_W.R).
	* [sim_analysis_W.R](code/sim_analysis_W.R): R-Script to reproduce the simulation while keeping all outcomes of kNNDM and establishing an link between (CV estimated - True) ~ W.
	* [sim_functions_W.R](code/sim_functions_W.R): contains the functions called from [sim_analysis_W.R](code/sim_analysis_W.R).
	* [knndm_W.R](code/knndm_W.R): contains the modified kNNDM-function used in [sim_analysis_W.R](code/sim_analysis_W.R).
	* [time_functions.R](code/time_functions.R): contains the functions called from [time_analysis.R](code/time_analysis.R).
	* [time_analysis.R](code/time_analysis.R):  R-Script to run the computational speed comparison.

* [data](data/): contains the data used in the simulation.
	* [species_stack.grd](data/species_stack.grd): the predictor variables, stacked with the true outcome.
	* [species_vdata.gpkg](data/species_vdata.gpkg): the sampling area (layer="sampling_area") and the prediction grid (layer="landscape_grid").

* [results](results/): contains the generated results.
	* [sim_res.csv](results/sim_res.csv): the results from the simulation.
	* [sim_res_W.csv](results/sim_res_W.csv): the results from the modified simulation.
	* [time_res_rand.csv](results/time_res_rand.csv): the results from the speed comparison for randomly distributed training data.
	* [time_res_sclust.csv](results/time_res_sclust.csv): the results from the speed comparison for strongly clustered distributed training data.

* [figures](figures/): contains the figures.
	* [1_example_ecdfs.png](figures/1_example_ecdfs.png): Figure 1 from the manuscript. Shows different sample configurations and the corresponding NND ECDFs.
	* [2_method_knndm.pdf](figures/2_method_knndm.pdf): Figure 3 from the manuscript. Shows the general workflow of kNNDM CV.
	* [3_results_sim.pdf](figures/3_results_sim.pdf): Figure 5 from the manuscript. The results of the simulations.
	* [4_time_analysis.pdf](figures/4_time_analysis.pdf): Figure 6 from the manuscript. The computational speed comparison.
	* [5_results_W_err.pdf](figures/5_results_W_err.pdf): Figure 7 from the manuscript. The association between the absolute difference between the CV and true map accuracy statistics and W statistic.
