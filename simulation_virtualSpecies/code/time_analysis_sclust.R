#-----------------------------------------------------------#
#             Simulation analysis: virtual species          #
#-----------------------------------------------------------#
.libPaths("/home/j/jlinnenb/r_packages/")
library("parallel")
library("doParallel")
library("pbapply")

setwd("~/kNNDM_paper/")

# Load utils, functions, and define number of iterations
source("code/time_functions.R")

# Read data
spoly <- st_read(dsn="data/species_vdata.gpkg", layer="sampling_polygon")
wclim <- rast("data/species_stack.grd")
wgrid <- st_read(dsn="data/species_vdata.gpkg", layer="landscape_grid")

# Launch simulation
set.seed(1234)
sims_rand <- sim_species(wgrid, wclim, spoly[1,], "sclust", ncores=30)
write_csv(sims, "results/time_res_sclust.csv")
