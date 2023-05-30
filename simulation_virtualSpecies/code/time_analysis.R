#-----------------------------------------------------------#
#             Simulation analysis: virtual species          #
#-----------------------------------------------------------#

# Load utils, functions, and define number of iterations
source("code/time_functions.R")

# Read data
spoly <- st_read(dsn="data/species_vdata.gpkg", layer="sampling_polygon")
wclim <- rast("data/species_stack.grd")
wgrid <- st_read(dsn="data/species_vdata.gpkg", layer="landscape_grid")

# Launch simulation for strongly clustered design
set.seed(1234)
sims <- sim_species(wgrid, wclim, spoly[1,], "sclust",interval=1:50)
write_csv(sims, "results/time_res_sclust.csv")

# Launch simulations for random design
set.seed(1234)
sims <- sim_species(wgrid, wclim, spoly[1,], "random",interval=1:50)
write_csv(sims, "results/time_res_random.csv")


