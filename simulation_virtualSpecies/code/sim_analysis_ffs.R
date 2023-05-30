#-----------------------------------------------------------#
#             Simulation analysis 2: virtual species        #
#-----------------------------------------------------------#
library("mcreplicate")

# Load utils, functions, and define number of iterations
source("code/sim_functions_ffs.R")
nsim <- 2
pboptions(type = "timer")

# Read data
spoly <- st_read(dsn="data/species_vdata.gpkg", layer="sampling_polygon")
wclim <- rast("data/species_stack_coords.tif")
wgrid <- st_read(dsn="data/species_vdata_coords.gpkg", layer="landscape_grid")

.wclim <- wrap(wclim)

# Launch simulation
set.seed(1234)
spoly=spoly[1,]

library(future.apply)
plan(multisession, workers=10)
sims <- future_replicate(n=nsim, sim_species(wgrid, rast(.wclim), spoly), simplify=FALSE)
sims_t <- apply(sims, 2, as.data.frame) 

write_csv(do.call(rbind, sims_t), "results/ffs.csv")
