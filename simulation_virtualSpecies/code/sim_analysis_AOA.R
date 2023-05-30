#-----------------------------------------------------------#
#             Simulation analysis 2: virtual species        #
#-----------------------------------------------------------#
#.libPaths("/home/j/jlinnenb/r_packages/")
library("parallel")
library("doParallel")
library("pbapply")
setwd("C:/0_Msc_Loek/M7_Fernerkundung/NNDMpaper-main/")

# Load utils, functions, and define number of iterations
source("./code/sim2_functions_global.R")
nsim <- 100
pboptions(type = "timer")


# Read data
spoly <- st_read(dsn="data/species_vdata.gpkg", layer="sampling_polygon")
wclim <- rast("data/species_stack.grd")
wgrid <- st_read(dsn="data/species_vdata.gpkg", layer="landscape_grid")

# Prepare parallelization
print(paste0("Process started with ", 11, " cores."))
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

# Launch simulation
set.seed(1234)
sims <- pbreplicate(nsim, sim_species(wgrid, wclim, spoly), simplify=FALSE)

# We're done
stopCluster(cl)
rm("cl")
write_csv(do.call(rbind, sims), "results/sim2_species/global.csv")
