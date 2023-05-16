#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# set the number of CPU cores per node
#SBATCH --ntasks-per-node 20

# How much memory is needed (per node)
#SBATCH --mem=120GB

# set a partition
#SBATCH --partition long

# set max wallclock time
#SBATCH --time=60:00:00

# set name of job
#SBATCH --job-name=kndm

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# set an output file
#SBATCH --output output_knndm.dat

# send mail to this address
#SBATCH --mail-user=jlinnenb@uni-muenster.de

# run the application
module add palma/2021a
module add foss R GDAL
R CMD BATCH --vanilla CV_knndm.R
