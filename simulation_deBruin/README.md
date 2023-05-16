# Add NNDM CV to deBruin et. al

Reproduce results of the paper "Dealing with clustered samples for assessing map accuracy by cross-validation", deBruin et. al, 2022, and add a Nearest Neighbour Distance Matching CV to compare the results.

## Run

The overall computational load is reduced here.

1. The original sampling scripts are run (100 samples for 5 different sampling designs with 5000 samples each)
2. 1000 samples are taken per sample (run `R/sample_from_samples.R`)
3. CV methods are run and saved in `CVresults`. Most methods use 5 instead of 100 CV repetitions per sample.
   1. exhaustive: `CV_exhaustive_single.R`. The script is started individually 100 times.
   2. random: `CV_random_jon.R`. The script is started once. 5*CV, 100 samples. Yields 2000 files (incl `pts*` files for the model based technique)
   3. spatial: `CV_spatial_jon.R`. The script is started once. 5*CV, 100 samples. Yields 1000 files
   4. model_based:
   5. 