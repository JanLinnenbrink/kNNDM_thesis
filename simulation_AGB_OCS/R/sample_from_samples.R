setwd("scratch/tmp/jlinnenb/deBruin_add_nndm/")
save_to <- "scratch/tmp/jlinnenb/deBruin_add_nndm/samples"
samples_root <- "scratch/tmp/jlinnenb/samples/"
samples <- 5000 # adjust to 5000

for (method in c("clusterGapped", "clusterMedium", "simpleRandom", "clusterStrong", "regular")) {
  real_AGB <- list.files(file.path(samples_root, method), glob2rx("AGBdata*.Rdata"))
  real_OCS <- list.files(file.path(samples_root, method), glob2rx("OCSdata*.Rdata"))
  real_pts <- list.files(file.path(samples_root, method), glob2rx("*_coords.Rdata"))
  for (i in 1:length(real_AGB)) {
    load(file.path(samples_root, method, real_AGB[i]))
    load(file.path(samples_root, method, real_OCS[i]))
    load(file.path(samples_root, method, real_pts[i]))
    
    row_numbers <- sample(1:nrow(AGBdata), samples)
    
    AGBdata <- AGBdata[row_numbers,]
    OCSdata <- OCSdata[row_numbers,]
    pts     <- pts[row_numbers,]
    if("ID" %in% names(AGBdata)) {AGBdata <- AGBdata[,! names(AGBdata) %in% c("ID")]}
    if("ID" %in% names(OCSdata)) {OCSdata <- OCSdata[,! names(OCSdata) %in% c("ID")]}
    
    save(AGBdata, file = file.path(save_to, method, real_AGB[i]))
    save(OCSdata, file = file.path(save_to, method, real_OCS[i]))
    save(pts, file = file.path(save_to, method, real_pts[i]))
  }
}
