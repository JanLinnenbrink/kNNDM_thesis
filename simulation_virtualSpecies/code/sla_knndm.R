library(caret)
library(CAST)
library(sf)
library(terra)
library(dplyr)
library(rnaturalearth)
library("tidyterra")
library("cowplot")
source("code/figures_utils.R")
source("./code/repr_functions.R")

ee <- st_crs("+proj=eqearth")
co <- ne_countries(returnclass = "sf")
co.ee <- st_transform(co, ee)

sla <- st_read("./data/sla.gpkg") |> 
  st_transform(ee)

#sla_rast <- rast("data/sla.tif")
#sla_rast[sla_rast <= 0] <- NA
#sla_rast <- project(sla_rast, "+proj=eqearth")

#(sla_map <- ggplot() +
#  geom_spatraster(data=sla_rast) +
#  scale_fill_viridis_c(expression(Specific~Leaf~Area~"["*m^2~kg^-1*"]"),
#                       option="A",
#                       na.value = NA,
#                       guide = guide_colourbar(
#                         direction = "horizontal",
#                         title.position = "top", 
#                         label.position = "bottom", 
#                         barwidth = 10
 #                      )) +
#  theme_map() +
#  theme(legend.position = c(0.05,0.3) ))

#ggsave("presentation/figures/sla_map.pdf",sla_map, width=12, height=6)


#sla_rast[sla_rast>0] <- 1
#sla_rast <- as.factor(sla_rast)
#writeRaster(sla_rast, "data/sla_mask.tif", overwrite=TRUE)

sla_rast <- rast("data/sla_mask.tif")
  
ppoints <- spatSample(sla_rast, 1000, as.points=TRUE, values=FALSE, na.rm=TRUE) |> 
  st_as_sf() |> 
  st_transform(ee)

#sla <- st_read("./data/nematodes.gpkg") |> 
#  st_transform(ee)

hyperparameter = expand.grid(mtry = 3,
                             splitrule = "variance",
                             min.node.size = 5)

trainDat <- st_drop_geometry(sla) 
r <- c("SLA","Dominant.class", "perc.dominant.class", "Total.class.coverage", "coord.ID")
predictor_names <- c("B5","EVIstd", "B6","EVImax","alt","B2","bio19","B4",
                     "bio4","bio9","B7", "NDWIstd","B3", "bio6", "bio3")
sla_kndm <- knndm(sla, ppoints = ppoints, k=4, maxp=0.8, clustering = "hierarchical")
plot(sla_kndm)

model_spat = train_model("reproduced_knndmcv", sla, predictors = predictor_names,
                         response = "SLA",folds = sla_kndm$clusters,
                         hyperparameter = hyperparameter)
err_spat <- global_validation(model_spat)

clust_rand <- sample(1:10, nrow(sla), replace=TRUE)
model_rand <- train_model("reproduced_random",sla, predictors = predictor_names, 
                          response = "SLA",
                          folds = clust_rand, hyperparameter = hyperparameter)
err_rand <- global_validation(model_rand)


distclust_geo <- function(distm, folds){
  alldist <- rep(NA, length(folds))
  for(f in unique(folds)){
    alldist[f == folds] <- apply(distm[f == folds, f != folds, drop=FALSE], 1, min)
  }
  alldist
}

distmat <- sf::st_distance(sla)
units(distmat) <- NULL
diag(distmat) <- NA
Gj <- apply(distmat, 1, function(x) min(x, na.rm=TRUE))
Gij <- sf::st_distance(ppoints, sla)
units(Gij) <- NULL
Gij <- apply(Gij, 1, min)
Gjstar <- distclust_geo(distmat, clust_rand)
sla_rand <- list(Gj=Gj, Gij=Gij, Gjstar=Gjstar)


ecdf_spat <- plot(sla_kndm) +
  scale_x_continuous(n.breaks = 3, limits=c(0,8*10^6)) +
  theme_bw(base_size=15) +
  theme(legend.position = NaN, aspect.ratio = 1)
  
ecdf_rand <- plot.knndm(sla_rand, overlay = TRUE)+
  scale_x_continuous(n.breaks = 3, limits=c(0,8*10^6)) +
  theme_bw(base_size=15) +
  theme(legend.position = NaN, aspect.ratio = 1)

#legend_ecdfs <- get_legend(ecdf_rand +
#                             guides(color = guide_legend(nrow = 1)) +
#                             theme(legend.position = "bottom"))
#bottom_row <- plot_grid(ecdf_rand,ecdf_spat, nrow=1)
#ecdfs <- plot_grid(bottom_row, NULL, legend_ecdfs, nrow=3,rel_heights = c(2,-0.4,1))
#ecdfs

#ecdfs <- ecdfs +
#  draw_text(paste0("RMSE = ", round(err_spat[[1]],2)), x = 0.78, y = 0.6, hjust = 0, size = 10) +
#  draw_text(paste0("R^2 = ", round(err_spat[[2]],2)), x = 0.78, y = 0.55, hjust = 0, size = 10) +
#  draw_text(paste0("RMSE = ", round(err_rand[[1]],2)), x = 0.28, y = 0.6, hjust = 0, size = 10) +
#  draw_text(paste0("R^2 = ", round(err_rand[[2]],2)), x = 0.28, y = 0.55, hjust = 0, size = 10) 
  

sla <- rbind(sla,sla)
sla$CV <- NaN
sla[1:(nrow(sla)/2),"CV"] <- clust_rand
sla[(nrow(sla)/2+1):nrow(sla),"CV"] <- sla_kndm$clusters
sla[1:(nrow(sla)/2),"what"] <- "random"
sla[(nrow(sla)/2+1):nrow(sla),"what"] <- "kNNDM"
sla$what <- forcats::fct_relevel(sla$what, "random","kNNDM")


sla_map_rand <- ggplot() +
    geom_spatraster(data=sla_rast)+
    scale_fill_manual(values="grey", na.value=NA) +
    geom_sf(data=sla[sla$what=="random",], size=0.5, color=sla[sla$what=="random",]$CV) +
    theme_minimal()+
    #facet_wrap(~what) +
    theme(strip.text.x = element_text(face="bold", size = 14),
          legend.position = NaN, axis.text.x = element_blank(),
          panel.grid = element_blank())


sla_map_kndm <- ggplot() +
    geom_spatraster(data=sla_rast)+
    scale_fill_manual(values="grey", na.value=NA) +
    geom_sf(data=sla[sla$what=="kNNDM",], size=0.5, color=sla[sla$what=="kNNDM",]$CV) +
    theme_minimal() +
    #facet_wrap(~what) +
    theme(strip.text.x = element_text(face="bold", size = 14),
          legend.position = NaN, axis.text.x = element_blank(),
          panel.grid = element_blank())

bp_rand <- plot_grid(sla_map_rand,ecdf_rand,nrow=1, rel_widths = c(1,0.5))
bp_kndm <- plot_grid(sla_map_kndm,ecdf_spat,nrow=1, rel_widths = c(1,0.5))

bp_rand <- bp_rand+
  draw_text(paste0("RMSE = ", round(err_rand[[1]],2)), x = 0.86, y = 0.5, hjust = 0, size = 15) +
  draw_text(paste0("R^2 = ", round(err_rand[[2]],2)), x = 0.86, y = 0.45, hjust = 0, size = 15) 

bp_kndm <- bp_kndm +
  draw_text(paste0("RMSE = ", round(err_spat[[1]],2)), x = 0.86, y = 0.5, hjust = 0, size = 15) +
  draw_text(paste0("R^2 = ", round(err_spat[[2]],2)), x = 0.86, y = 0.45, hjust = 0, size = 15) 


ggsave("presentation/figures/sla_plot_rand.pdf",bp_rand, width = 12, height = 6)
ggsave("presentation/figures/sla_plot_kndm.pdf",bp_kndm, width = 12, height = 6)


