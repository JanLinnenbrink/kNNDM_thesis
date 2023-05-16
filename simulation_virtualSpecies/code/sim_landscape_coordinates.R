#-----------------------------------------------------------#
#             virtualspecies landscape generation           #
#-----------------------------------------------------------#

# Libraries and utils ----
library("tidyverse")
library("rnaturalearth")
#library("rasterVis")
library("raster")
library("sf")
library("virtualspecies")
library("terra")

# No need for proj4 warnings
options("rgdal_show_exportToProj4_warnings"="none")

# AOI: Iberian peninsula ----
AOI <- ne_countries(country=c("Spain", "Portugal"), scale="large",
                    returnclass="sf") %>%
  summarise() %>%
  st_cast("POLYGON") %>%
  mutate(area=st_area(.)) %>%
  dplyr::filter(area==max(.$area)) %>%
  summarise()

# Fetch climate data ----
wclim <- raster::getData('worldclim', var='bio', res=2.5)
wclim <- mask(crop(wclim, AOI), AOI)
# levelplot(stretch(wclim))

# add coordinates
template_ll <- wclim$bio1 # Kopie erstellen
coords <- coordinates(template_ll) # Koordinaten für jedes Pixel als tabelle bekommen
lat <- template_ll # Kopie erstellen
lon <- template_ll# Kopie erstellen
values(lat) <- coords[,2] # Koordinaten auf das Raster schreiben
values(lon) <- coords[,1] # Koordinaten auf das Raster schreiben
coords <- stack(lat,lon) # neuen Steck erstellen
names(coords) <- c("lat","lon") #...und benennen
coords <- resample(coords,wclim) # und auf die gleiche Geometrie
#plot(coords)
wclim <- stack(wclim,coords) # alles zusammenführen
# Simulate virtual species ----
set.seed(1234)
meansPCA <- c(0, 0)
sdPCA <- c(2, 2)
covs <- c("bio1", "bio2", "bio4", "bio5", "bio6",
          "bio12", "bio13", "bio14", "bio15")
vspecies <- generateSpFromPCA(subset(wclim, covs),
                              means = meansPCA, sds = sdPCA, plot=T)
# png("figures/vspecies_pca.png", width=2000, height=1500, res=300)
plotResponse(vspecies)
# dev.off()
outcome <- vspecies$suitab.raster
names(outcome) <- "outcome"
# levelplot(outcome, margin = FALSE)
wclim <- stack(wclim, outcome)

# Project to EPSG:25830 (ETRS89 / UTM zone 30N) ----
AOI <- st_transform(AOI, crs = 25830)
spoly <- st_buffer(AOI, -5000) # Shrink boundaries to avoid NAs
wclim <- rast(wclim)
wclim <- project(wclim, "epsg:25830", res = 5000)
wclim <- mask(wclim, AOI)
# levelplot(stretch(wclim))
wgrid <- st_as_sf(rasterToPoints(stack(wclim), spatial = TRUE)) # derive grid
wgrid <- st_transform(wgrid, crs=st_crs(spoly))


# Export elements for simulation
st_write(spoly, dsn="data/species_vdata.gpkg", layer="sampling_polygon")
st_write(wgrid, dsn="data/species_vdata_coords.gpkg", layer="landscape_grid")
terra::writeRaster(wclim,  "data/species_stack_coords.tif")


