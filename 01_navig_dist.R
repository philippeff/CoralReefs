### CALCULATE NAVIGABLE DISTANCE BY WATER ###

## Author: Philippe Fernandez-Fournier, 2025

#setwd("D:/sfuvault/SFU/MooersLab/Ch4_CoralReefDiv/Analysis/")
rm(list=ls())

library(raster)
library(sf)
library(gdistance)
library(maptools)
library(rgdal)
library(ggplot2)
library(spdep)

## Data
transects <- read.csv("./data/data_transects_final.csv")
sites <- read.csv("./data/data_sites_final.csv")
villages <- read.csv("./data/data_villages.csv")

ocean <- shapefile("./data/maps/geoBoundaries-IDN-ADM0.shp") # download here: https://www.geoboundaries.org/countryDownloads.html

## Parameters
res1 <- 0.001
site.names <- sites$Site.name
site.xy <- sites[,c(6,5)]

xmin <- 130.2
xmax <- 130.9
ymin <- -0.9
ymax <- -0.2

# CREATE RASTER -----------------------------------------------------------

ocean <- raster::crop(ocean, y = raster::extent(xmin, xmax, ymin, ymax)) 
region.as.table <- matrix(NA, nrow = ((ymax-ymin)/res1), ncol = ((xmax-xmin)/res1))
region.as.raster <- raster(region.as.table)
extent(region.as.raster) <- c(xmin,xmax,ymin,ymax)

ocean.r <- rasterize(ocean, region.as.raster)
ocean.r[!is.na(ocean.r)] <- 0
ocean.r[is.na(ocean.r)] <- 1


# RELOCATE VILLAGES INTO WATER IF NEEDED ---------------------------------------------

sites.to.relocate <- extract(ocean.r, villages[, c(3,2)]) == 0
sites.to.relocate.xy <- villages[sites.to.relocate,c(3,2)]

if(nrow(sites.to.relocate.xy) > 0 ) {

  ocean.r.sites <- as.data.frame(ocean.r, xy=TRUE)[,1:2]
  ocean.r.sites <- cbind(ocean.r.sites,extract(ocean.r,ocean.r.sites))
  ocean.r.sites <- ocean.r.sites[ocean.r.sites[,3] == 1 , 1:2 ]
  
  for (k in 1:nrow(sites.to.relocate.xy)) {
    near.cells <- as.data.frame(sort(spDistsN1(as.matrix(ocean.r.sites), as.matrix(sites.to.relocate.xy[k,1:2]), longlat = TRUE), 
                                       decreasing = FALSE, index.return = TRUE)) 
    villages[which(sites.to.relocate)[k], c(3,2)] <- ocean.r.sites[ near.cells[1,2] , 1:2 ]
  }
}

# West Gam reef is on land, jitter to the west
sites[8,]
sites$Longitude[8] <- sites$Longitude[8] -0.001

# CALCULATE DISTANCES -----------------------------------------------------

cost.surface <- ocean.r
projection(cost.surface) <- CRS("+proj=longlat +datum=WGS84")

raster_tr <- transition(cost.surface, mean, directions=16)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

# Create points
dist.mat <- matrix(nrow = nrow(sites), ncol = nrow(villages))
row.names(dist.mat) <- sites$Site.name
colnames(dist.mat)  <- villages$village

for(i in 1:nrow(sites)){
  
  # tracker
  cat("\r", paste0("Progress: ", i, "/", nrow(sites), ' reefs\n'))
  
  # Create coordinates
  from_coords <- as.matrix(sites[i, c(6,5)])
  to_coords <- as.matrix(villages[, c(3,2)])

  # Calculate distances (first reef to all villages)
  raster_cosDist <- costDistance(raster_tr_corrected, fromCoords = from_coords, toCoords = to_coords)
  distance.marine <- as.matrix(raster_cosDist)/1000

  dist.mat[i,] <- distance.marine
  
}

## Check: is any distance not calculated?
any(is.na(dist.mat))

## Check: Any zero distances? Convert to 10m
dist.mat[dist.mat == 0] <- 0.01


## Save the distance matrix
write.csv(dist.mat, "./results/dist_navig.csv", row.names = TRUE)

# Plot polygons, points and lines
i <- 8
plot(ocean.r , col = c("#fffafa","lightblue"), legend = FALSE)
points(from_coords, pch = 20, col = "blue", cex = 2)
points(to_coords, pch = 20, col = "red", cex = 2)

# Plot all lines (takes a while)
for(ii in 1:nrow(villages)){
  lines(shortestPath(raster_tr_corrected, from_coords, to_coords[ii,], output="SpatialLines"))
}


