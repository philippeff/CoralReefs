### CALCULATE SITE-LEVEL WAVE EXPOSURE ###

## Author: Philippe Fernandez-Fournier, 2024

#setwd("D:/sfuvault/SFU/MooersLab/Ch4_CoralReefDiv/Analysis/")
rm(list=ls())

## Load packages
library(sf)
library(ggplot2)
library(ggrepel)
library(ggmap)
library(terra)
library(raster)
library(waver)
# Try also stagg: https://github.com/tcarleton/stagg

## Load data
#data_reef <- read.csv("D:/sfuvault/SFU/MooersLab/Ch4_CoralReefDiv/Analysis/data/data_reef_final.csv")
land <- shapefile("./data/maps/geoBoundaries-IDN-ADM0.shp")
sites <- read.csv("./data/data_sites_final.csv")

# Wind data from the Global Wind Atlas (GWA)
# https://globalwindatlas.info/
# Wind frequency rose: relative wind frequency from each direction 
# At each site


# WAVE EXPOSURE PER SITE -----------------------------------------------------------

# All SITES coordinates
coords <- data.frame(lon = sites$Longitude, 
                     lat = sites$Latitude)

coords_sf <- st_as_sf(coords, 
                coords = c("lon", "lat"), 
                crs = 4326)

# Bearings (every 30 degrees)
brgs <- seq(0, 330, by = 30)

# Shoreline as sf object
land_sf <- st_as_sf(land)
shoreline_sf <- st_cast(land_sf, to = "MULTILINESTRING") # Needs this format for wind fetch

# Calculate wind fetch lengths
fetch_lengths <- fetch_len_multi(pts = coords_sf,      # coordinates of sites
                                 bearings = brgs,      # Wind directions
                                 shoreline = shoreline_sf, # sf of shoreline
                                 dmax = 100000,        # max fetch length if no land (100km)
                                 spread = 0)           # Vector of relative bearings for each main bearing, e.g. c(-10, 0, 10)

# Average wind fetch for each site (use wind frequency rose by ZONE for weighted mean)

path_gwa <- "./data/Wind_GWA/sites/"

wind_fetch_avg <- vector()

for (i in 1:nrow(fetch_lengths)) {
  windrose_i <- read.csv(paste0(path_gwa, "windFrequencyRose_", sites$Site.code[i], ".csv"))
  
  # direction weights for weighted mean relative to max wind frequency
  weights_i <- windrose_i$value / max(windrose_i$value)

  wind_fetch_avg[i] <- mean(fetch_lengths[i, ] * weights_i)
}



# Calculate wave exposure as wave flux in power per meter of wave crest (kW/m)

# Wind speed (m/s) of windiest 10% of area for each site (manually from GWA website)
wind_speed <- c(3.05, 3.02, 3.00, 3.04, 2.76, 2.72, 3.10, 2.75, 2.84, 2.82, 
                2.65, 2.67, 2.47, 2.66, 2.50, 2.69, 2.42, 2.41, 2.41, 2.51, 
                2.31, 2.77, 2.28, 2.34, 3.18, 3.08, 2.44, 2.94, 2.70, 2.48, 
                2.41, 2.50, 2.38, 2.39, 2.86, 2.50, 2.56, 2.39, 2.49, 2.43)
names(wind_speed) <- sites$Site.code

hist(wind_speed, breaks = 20)

wave_exp <- vector()

for(ii in 1:nrow(fetch_lengths)){
  
  # direction weights for weighted mean relative to max wind frequency
  wave_exp[ii] <- wave_energy(wind = wind_speed[ii], 
                              fetch = wind_fetch_avg[ii], 
                              depth = 10)   # Water depth (in meters)
}

range(wave_exp)

# West Mansuar is outlier?
hist(wave_exp, breaks = 20)
hist(log(wave_exp), breaks = 20) # better
hist(sqrt(wave_exp), breaks = 20) # best

hist(sites$wave_exp_norm2, breaks = 20)
hist(log(sites$wave_exp_norm2), breaks = 20) # better
hist(sqrt(sites$wave_exp_norm2), breaks = 20) # best

# ADD TO SITES ------------------------------------------------------------

# Add wave exposure
sites$wave_exposure2 <- wave_exp

# Normalized wave exposure
sites$wave_exp_norm2 <- (sites$wave_exposure2 - min(sites$wave_exposure2)) / (max(sites$wave_exposure2) - min(sites$wave_exposure2))

# Compare to other method
plot(sites$wave_exp_norm, sites$wave_exp_norm2)
abline(a = 0, b = 1, col = "red")

# Save
write.csv(sites, "./data/data_sites_final.csv", row.names = FALSE)


# PLOT --------------------------------------------------------------------

# Plot each site and its wave exposure
register_google(key = "AIzaSyApiB1UjFPZiNfBYzEAka4jiEWHjlKpJcI")
map <- get_googlemap(center = c(lon = mean(range(sites$Longitude)), lat = mean(range(sites$Latitude))), 
                     zoom = 10, maptype = "satellite")

jpeg(filename = "../Graphs/map_wave_exposure_bysite.jpeg", 
     res = 300, width = 8, height = 6, units = "in")

ggmap(map) +   
  geom_point(data = sites, aes(x = Longitude, y = Latitude, color = wave_exp_norm2),
             size = 3, shape = 19) +
  scale_color_gradientn(name = "Wave\nexposure", colors = c("green", "greenyellow", "yellow", "orange", "firebrick2")) +
  xlim(c(130.3, 130.85)) +
  ylim(c(-0.70, -0.33)) +
  labs(x = "Longitude", y = "Latitude") +
  ggsn::scalebar(x.min = 130.75, x.max = 130.85,
                 y.min = -0.69,  y.max = -0.43, transform = TRUE,
                 dist = 5, dist_unit = "km", model = 'WGS84',
                 st.size = 2.5, st.color = "white")

dev.off()



