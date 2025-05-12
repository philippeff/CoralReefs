### CREATE VARIABLES FOR FURTHER ANALYSES ###

## Author: Philippe Fernandez-Fournier, 2025

#setwd("D:/sfuvault/SFU/MooersLab/Ch4_CoralReefDiv/Analysis/")
rm(list=ls())

#devtools::install_github("jmadinlab/habtools") # Re-install often, since it is in development

## Load packages
library(ape)
library(phytools)
library(picante)
library(tidyverse)
library(geosphere)
library(habtools) 
library(GeoDist)
library(raster)
library(sf)
library(ggmap)
library(ggsn)
library(ggrepel)
library(missForest)
library(vegan)
library(reshape2)
library(FD)
library(dplyr)
library(ggplot2)

## Fish data
fish.data <- read.csv("./data/data_fish_final.csv")
fish.traits <- read.csv("./data/data_fishtraits_PFF.csv")

## Coral data
coral.data <- read.csv("./data/data_coral_final.csv") # about half of the ID done
coral.cover <- read.csv("./data/data_coral_reefcloud_2024-05-26.csv")

## Site data
transects <- read.csv("./data/data_transects_final.csv")
sites <- read.csv("./data/data_sites_final.csv")
villages <- read.csv("./data/data_villages.csv")
dist.mat <- read.csv("./results/dist_navig.csv", row.names = 1) # Navigable distance in km between reefs (rows) and villages (cols)

## Map data
r4_land <- read_sf("./data/maps/geoBoundaries-IDN-ADM0.shp") 
r4_mpa <- st_read("./data/maps/Selat_Dampier.shp") 
r4_mpa <- st_zm(r4_mpa) # remove z-axis
r4_mpa <- st_transform(r4_mpa, crs=4326) # change CRS to WGS84
r4_mpa <- r4_mpa[-c(2,5),] # Remove zona Inti
r4_mpa$Use_Zone <- c("Commercial", "No take zone", "Use zone") # English names


## Function to normalize a vector
norm_var <- function(var){ return(((var-min(var))/(max(var)-min(var)))) }


# FISHING PRESSURE --------------------------------------------------------

## Create pressure index based on distance that is not linear (sasi!)
# Named it reef accessibility
lnorm.mean <- 2
lnorm.sd   <- 1
lnorm.dens <- dlnorm(seq(0.01, 10, length.out = 1000), meanlog = lnorm.mean, sdlog = lnorm.sd)

# Reef accessibility plot
access.data <- data.frame(x = seq(0, 50, length.out = 1000),
                          y = dlnorm(seq(0, 50, length.out = 1000), meanlog = lnorm.mean, sdlog = lnorm.sd))

access.data$rel_press <- (access.data$y - min(lnorm.dens))/(max(lnorm.dens) - min(lnorm.dens))

# Plot
plot_access <- ggplot(access.data, aes(x = x, y = rel_press)) +
  geom_line(color = "grey30", linewidth = 1) +
  labs(x = "Navigable distance from settlement to reef (km)",
       y = "Relative accessibility of reef") +
  theme_bw(); plot_access

ggsave("../Graphs/fishing_pressure_access.jpeg", plot_access, width = 7, height = 4, dpi = 300)


# Create function to output lognormal density probability of a given vector of distances (i.e. pressure from distance to village)
lnorm_prob_dist <- function(d){
  lnorm.dens <- dlnorm(seq(0.01, 10, length.out = 1000), meanlog = lnorm.mean, sdlog = lnorm.sd)
  out <- dlnorm(d, meanlog = lnorm.mean, sdlog = lnorm.sd)
  out.normed <- (out - min(lnorm.dens))/(max(lnorm.dens) - min(lnorm.dens))
  return(out.normed)
}

## Function for accessibility: index based on combined vicinity of villages
# Searches for lat and lon columns in "reefs" and calculate accessibility index based on "towns"
access_index2 <- function(reefs, towns, sasi = TRUE){
  # reefs <- sites
  # towns <- villages
  
  reefs.lat <- which(colnames(reefs) == "Latitude" | colnames(reefs) == "latitude" | 
                       colnames(reefs) == "Lat" | colnames(reefs) == "lat")
  reefs.lon <- which(colnames(reefs) == "Longitude" | colnames(reefs) == "longitude" | 
                       colnames(reefs) == "Lon" | colnames(reefs) == "lon")
  towns.lat <- which(colnames(towns) == "Latitude" | colnames(towns) == "latitude" | 
                       colnames(towns) == "Lat" | colnames(towns) == "lat")
  towns.lon <- which(colnames(towns) == "Longitude" | colnames(towns) == "longitude" | 
                       colnames(towns) == "Lon" | colnames(towns) == "lon")
  
  ## Calculate reef accessibility for each site
  pressure <- vector()
  for(i in 1:nrow(reefs)){
    ## 1. Pressure from village population
    pop_press <- sqrt(towns$population)
    
    ## 2. Pressure from distance to village
    #d.town <- round(as.vector(distm(reefs[i, c(reefs.lon, reefs.lat)], towns[,c(towns.lon, towns.lat)])/1000), 2)
    d.town <- as.numeric(dist.mat[i,])
    
    if(sasi){
      dist_press <- lnorm_prob_dist(d.town)
      pressure[i] <- sum(pop_press * dist_press)
    } else {
      pressure[i] <- sum(pop_press/(d.town+1))
    }
  }
  return(pressure)
}


(sites$fp.sasi <- norm_var(access_index2(sites, villages, sasi = TRUE))) 
(sites$fp.nosasi <- norm_var(access_index2(sites, villages, sasi = FALSE)))


## Save sites data
write.csv(sites, "./data/data_sites_final.csv", row.names = FALSE)

# Map the indices
register_google(key = "AIzaSyApiB1UjFPZiNfBYzEAka4jiEWHjlKpJcI")
map <- get_googlemap(center = c(lon = mean(range(sites$Longitude)), lat = mean(range(sites$Latitude))), 
                     zoom = 10, maptype = "satellite")

## With sasi
jpeg(filename = "../Graphs/map_sites_FP_sasi.jpeg", 
     res = 300, width = 8, height = 6, units = "in")

ggmap(map) +   
  geom_point(data = villages, aes(x = Longitude, y = Latitude), 
             size = 3, shape = 15, color = "white") +
  geom_point(data = sites, aes(x = Longitude, y = Latitude, color = fp.sasi),
             size = 3, shape = 19) +
  scale_color_gradientn(name = "Fishing\npressure", colors = c("green", "yellow", "orange", "firebrick2")) +
  xlim(c(130.3, 130.85)) +
  ylim(c(-0.70, -0.33)) +
  labs(x = "Longitude", y = "Latitude") +
  ggsn::scalebar(x.min = 130.75, x.max = 130.85,
                 y.min = -0.69,  y.max = -0.43, transform = TRUE,
                 dist = 5, dist_unit = "km", model = 'WGS84',
                 st.size = 2.5, st.color = "white")

dev.off()


# STRUCTURAL COMPLEXITY ---------------------------------------------------

# Metrics based on orthomosaics with package "habtools"

# FD (fractal dimension): captures how much a surface fills volume. Better for fish? 
# RG (rugosity):          captures surface area per unit planar area, i.e. "surface rugosity"

dem.list <- list.files(pattern = "dem_", path = "./DEMs/")

## Structural complexity metrics given a list of DEM file names
sc_dem <- function(dem_path = "./DEMs/",  
                   dem.list = list.files(pattern = "dem_", path = dem_path, full.names = TRUE),
                   n.samp = 100) {
  #dem.list <- list.files(pattern = "dem_", path = dem_path, full.names = TRUE)
  
  tr      <- vector()
  dem.rg  <- vector()
  dem.fd  <- vector()
  dem.res <- vector()
  
  for(i in seq_along(dem.list)) {
    dem.i <- raster(dem.list[i])
    
    # Set CRS
    crs(dem.i) <- "+proj=tmerc +lat_0=0 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    
    # Extract transect name
    tr[i] <- sapply(strsplit(names(dem.i), split = "_"), function(x) x[2])
    
    dem.rg.samp <- vector()
    dem.fd.samp <- vector()
    
    for(j in 1:n.samp) {
      dem.i.samp <- dem_sample(dem.i, L = 1, allow_NA = FALSE, plot = (j == 1), max_iter = 1000)
      dem.rg.samp[j] <- rg(dem.i.samp, method = "area", L0 = 0.01)
      dem.fd.samp[j] <- fd(dem.i.samp, method = "area", lvec = c(0.01, 0.05, 0.15, 0.3, 0.6), plot = FALSE, diagnose = FALSE)
    }
    
    dem.rg[i]  <- mean(dem.rg.samp, na.rm = TRUE)
    dem.fd[i]  <- mean(dem.fd.samp, na.rm = TRUE)
    dem.res[i] <- res(dem.i)[1]
  }
  
  SC.transects <- data.frame(transect = tr, rugosity = dem.rg, fractal.d = dem.fd, dem.res = dem.res)
  return(SC.transects)
}

SC.transects <- sc_dem()
write.csv(SC.transects, "./results/SC_transects.csv", row.names = FALSE)


# FISH TRAITS IMPUTATION ----------------------------------------------

# Are all species in fish.data found in fish.traits?
all(unique(fish.data$Genus_species) %in% fish.traits$Genus_species)

## trait data to calculate FD metrics: species as row names and traits as columns
fish.traits.matrix <- fish.traits %>%
  dplyr::select(Genus_species, Size, Mobility, Activity, Schooling, Position, Diet_Parravicini_2020) %>%
  filter(Genus_species %in% fish.data$Genus_species) %>%
  filter(!is.na(Diet_Parravicini_2020)) %>%
  arrange(Genus_species) %>%
  column_to_rownames("Genus_species") 

# Convert variables to factors (ordinal & categorical)
fish.traits.matrix$Size <- ordered(fish.traits.matrix$Size)
fish.traits.matrix$Mobility <- ordered(fish.traits.matrix$Mobility)
fish.traits.matrix$Activity <- ordered(fish.traits.matrix$Activity)
fish.traits.matrix$Schooling <- ordered(fish.traits.matrix$Schooling)
fish.traits.matrix$Position <- ordered(fish.traits.matrix$Position)
fish.traits.matrix$Diet_Parravicini_2020 <- factor(fish.traits.matrix$Diet_Parravicini_2020)

#write.csv(fish.traits.matrix, "./results/matrix_traits_fish.csv", row.names = TRUE)
#fish.traits.matrix <- read.csv("./results/matrix_traits_fish.csv", row.names = 1)

## Species occurrence matrix: transects as rows and species as columns
fish.occ.matrix <- dcast(fish.data, Transect ~ Genus_species, value.var = "Total", fun.aggregate = sum, fill = 0) %>%
  column_to_rownames("Transect")

# Remove species with full COLUMN of NAs
na.species <- fish.traits$Genus_species[is.na(fish.traits$Diet_Parravicini_2020)]
fish.occ.matrix <- fish.occ.matrix[, -which(colnames(fish.occ.matrix) %in% na.species)]

#write.csv(fish.occ.matrix, "./results/matrix_occurrence_fish.csv", row.names = TRUE)
#fish.occ.matrix <- read.csv("./results/matrix_occurrence_fish.csv", row.names = 1)

# Check if same order of species
all(colnames(fish.occ.matrix) == rownames(fish.traits.matrix))

## Deal with NAs in traits matrix (dbFD function does not allow NAs): use missForest without the phylogeny
#traits.genus <- unlist(lapply(strsplit(row.names(fish.traits.matrix), "_"), function(x) x[1]))

# Calculate proportion of NAs for each column (very few)
na_prop <- sapply(fish.traits.matrix, function(col) mean(is.na(col)))

ftm_imp <- missForest(fish.traits.matrix, 
                      maxiter = 100, ntree = 1000, 
                      variablewise = TRUE, verbose = FALSE)

fish.traits.matrix_imp <- ftm_imp$ximp

# Check error
ftm_imp$OOBerror

# Save imputation stats

imp_table <- data.frame(trait = names(na_prop),
                        prop_NA = as.numeric(na_prop),
                        imp_err = ftm_imp$OOBerror)

write.csv(imp_table, file = "./results/fishtraits_imputation_err.csv", row.names = FALSE)

# FISH FUNCTIONAL DIVERSITY ----------------------------------------------
## Calculate Functional Diversity indices (presence/absence)

# FRic: (Functional Richness) is measured as the number of unique trait combinations when only categorical and/or ordinal traits are present (NOT as the convex hull volume). 
# FDiv: (Functional divergence) quantifies the degree to which the abundances of species in a community 
#       are distributed towards the boundaries of the occupied functional space.
# FDis: (Functional dispersion) is estimated as the mean distance of all species to the weighted centroid of the community in the trait space.
# RaoQ: (Rao’s quadratic entropy) is a measure of diversity of ecological communities defined by Rao (1982) and is based on 
#       the proportion of the abundance of species present in a community and some measure of dissimilarity among them.
# FGR:  (a posteriori functional group richness) is  the number of functional groups. 
# CWM:  (community-level weighted means of trait values) is an index of functional composition.

# FEs:  (Functional Entities) are unique combinations of trait values (Mouillot et al. 2014). Measured as number of FEs, i.e. same as FRic on categorical-only traits...

# Hellinger-transformed matrix
fish.occ.matrix_hell <- decostand(fish.occ.matrix, method = "hellinger")

fish.fd <- dbFD(fish.traits.matrix_imp, fish.occ.matrix_hell, 
                w.abun = TRUE, ord = "podani", corr = "none",
                m = 5,
                calc.FRic = TRUE,
                calc.CWM = TRUE, CWM.type = "dom",
                calc.FGR = FALSE, 
                calc.FDiv = FALSE,
                print.pco = TRUE)

fish.fd$FRic
range(fish.fd$FRic)
hist(fish.fd$FRic, breaks = 25)

fish.fd$FDis
range(fish.fd$FDis)
hist(fish.fd$FDis, breaks = 25)

head(fish.fd$CWM)

fish.fd$RaoQ
range(fish.fd$RaoQ)
hist(fish.fd$RaoQ, breaks = 25)

# Calculate the proportion of variance explained by each axis
var_exp <- fish.fd$x.values/sum(fish.fd$x.values)

# Extract the proportion of variance explained for PCoA1 and PCoA2
pcoa1_var <- round(var_exp[1] * 100, 2)
pcoa2_var <- round(var_exp[2] * 100, 2)

# Create a data frame with PCoA results
pcoa_df <- as.data.frame(fish.fd$x.axes) 

# Rename the columns for convenience
colnames(pcoa_df) <- gsub("A", "PCoA", colnames(pcoa_df))

# Add the species information to the PCoA data
pcoa_df$species_name <- rownames(fish.traits.matrix_imp) 
pcoa_df$Diet <- fish.traits.matrix_imp$Diet_Parravicini_2020

# Replace combined category with "omnivore"
pcoa_df$Diet <- gsub("herbivores microvores detritivores", "omnivore", pcoa_df$Diet)

# Plot the PCoA results
# Define color palette
diet_colors <- c("corallivores" = "#FF7F50",
                 "crustacivores" = "#DC143C",
                 "omnivore" = "gold2",
                 "macroinvertivores" = "#008080",
                 "microinvertivores" = "#00CED1",
                 "piscivores" = "#1E90FF",
                 "planktivores" = "#9370DB",
                 "sessile invertivores" = "#FF00FF")

plot_pcoa_fish <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Diet)) +
  geom_point(size = 2) +  
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(title = "PCoA of fish traits (Gower distance)",
       x = paste("PCoA 1 (", pcoa1_var, "%)", sep = ""),
       y = paste("PCoA 2 (", pcoa2_var, "%)", sep = "")) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_color_manual(values = diet_colors); plot_pcoa_fish


ggsave(filename = "../Graphs/fish_PCoA_FD.jpeg", plot = plot_pcoa_fish,
       width = 7.5, height = 5, units = "in", dpi = 300)


# COMPILE FISH DIVERSITY METRICS ------------------------------------------
all(row.names(fish.occ.matrix) == fish.pd$transect) # same order

fish.div <- data.frame(transect = row.names(fish.occ.matrix),
                       fish.FRic     = fish.fd$FRic,
                       fish.FDis     = fish.fd$FDis,
                       fish.RaoQ     = fish.fd$RaoQ,
                       fish.PD       = fish.pd$PD,
                       fish.SR       = fish.pd$SR)

write.csv(fish.div, "./data/fish_div.csv", row.names = FALSE)


# RELATIVE CORAL COVER METRICS --------------------------------------------------
# Combine quadrat measurements of percent cover to transect-level and site-level
colnames(coral.cover)

coral.cover <- coral.cover %>%
  group_by(survey_title) %>% # group by transect
  dplyr::select(c(11,25,40,41))

# Combine Mikes point
coral.cover$survey_title[which(coral.cover$survey_title %in% c("MKP02", "MKP03"))] <- "MKP04"

# Replace non-annotated human points by machine points
coral.cover$point_human_classification[coral.cover$point_human_classification == ""] <- coral.cover$point_machine_classification[coral.cover$point_human_classification == ""]
coral.cover <- coral.cover[-which(coral.cover$point_human_classification == ""),]

point.types <- unique(coral.cover$point_human_classification)

df <- data.frame(matrix(ncol = length(point.types)+1, nrow = 0))
colnames(df) <- c("transect", point.types)

for(i in 1:length(unique(coral.cover$survey_title))){
  points.i <- table(coral.cover$point_human_classification[coral.cover$survey_title == unique(coral.cover$survey_title)[i]])
  df[i, names(points.i)] <- points.i
  df$transect[i] <- unique(coral.cover$survey_title)[i]
}

df[is.na(df)] <- 0

df <- df %>%
  mutate_if(is.numeric, ~ . / rowSums(dplyr::select(df, where(is.numeric))))

coral.cover.tr <- df

write.csv(coral.cover.tr, "./data/data_coralcover.csv", row.names = F)



#  COMBINE ALL VARIABLES --------------------------------------------------

fish.div <- read.csv("./data/fish_div.csv")
coral.div <- read.csv("./data/coral_div.csv")
coral.data <- read.csv("./data/data_coral_final.csv") 
coral.cover <- read.csv("./data/data_coralcover.csv")
SC.transects <- read.csv("./results/SC_transects.csv")
transects <- read.csv("./data/data_transects_final.csv")
sites <- read.csv("./data/data_sites_final.csv")

# rename columns of coral.cover
colnames(coral.cover)[2:ncol(coral.cover)] <- paste0("cover_", colnames(coral.cover)[2:ncol(coral.cover)])

# Sub transects
tr_sub <- transects %>%
  rename(transect = Transect.code) %>%
  dplyr::select(transect, Site.name, zone, Latitude, Longitude, Diver, Buddy, Notes)

# Sub sites
sites_sub <- sites %>%
  dplyr::select(Site.name, Dive.site, Sasi, fp.sasi, fp.nosasi)

# Match site names?
all(transects$Site.name %in% sites$Site.name)

# Combine
data_reef <- left_join(tr_sub, SC.transects, by = join_by(transect))
data_reef <- left_join(data_reef, sites_sub, by = join_by(Site.name))
data_reef <- left_join(data_reef, fish.div, by = join_by(transect))
data_reef <- left_join(data_reef, coral.div, by = join_by(transect))
data_reef <- left_join(data_reef, coral.cover, by = join_by(transect))

# Dead and live cover
data_reef$cover_live.coral <- data_reef$cover_SC + data_reef$cover_HC # Live coral is soft + hard coral
data_reef$cover_dead.coral <- data_reef$cover_RB + data_reef$cover_DC# Dead coral is rubble + dead coral

# Simpson and Shannon indices for fish
simpson_ind <- diversity(fish_comm, index = "simpson") # (1-D) 0 = no div, 1 = all spp equal div
names(simpson_ind) <- rownames(fish_comm)
hist(simpson_ind)

shannon_ind <- diversity(fish_comm, index = "shannon") 
names(shannon_ind) <- rownames(fish_comm)
hist(shannon_ind)

# Add to data_reef
data_reef$fish_simpson <- simpson_ind[match(data_reef$transect, names(simpson_ind))]
data_reef$fish_shannon <- shannon_ind[match(data_reef$transect, names(shannon_ind))]

# Combine DEM complexity for Mikes point
data_reef$rugosity[data_reef$transect == "MKP04"] <- mean(SC.transects$rugosity[SC.transects$transect %in% c("MKP02", "MKP03")])
data_reef$fractal.d[data_reef$transect == "MKP04"] <- mean(SC.transects$fractal.d[SC.transects$transect %in% c("MKP02", "MKP03")])

# Coral.SR for MKP04
data_reef$coral.SR[data_reef$transect == "MKP04"] <- length(unique(coral.data$Genus_morpho[coral.data$Transect == "MKP04"]))

# Recalculate structural complexity for transects with missing DEMs (based on other transect...)
dem_all <- list.files(pattern = "dem_", path = dem_path, full.names = TRUE)
tr_redo <- c("AGU01", # AGU02 missing
             "MER01", # MER02 missing
             "COV02", # COV01 missing
             "PAP01", # PAP02 missing
             "YEB01") # YEB02 missing
dem_miss <- grep(paste(tr_redo, collapse = "|"), dem_all, value = TRUE)

new_sc <- sc_dem(dem_path = "./DEMs/", dem.list = dem_miss, n.samp = 100) # takes a while

data_reef$rugosity[data_reef$transect == "AGU02"] <- new_sc$rugosity[new_sc$transect == "AGU01"]
data_reef$fractal.d[data_reef$transect == "AGU02"] <- new_sc$fractal.d[new_sc$transect == "AGU01"]

data_reef$rugosity[data_reef$transect == "MER02"] <- new_sc$rugosity[new_sc$transect == "MER01"]
data_reef$fractal.d[data_reef$transect == "MER02"] <- new_sc$fractal.d[new_sc$transect == "MER01"]

data_reef$rugosity[data_reef$transect == "COV01"] <- new_sc$rugosity[new_sc$transect == "COV02"]
data_reef$fractal.d[data_reef$transect == "COV01"] <- new_sc$fractal.d[new_sc$transect == "COV02"]

data_reef$rugosity[data_reef$transect == "PAP02"] <- new_sc$rugosity[new_sc$transect == "PAP01"]
data_reef$fractal.d[data_reef$transect == "PAP02"] <- new_sc$fractal.d[new_sc$transect == "PAP01"]

data_reef$rugosity[data_reef$transect == "YEB02"] <- new_sc$rugosity[new_sc$transect == "YEB01"]
data_reef$fractal.d[data_reef$transect == "YEB02"] <- new_sc$fractal.d[new_sc$transect == "YEB01"]

# Add column indicating whether SC was imputed
tr_miss <- c("AGU02", "MER02", "COV01", "PAP02", "YEB02")
data_reef$sc_imputed <- ifelse(data_reef$transect %in% tr_miss, 1, 0)

# GAB coral SR based on visual inspection of quadrats (no closeup taken...)
data_reef$coral.SR[data_reef$transect == "GAB01"] <- 10
data_reef$coral.SR[data_reef$transect == "GAB02"] <- 13

# AGU02 and MER02 coral SR copied from other transect...
data_reef$coral.SR[data_reef$transect == "AGU02"] <- data_reef$coral.SR[data_reef$transect == "AGU01"]
data_reef$coral.SR[data_reef$transect == "MER02"] <- data_reef$coral.SR[data_reef$transect == "MER01"]

# AGU02 and MER02 dead coral cover copied from other transect...
data_reef$cover_dead.coral[data_reef$transect == "AGU02"] <- data_reef$cover_dead.coral[data_reef$transect == "AGU01"]
data_reef$cover_dead.coral[data_reef$transect == "MER02"] <- data_reef$cover_dead.coral[data_reef$transect == "MER01"]

# Save
write.csv(data_reef, "./data/data_reef_final.csv", row.names = FALSE)


# SAVE FISH COMM MATRIX ---------------------------------------------------

fish_comm <- read.csv("./results/matrix_occurrence_fish.csv", row.names = 1)
fish_comm <- fish_comm[rownames(fish_comm) %in% data_reef$transect, ] # only transects from data_reef
fish_comm <- fish_comm[, colSums(fish_comm) > 0] # remove species with zero obs
write.csv(fish_comm, "./results/matrix_comm_fish.csv", row.names = TRUE)
