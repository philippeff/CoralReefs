### MAIN ANALYSES AND PLOTS ###

## Author: Philippe Fernandez-Fournier, 2025

#setwd("D:/sfuvault/SFU/MooersLab/Ch4_CoralReefDiv/Analysis/")
rm(list=ls())

## Load packages
library(sf)
library(ggplot2)
library(ggmap)
library(rnaturalearth)
library(raster)
library(tidyr)
library(dplyr)
library(tibble)
library(ggrepel)
library(GGally)
library(ggridges)
library(ggforce)
library(ggpp)
library(ggConvexHull)
library(colorRamps)
library(patchwork)
library(vegan)
library(lme4)
library(brms)
library(piecewiseSEM)


## Load transect-level data (all combined)
data_reef <- read.csv("./data/data_reef_final.csv")

## Fish data
fish.data <- read.csv("./data/data_fish_final.csv")
fish.species <- read.csv("./data/meta_fish_species.csv")
fish_comm <- read.csv("./results/matrix_comm_fish.csv", row.names = 1)

## Coral data
coral.traits <- read.csv("./results/matrix_traits_coral.csv", row.names = 1)
coral.data <- read.csv("./data/data_coral_final.csv")

## Site data
transects <- read.csv("./data/data_transects_final.csv")
sites <- read.csv("./data/data_sites_final.csv")
villages <- read.csv("./data/data_villages.csv")

## Map data
r4_land <- read_sf("./data/maps/geoBoundaries-IDN-ADM0.shp") 
r4_mpa <- st_read("./data/maps/Selat_Dampier.shp") 
r4_mpa <- st_zm(r4_mpa) # remove z-axis
r4_mpa <- st_transform(r4_mpa, crs=4326) # change CRS to WGS84
r4_mpa <- r4_mpa[-c(2,5),] # Remove zona Inti
r4_mpa$Use_Zone <- c("Commercial", "No-take zone", "Use zone") # English names


#  NEW DATASETS -----------------------------------------------------------

## Create data for SEMs
data_nona <- data_reef %>%
  filter(!is.na(coral.SR) & !is.na(fractal.d) & !is.na(cover_dead.coral)) 

# convert zone to binary (1 = NTZ, 0 = UZ)
data_nona$zone_bin <- ifelse(data_nona$zone == "NTZ", 1, 0) 

# Add wave exposure (site-level GWA data)
data_nona$wave_exp <- sites$wave_exp_norm2[match(data_nona$Site.name, sites$Site.name)]

# Transform 0:1 vars (Arcsine square-root)
data_nona <- data_nona %>% 
  mutate(fp.sasi_transf         = asin(sqrt(fp.sasi)),          # relative fishing pressure (0-1)
         cover_deadcoral_transf = asin(sqrt(cover_dead.coral)), # relative dead coral cover (0-1)
         cover_livecoral_transf = asin(sqrt(cover_live.coral)), # relative live coral cover (0-1)
         cover_malgae_transf    = asin(sqrt(cover_MA)),         # relative macroalgae cover (0-1)
         cover_sponge_transf    = asin(sqrt(cover_SP)),         # relative sponge cover (0-1)
         fish_simpson_transf    = asin(sqrt(fish_simpson)),     # Simpson index
         wave_exp_transf        = asin(sqrt(wave_exp)))         # relative dead coral cover (0-1)

# Standardize all vars
data_nona <- data_nona %>% 
  mutate(fish_SR          = as.numeric(scale(fish.SR)),
         fish_FRic        = as.numeric(scale(fish.FRic)),
         fish_FDis        = as.numeric(scale(fish.FDis)),
         fish_RaoQ        = as.numeric(scale(fish.RaoQ)),
         coral_SR         = as.numeric(scale(coral.SR)),
         fractal_D        = as.numeric(scale(fractal.d)),
         fishing_pressure = as.numeric(scale(fp.sasi_transf)),
         wave_exp         = as.numeric(scale(wave_exp_transf)),
         cover_deadcoral  = as.numeric(scale(cover_deadcoral_transf)),
         cover_malgae     = as.numeric(scale(cover_malgae_transf)),
         cover_sponge     = as.numeric(scale(cover_sponge_transf)),
         fish_simpson     = as.numeric(scale(fish_simpson_transf)),
         fish_shannon     = as.numeric(scale(fish_shannon)))

# Add PC1 fish species comp
data_nona$fish_comp_PC1 <- as.numeric(scale(data_reef$fish_comp_PC1))


# DATA STATS --------------------------------------------------------------

# How many individual fish?
sum(fish_comm)

# How many species of fish?
ncol(fish_comm)

# How many individual corals?
nrow(coral.data)

# How many genera of coral?
length(unique(coral.data$Genus))

# Histograms of variables used in SEMs
data_nona_num <- select(data_nona, where(is.numeric))
data_nona_num <- data_nona_num %>%
  select(fish_SR, fish_FDis, fish_comp_PC1, 
         coral_SR, fractal_D, cover_deadcoral, 
         fishing_pressure, zone_bin, wave_exp)

vars_gg <- pivot_longer(data_nona_num, cols = everything(), names_to = "variable", values_to = "value")

plot_semvars <- ggplot(vars_gg, aes(x = value)) +
  geom_histogram(bins = 30, fill = "grey40", color = "black", alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  theme_test() +
  labs(x = "Value", y = "Count", title = "Scaled variables used in SEMs"); plot_semvars

ggsave("../Graphs/SEM_variables_hist.jpeg", plot_semvars, width = 9, height = 9.5, dpi = 300)



# MAP OF SITES & ZONES ------------------------------------------------------------

# Plot each site and its wave exposure
register_google(key = "AIzaSyApiB1UjFPZiNfBYzEAka4jiEWHjlKpJcI")
map <- get_googlemap(center = c(lon = mean(range(sites$Longitude)), lat = mean(range(sites$Latitude))), 
                     zoom = 10, maptype = "satellite")

# Map of sites only
jpeg(filename = "../Graphs/map_sites.jpeg", 
     res = 300, width = 8, height = 6, units = "in")

ggmap(map) +   
  geom_point(data = sites, aes(x = Longitude, y = Latitude),
             color = "green", size = 3, shape = 19) +
  xlim(c(130.3, 130.85)) +
  ylim(c(-0.70, -0.33)) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel(data = sites, aes(x = Longitude, y = Latitude, label = Site.name),
                  size = 3, color = "white", box.padding = 0.1, 
                  point.padding = 0.05, max.overlaps = Inf)
dev.off()

# Map of sites with codes
jpeg(filename = "../Graphs/map_sites_codes.jpeg", 
     res = 300, width = 8, height = 6, units = "in")

ggmap(map) +   
  geom_point(data = sites, aes(x = Longitude, y = Latitude),
             color = "green", size = 3, shape = 19) +
  xlim(c(130.3, 130.85)) +
  ylim(c(-0.70, -0.33)) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel(data = sites, aes(x = Longitude, y = Latitude, label = Site.code),
                  size = 3, color = "white", box.padding = 0.1, 
                  point.padding = 0.05, max.overlaps = Inf)
dev.off()

# Map of sites and whether they are dive sites
sites$Dive.site <- factor(sites$Dive.site)

jpeg(filename = "../Graphs/map_sites_divesite.jpeg", 
     res = 300, width = 7, height = 5.5, units = "in")

ggmap(map) +   
  geom_point(data = sites, aes(x = Longitude, y = Latitude, color = Dive.site),
             size = 3, shape = 19) +
  xlim(c(130.3, 130.85)) +
  ylim(c(-0.70, -0.33)) +
  labs(x = "Longitude", y = "Latitude")

dev.off()

hist(sites$wave_exp_norm2, breaks = 20)

# Site and village data for ggplot
site2 <- sites %>%
  dplyr::select(Site.name, Latitude, Longitude) %>%
  mutate(type = "Site")

village2 <- villages %>%
  rename(Site.name = village) %>%
  dplyr::select(Site.name, Latitude, Longitude) %>%
  mutate(type = "Village")

points.gg <- rbind(site2, village2) %>%
  arrange(desc(type))

# blue background
bbox <- st_as_sf(st_sfc(st_polygon(list(matrix(
  c(130, -0.75,
    131, -0.75,
    131, -0.30,
    130, -0.30,
    130, -0.75), 
  ncol = 2, byrow = TRUE))), crs = 4326))

# Whole MPA
r4_mpa_comb <- st_union(r4_mpa)

# Plot 
jpeg(filename = "../Graphs/map_sites_zones2.jpeg", 
     res = 600, width = 7.5, height = 5, units = "in")

ggplot() +
  geom_sf(data = bbox, fill = "lightblue1", color = NA) +  # Background rectangle
  geom_sf(data = r4_land, fill = "grey60", color = NA) +
  geom_sf(data = r4_mpa, aes(fill = Use_Zone)) +
  #geom_sf(data = r4_mpa_comb, fill = NA, linewidth = 1) +
  scale_fill_manual(name = "MPA zoning", values = c("lightblue1", "#d3d953", "#f59178")) +
  geom_point(data = points.gg, aes(x = Longitude, y = Latitude, shape = type, color = type), size = 2) +
  scale_color_manual(name = "", values = c("grey20", "red3")) +
  scale_shape_manual(name = "", values = c(19, 15)) +
  xlim(c(130.3, 130.85)) +
  ylim(c(-0.70, -0.33)) +
  theme_test() +
  ggsn::scalebar(x.min = 130.30, x.max = 130.40,
                 y.min = -0.69,  y.max = -0.45, transform = TRUE,
                 dist = 5, dist_unit = "km", model = 'WGS84',
                 st.size = 2.5, st.color = "black")

dev.off()

# Plot Indonesia
jpeg(filename = "../Graphs/map_Indo.jpeg", 
     res = 600, width = 2.5, height = 1, units = "in")

ggplot() +
  geom_sf(data = r4_land, fill = "grey40", color = NA) +
  theme_test() +
  theme(axis.title = element_blank(),  # Remove axis titles
        axis.text = element_blank(),   # Remove axis text (tick labels)
        axis.ticks = element_blank())   # Remove axis ticks

dev.off()



# STEPWISE MODEL SELECTION (dbRDA) -------------------------------------------
data_rda <- data_nona

# Hellinger-transformed matrix
fish_comm_hell <- decostand(fish_comm, method = "hellinger")

# Same order as data
fish_comm_hell <- fish_comm_hell[match(data_nona$transect, rownames(fish_comm_hell)), ]

fish_rda_null <- dbrda(fish_comm_hell ~ 1, data = data_rda, distance = "bray")

fish_rda_all <- dbrda(fish_comm_hell ~ coral_SR + fractal_D + cover_deadcoral  +
                        fishing_pressure + zone_bin + wave_exp, 
                      data = data_rda, distance = "bray")

# Forward selection of spatial variables
fish_rda_sel <- ordistep(fish_rda_all, scope = formula(fish_rda_null), 
                         Pin = 0.05, Pout = 0.25, direction = "backward", permutations = 999)

# retained explanatory variables
fish_rda_sel$terms[[3]] 
  # coral_SR + fractal_D + cover_deadcoral + fishing_pressure + zone_bin + wave_exp


# dbRDA ----------------------------------------------------------------------
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
# RDA of fish communities at each transect and how they are affected by certain variables
data_rda <- data_nona

# Hellinger-transformed matrix
fish_comm_hell <- decostand(fish_comm, method = "hellinger")

# Same order as data
fish_comm_hell <- fish_comm_hell[match(data_nona$transect, rownames(fish_comm_hell)), ]

fish_rda <- dbrda(fish_comm_hell ~ coral_SR + fractal_D + cover_deadcoral  +
                    fishing_pressure + zone_bin + wave_exp, 
                  data = data_rda, distance = "bray")


rda_summ <- summary(fish_rda)
anova.cca(fish_rda, by = "term")  
vif.cca(fish_rda)

## Plot
#Scaling 2 shows the effects of explanatory variables
ordiplot(fish_rda, scaling = 2, type = "text")

# Extract variance explained (overall)
var_exp <- round(rda_summ$cont$importance[2, 1:2] * 100, 2)

# Extract variance explained (by constraining factors)
var_exp_c <- round(rda_summ$concont$importance[2, 1:2] * 100, 2)

# Extract scores for sites and variables
# Scaling 1 shows similarities between objects in the response matrix (fish comm matrix).
site_scores <- as.data.frame(vegan::scores(fish_rda, display="sites", choices=c(1,2), scaling=1))
site_scores$zone <- data_rda$zone  # Add zone information
species_scores   <- as.data.frame(vegan::scores(fish_rda, display="species", choices=c(1,2), scaling=1))
variable_scores  <- as.data.frame(vegan::scores(fish_rda, display="bp", choices=c(1, 2), scaling=1))
variable_scores$variable <- rownames(vegan::scores(fish_rda, display = "bp"))

# Remove categorical var "zone"
variable_scores <- variable_scores %>%
  filter(variable != "zoneUZ")

# Rename
variable_scores$variable <- c("Coral SR", "Structural complexity", "Dead coral cover", "Fishing pressure", "No-take", "Wave exposure")

# Create ggplot
scale_arr <- 2

plot_rda1 <- ggplot(site_scores, aes(x = dbRDA1, y = dbRDA2, color = zone)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "grey20") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "grey20") +
  geom_point(size = 3, alpha = 0.7) +
  geom_segment(data = variable_scores, aes(x = 0, y = 0, xend = dbRDA1 * scale_arr, yend = dbRDA2 * scale_arr),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1, color = "black") +
  geom_text_repel(data = variable_scores, aes(x = dbRDA1 * scale_arr, y = dbRDA2 * scale_arr, label = variable),
                  nudge_x = c(0.02, 0, 0, -0.05, 0.01, 0),
                  nudge_y = c(0, 0.02, -0.02, 0, -0.01, -0.02), 
                  segment.color = NA, color = "gray20", size = 4, fontface = "bold", box.padding = 0.2, max.overlaps = Inf)+
  scale_colour_manual(values = c("UZ" = "#f59178", "NTZ" = "yellowgreen"), 
                      labels = c("NO-TAKE", "USE")) +
  #xlim(c(-0.17, 0.17)) +
  #ylim(c(-0.14, 0.14)) +
  theme_bw() +
  labs(x = paste0("dbRDA1: ", var_exp[1], "% of total variation (", var_exp_c[1], "% of fitted variation)"),
       y = paste0("dbRDA2: ", var_exp[2], "% of total variation (", var_exp_c[2], "% of fitted variation)"),
       color = "Zone"); plot_rda1

ggsave(filename = paste0('../Graphs/RDA_fishcomm.jpeg'), plot_rda1,
       dpi = 300, width = 8, height = 6, units = "in")  


# PERMDISP ----------------------------------------------------------------

fish_comm_dist <- vegdist(fish_comm_hell, method = "bray")

## MPA zoning
pdisp_zone <- betadisper(fish_comm_dist, group = data_nona$zone_bin, type = "centroid")

plot(pdisp_zone)
anova(pdisp_zone) # no diff

## Dive site
pdisp_dive <- betadisper(fish_comm_dist, group = data_nona$Dive.site, type = "centroid")

plot(pdisp_dive)
anova(pdisp_dive) # no diff



# SEM ---------------------------------------------------------------------
# https://kevintshoemaker.github.io/NRES-746/SEM.RMarkdown.html
# Structural equation model for drivers of fish diversity
# Five steps in SEM: model specification, identification, parameter estimation, model evaluation, and model modification
# Keep MPA zoning as binary as per: https://jslefche.github.io/sem_book/categorical-variables.html#introduction-to-exogenous-categorical-variables

data_psem <- data_nona %>% 
  select(fish_SR, fish_FDis, fish_RaoQ, 
         coral_SR, fractal_D, cover_deadcoral,
         fishing_pressure, zone_bin, wave_exp, Site.name, sc_imputed)

# Sensitivity analyses:
#data_psem$fish_FDis <- data_psem$fish_RaoQ
#data_psem <- data_psem %>% filter(sc_imputed == 0)

model_psem <- psem(
  lmer(fish_SR ~ coral_SR + fractal_D + cover_deadcoral + fishing_pressure + zone_bin + fishing_pressure:zone_bin + 
         (1 | Site.name), data = data_psem),
  
  lmer(fish_FDis ~ coral_SR + fractal_D + cover_deadcoral + fishing_pressure + zone_bin + 
         fish_SR + # Suggested by D-Sep
         (1 | Site.name), data = data_psem),

  lmer(coral_SR ~ fishing_pressure + zone_bin + wave_exp + 
         (1 | Site.name), data = data_psem),
  
  lmer(fractal_D ~ fishing_pressure + zone_bin + wave_exp + 
         coral_SR + # Suggested by D-Sep
         (1 | Site.name), data = data_psem),
  
  lmer(cover_deadcoral ~ fishing_pressure + zone_bin + wave_exp + 
         coral_SR + fractal_D + # Suggested by D-Sep
         (1 | Site.name),  data = data_psem)
)


summary(model_psem) 
AIC(model_psem) 
r2_psem <- rsquared(model_psem, method = NULL)
sem_coeff <- coefs(model_psem)
#write.csv(sem_coeff, "./results/sem_coefficients.csv", row.names = FALSE)

# Tests of directed separation
(dsep_psem <- piecewiseSEM::dSep(model_psem))
# Suggestions added in model above

# Plot (too basic)
plot(model_psem)

## Plot coefficients
coef_df <- coefs(model_psem)[, 1:(ncol(coefs(model_psem))-1)] # Last column is significance code ( *, **, *** )

# Covariances (report in text)
covar_df <- coef_df %>%
  filter((Predictor == "coral_SR" & Response == "cover_deadcoral") |
           (Predictor == "coral_SR" & Response == "fractal_D") |
           (Predictor == "fish_SR" & Response == "fish_FDis") |
           (Predictor == "fractal_D" & Response == "cover_deadcoral"))

#write.csv(covar_df, "./results/sem_coeff_covars.csv", row.names = FALSE)

# Direct effects 
coef_df <- coef_df %>% 
  filter(!(Predictor == "coral_SR" & Response == "cover_deadcoral")) %>% 
  filter(!(Predictor == "coral_SR" & Response == "fractal_D")) %>% 
  filter(!(Predictor == "fish_SR" & Response == "fish_FDis")) %>% 
  filter(!(Predictor == "fractal_D" & Response == "cover_deadcoral"))

#write.csv(coef_df, "./results/sem_coeff_direct.csv", row.names = FALSE)

# Indirect effects
dir_eff <- coef_df %>%
  select(Predictor, Response, Estimate, Std.Error, P.Value) %>% 
  filter(P.Value < 0.05) # keep only significant paths

indirect_paths <- dir_eff %>%
  rename(Mediator = Response) %>%
  inner_join(dir_eff, by = c("Mediator" = "Predictor")) %>% # (ignore warning)
  transmute(Predictor = Predictor,
            Mediator = Mediator,
            Response = Response,
            Indirect_Effect = Estimate.x * Estimate.y) %>%
  arrange(Response)

#write.csv(indirect_paths, "./results/sem_coeff_indirect.csv", row.names = FALSE)

# Uncertainty and Significance
coef_df <- coef_df %>%
  mutate(lower_95 = Estimate - 1.96 * Std.Error,
         upper_95 = Estimate + 1.96 * Std.Error,
         sig_color = case_when(P.Value > 0.05 ~ "nonsig", 
                               P.Value < 0.05 & Estimate < 0 ~ "neg_sig", # Negative and significant
                               TRUE ~ "sig")) # Positive and significant

# Clean names
coef_df <- coef_df %>%
  mutate(
    Predictor = case_when(
      Predictor == "coral_SR" ~ "Coral SR",
      Predictor == "fractal_D" ~ "SC",
      Predictor == "cover_deadcoral" ~ "DCC",
      Predictor == "fishing_pressure" ~ "FP",
      Predictor == "fishing_pressure_UZ" ~ "Fishing pressure (UZ)",
      Predictor == "fishing_pressure_NTZ" ~ "Fishing pressure (NTZ)",
      Predictor == "fishing_pressure:zone_bin" ~ "FP:NTZ",
      Predictor == "zone_bin" ~ "NTZ",
      Predictor == "wave_exp" ~ "Wave Exp",
      Predictor == "fish_SR" ~ "Fish SR",
      TRUE ~ Predictor  # Keep the original value if it doesn't match any condition
    ),
    Response = case_when(
      Response == "coral_SR" ~ "Coral SR",
      Response == "fractal_D" ~ "Structural complexity",
      Response == "cover_deadcoral" ~ "Dead coral cover",
      Response == "fish_SR" ~ "Fish SR",
      Response == "fish_FDis" ~ "Fish FDis",
      TRUE ~ Response  # Keep the original value if it doesn't match any condition
    )
  )

# Re-order predictor from seascape to habitat
coef_df$Predictor <- factor(coef_df$Predictor, 
                            levels = c("Wave Exp", "NTZ", "FP", 
                                       "FP:NTZ",
                                       "DCC", "SC", "Coral SR",
                                       "Fish SR"))

# Re-order response factors from habitat to fish
coef_df$Response <- factor(coef_df$Response, levels = c("Coral SR", "Structural complexity", "Dead coral cover", "Fish SR", "Fish FDis", "Fish PC1"))

# Facet labels
r2_labels <- c(
  "Coral SR" = paste0("Coral SR (R² = ", round(r2_psem$Conditional[r2_psem$Response == "coral_SR"], 2), ")"),
  "Structural complexity" = paste0("Structural complexity (R² = ", round(r2_psem$Conditional[r2_psem$Response == "fractal_D"], 2), ")"),
  "Dead coral cover" = paste0("Dead coral cover (R² = ", round(r2_psem$Conditional[r2_psem$Response == "cover_deadcoral"], 2), ")"),
  "Fish SR" = paste0("Fish SR (R² = ", round(r2_psem$Conditional[r2_psem$Response == "fish_SR"], 2), ")"),
  "Fish FDis" = paste0("Fish FDis (R² = ", round(r2_psem$Conditional[r2_psem$Response == "fish_FDis"], 2), ")")
)

# Plot coefficients 
plot_coeff_psem <- ggplot(coef_df, aes(x = Estimate, y = Predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey30") +
  geom_point(aes(color = sig_color), size = 4, shape = 19) +
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95, color = sig_color), linewidth = 0.8, height = 0) +
  facet_wrap(~Response, scales = "free_y", ncol = 1, labeller = labeller(Response = r2_labels)) +
  scale_color_manual(values = c("nonsig" = "grey70", "sig" = "black", "neg_sig" = "red2")) +
  xlim(c(-1.13, 1.13)) +
  labs(x = "Estimate", y = NULL) +
  theme_test() +
  theme(legend.position = "none"); plot_coeff_psem

ggsave(filename = '../Graphs/PSEM_coeff.jpeg', plot_coeff_psem,
       dpi = 300, width = 6.75, height = 18, units = "cm")




###-######################-###
###### ***APPENDIX*** ########
###-######################-###

# REEF RUGOSITY ~ FRACTAL D -----------------------------------------------

lm_sc <- lm(rugosity ~ fractal.d, data = data_reef)
summary(lm_sc)
anova(lm_sc)

jpeg(filename = "../Graphs/SC_rugosity-FractalD.jpeg", 
     res = 300, width = 6, height = 5, units = "in")

ggplot(data_reef, aes(x = fractal.d, y = rugosity)) +
  geom_point() + 
  geom_text_repel(aes(label = transect)) +
  stat_smooth(method = "lm", linewidth = 1, fill = "grey80") +
  xlab("Reef fractal dimension") +
  ylab("Reef rugosity") +
  theme_bw()

dev.off()


# FISH FRic ~ SR ---------------------------------------------------

lm_sf <- lm(fish.FRic ~ fish.SR, data = data_reef)
summary(lm_sf)
anova(lm_sf)

jpeg(filename = "../Graphs/fish_SR-FRic2.jpeg", 
     res = 300, width = 6, height = 5, units = "in")

ggplot(data_reef, aes(x = fish.SR, y = fish.FRic)) +
  geom_point() + 
  geom_text_repel(aes(label = transect)) +
  geom_smooth(method = "lm") +
  xlab("Fish species richness") +
  ylab("Fish functional richness (FRic)") +
  theme_bw()

dev.off()


# FISH RaoQ ~ FDis ---------------------------------------------------

lm_rao_fdis <- lm(fish.RaoQ ~ fish.FDis, data = data_reef)
summary(lm_rao_fdis)
anova(lm_rao_fdis)

jpeg(filename = "../Graphs/fish_RaoQ-FDis.jpeg", 
     res = 300, width = 6, height = 5, units = "in")

ggplot(data_reef, aes(x = fish.FDis, y = fish.RaoQ)) +
  geom_point() + 
  geom_text_repel(aes(label = transect)) +
  geom_smooth(method = "lm") +
  xlab("Fish functional dispersion (FDis)") +
  ylab("Fish Rao's Q") +
  theme_bw()

dev.off()



# FISH FDis ~ DCC ---------------------------------------------------

lm_fdis_dcc <- lm(fish_FDis ~ cover_deadcoral, data = data_nona)
summary(lm_fdis_dcc)
anova(lm_fdis_dcc)
AIC(lm_fdis_dcc)

lm_fdis_dcc2 <- lm(fish_FDis ~ cover_deadcoral + I(cover_deadcoral^2), data = data_nona)
summary(lm_fdis_dcc)
AIC(lm_fdis_dcc2)

anova(lm_fdis_dcc, lm_fdis_dcc2)


jpeg(filename = "../Graphs/fish_FDis-DCC_quad.jpeg", 
     res = 300, width = 6, height = 5, units = "in")

ggplot(data_nona, aes(x = cover_deadcoral, y = fish_FDis)) +
  geom_point() + 
  #geom_text_repel(aes(label = transect)) +
  #geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE, color = "firebrick2") +  # Quadratic
  xlab("Dead coral cover (scaled)") +
  ylab("Fish functional dispersion (scaled)") +
  theme_bw()

dev.off()

# Range in relative dead coral cover
range(data_reef$cover_dead.coral)

jpeg(filename = "../Graphs/fish_FDis-DCC_quad_notscaled.jpeg", 
     res = 300, width = 6, height = 5, units = "in")

ggplot(data_reef, aes(x = cover_dead.coral, y = fish.FDis)) +
  geom_point() + 
  #geom_text_repel(aes(label = transect)) +
  #geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE, color = "firebrick2") +  # Quadratic
  xlab("Dead coral cover (proportion)") +
  ylab("Fish functional dispersion (FDis)") +
  theme_bw()

dev.off()

# FISH FDis ~ coral SR ---------------------------------------------------

lm_fdis_csr <- lm(fish_FDis ~ coral_SR, data = data_nona)
summary(lm_fdis_csr)
anova(lm_fdis_csr)
AIC(lm_fdis_csr)

lm_fdis_csr2 <- lm(fish_FDis ~ coral_SR + I(coral_SR^2), data = data_nona)
summary(lm_fdis_csr)
AIC(lm_fdis_csr2)

anova(lm_fdis_csr, lm_fdis_csr2)


jpeg(filename = "../Graphs/fish_FDis-CoralSR_quad.jpeg", 
     res = 300, width = 6, height = 5, units = "in")

ggplot(data_nona, aes(x = coral_SR, y = fish_FDis)) +
  geom_point() + 
  #geom_text_repel(aes(label = transect)) +
  #geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE, color = "firebrick2") +  # Quadratic
  xlab("Coral SR (scaled)") +
  ylab("Fish functional dispersion (scaled)") +
  theme_bw()

dev.off()

# FISH FDis ~ fractal D ---------------------------------------------------

lm_fdis__frd <- lm(fish_FDis ~ fractal_D, data = data_nona)
summary(lm_fdis__frd)
anova(lm_fdis__frd)
AIC(lm_fdis__frd)

lm_fdis__frd2 <- lm(fish_FDis ~ fractal_D + I(fractal_D^2), data = data_nona)
summary(lm_fdis__frd)
AIC(lm_fdis__frd2)

anova(lm_fdis__frd, lm_fdis__frd2)


jpeg(filename = "../Graphs/fish_FDis-fractalD_quad.jpeg", 
     res = 300, width = 6, height = 5, units = "in")

ggplot(data_nona, aes(x = fractal_D, y = fish_FDis)) +
  geom_point() + 
  #geom_text_repel(aes(label = transect)) +
  #geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE, color = "firebrick2") +  # Quadratic
  xlab("Fractal dimension (scaled)") +
  ylab("Fish functional dispersion (scaled)") +
  theme_bw()

dev.off()

# FISH NMDS ---------------------------------------------------------------
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/nmds/
# 
#fish_comm <- fish_comm[rownames(fish_comm) %in% data_reef$transect,]

# Transform comm data and calculate dist
fish_comm_tr <- fish_comm^0.25 

#fish_comm_dist <- vegdist(fish_comm_tr, method = "bray")

# NMDS 
fish.nmds <- metaMDS(comm = fish_comm_tr, distance = "bray",
                     autotransform = FALSE, engine = "monoMDS", 
                     k = 3,
                     weakties = TRUE,
                     model = "global",
                     maxit = 999,
                     try = 40,
                     trymax = 500,
                     wascores = TRUE)

# Diagnostics
(stress <- round(fish.nmds$stress, 3))

(gof <- goodness(object = fish.nmds)) # proportion of the overall variance that is explained by each sample unit
barplot(gof) # Should be even-ish

plot(fish.nmds, display = "sites", type = "none")
points(fish.nmds, display = "sites", cex = 2*gof/mean(gof)) # Sample units shown with larger circles account for more of the variation

stressplot(object = fish.nmds, lwd = 5) # Evaluate fit (Shepard plot)

## Organize for ggplot
nmds_scores <- as.data.frame(fish.nmds$points)  
nmds_scores$transect <- rownames(nmds_scores)  
nmds_scores$zone <- transects$zone[match(nmds_scores$transect, transects$Transect.code)]
nmds_scores$site <- strtrim(nmds_scores$transect, 3)
nmds_scores$fp.sasi <- sites$fp.sasi[match(nmds_scores$site, sites$Site.code)] 

# Fit species vectors to the NMDS ordination
spp_cor <- cor(fish_comm_tr,
               fish.nmds$points,
               use = "complete.obs",
               method = "pearson")

spp_cor <- spp_cor[,1:2] # Keep first two axes

spp_cor_f <- spp_cor[apply(abs(spp_cor) > 0.5, 1, any), ] # Keep species with Pearson r > 0.5

# Convert to dataframe for ggplot
spp_cor_df <- as.data.frame(spp_cor_f)
spp_cor_df$species <- gsub("_", " ", rownames(spp_cor_df))
spp_cor_df$species_italic <- paste0("italic('", spp_cor_df$species, "')")


# Plot NMDS
plot_nmds1 <- ggplot() + 
  geom_point(data = nmds_scores, aes(x = MDS1, y = MDS2, shape = zone, colour = fp.sasi), size = 4) + 
  #geom_text_repel(data = nmds_scores, aes(x = MDS1, y = MDS2, label = transect), 
  #                size = 1.5, vjust = 0.5, hjust = 0.5) +
  scale_shape_manual(name = "MPA zoning", values = c("UZ" = 16, "NTZ" = 17)) +  
  scale_color_gradientn(name = "Fishing\npressure", colors = c("green", "greenyellow", "yellow", "orange", "firebrick2")) +
  annotate("text", x = max(nmds_scores$MDS1), y = max(nmds_scores$MDS2), 
           label = paste("Stress =", stress), hjust = 1, vjust = 1, size = 5) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "left",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()); plot_nmds1

# Plot species scores
scale_spp <- 0.7
plot_nmds2 <- ggplot() + 
  geom_point(data = nmds_scores, aes(x = MDS1, y = MDS2), alpha = 0) + 
  geom_segment(data = spp_cor_df, 
               aes(x = 0, y = 0, xend = MDS1 * scale_spp, yend = MDS2 * scale_spp), 
               #arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +
  geom_text_repel(data = spp_cor_df, position = position_nudge_center(0.3, 0.1, 0, 0, direction = "radial"),
                  aes(x = MDS1 * scale_spp, y = MDS2 * scale_spp, label = species_italic), 
                  parse = TRUE, size = 3, max.overlaps = Inf, min.segment.length = Inf) +
  coord_equal() +
  theme_void(); plot_nmds2

plot_nmds_comb <- plot_nmds1 + plot_nmds2

# Save 
ggsave("../Graphs/fish_NMDS_zone_fp_4root.jpeg", plot_nmds_comb, width = 10, height = 5, dpi = 300)

# Check some sites and secies to see if they match species scores
sort(fish_comm["SGK01", ], decreasing = TRUE)[1:10]
rownames(fish_comm)[order(fish_comm[, "Chelmon_rostratus"], decreasing = TRUE)][1:10]

## Check similarity of transects within sites 
# Draw lines between pairs of transects for each site, color by fishing pressure

jpeg(filename = "../Graphs/fish_NMDS_pairs_FPsasi2.jpeg", res = 300, width = 6, height = 5, units = "in")

ggplot() + 
  geom_point(data = nmds_scores, aes(x=MDS1,y=MDS2,colour=fp.sasi),size=4) + 
  geom_line(data = nmds_scores, aes(x=MDS1, y=MDS2, group = site)) +
  geom_text(data=nmds_scores,aes(x=MDS1,y=MDS2,label=transect),size=2,vjust=0,hjust=0) +  
  scale_color_gradientn(name = "Fishing\npressure", colors = c("green", "greenyellow", "yellow", "orange", "firebrick2")) +
  coord_equal() +
  theme_bw()

dev.off()



# COMPARE MATRIX TRANSFORMATIONS ------------------------------------------

# Fourth root is better for PCoA and NMDS (works well with Bray-Curtis dist)
dist1 <- vegdist(fish_comm^0.25, method = "bray")

# Hellinger is better for PCA and RDA (works well with Euclidean dist)
dist2 <- vegdist(decostand(fish_comm, method = "hellinger"), method = "euclidean")

mantel(dist1, dist2, method = "pearson", permutations = 999)


# PLOT FP:ZONE INTERACTION ------------------------------------------------

p1 <- ggplot(data_reef, aes(x = fp.sasi, y = fish.SR)) +
  geom_point(aes(color = zone)) + 
  stat_smooth(aes(color = zone), fill = "grey90", method = "lm", linewidth = 1, se = TRUE) +
  scale_colour_manual(values = c("UZ" = "#f59178", "NTZ" = "yellowgreen"),
                      labels = c("NO-TAKE", "USE")) +
  xlab(NULL) +
  ylab("Fish species richness") +
  labs(color = "MPA zoning:") +
  theme_bw() +
  theme(legend.position = "top")

p2 <- ggplot(data_reef, aes(x = fp.sasi, y = fish.FDis)) +
  geom_point(aes(color = zone)) + 
  stat_smooth(aes(color = zone), fill = "grey90", method = "lm", linewidth = 1, se = TRUE) +
  scale_colour_manual(values = c("UZ" = "#f59178", "NTZ" = "yellowgreen"),
                      labels = c("NO-TAKE", "USE")) +
  xlab("Reef fishing pressure index") +
  ylab("Fish functional dispersion") +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork with annotations
inter_plot <- (p1 / p2) +
  plot_annotation(tag_levels = "a")

# Save the combined plot
ggsave(filename = "../Graphs/Supp_interaction_plots.jpeg", 
       plot = inter_plot, dpi = 300, width = 4, height = 6, units = "in")

