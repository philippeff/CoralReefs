# Reef Fish Taxonomic and Functional Diversity Respond Differently to Habitat Degradation and Protection

This repository contains data and R scripts for analyses presented in a chapter of my PhD thesis.

---

## Abstract

Coral reef ecosystems are impacted by human-induced pressures and are shaped by complex interactions between habitat characteristics and environmental stressors, yet the relative influence of these factors on reef fish diversity remains unclear. Here, we examine the effects of coral species richness (SR), structural complexity, fishing pressure, marine protection and wave exposure on reef fish taxonomic and functional diversity in a multi-use marine protected area (MPA) in Raja Ampat, Indonesia. We found that fish SR was strongly associated with both coral SR and structural complexity, while fish functional diversity showed no significant relationship with these habitat variables. Instead, reef fish functional dispersion, which is a measure of a community’s distribution in functional trait space, was (unexpectedly) higher in reefs with higher relative dead coral cover and decreased in reefs that were closer to settlements – and so putatively experiencing higher fishing pressure. Moreover, higher fishing pressure reduced structural complexity and increased dead coral cover, thereby indirectly decreasing fish SR. We found that marine protection mitigated some of these impacts and supported more structurally complex reefs. Our results highlight the need to incorporate measures of reef accessibility and fine-scale habitat variability when designing conservation strategies for complex and biodiverse coral reef ecosystems. We emphasize the need to understand the interactions between habitat features and environmental stressors for improving conservation management, particularly as human activities continue to impact coral reefs worldwide.

---

## R Scripts

- `01_navig_dist.R` — Calculates navigable distances from reefs to villages (for index of fishing pressure).
- `02_variables.R` — Processes, calculates and combines fishing pressure index, structural complexity, fish functional diversity, and relative cover.
- `03_wave_exposure_site.R` — Estimates wave exposure per site using the Global Wind Atlas.
- `04_analysis.R` — Main analyses and plots.

## Data

- `data/data_coral_final.csv` — Coral identification data.
- `data/data_coral_reefcloud_2024-05-26.csv` — Coral annotation data from ReefCloud.
- `data/data_fish_final.csv` — Reef fish survey data (species, abundance and richness).
- `data/data_fishtraits_PFF.csv` — Functional trait data for fish species.
- `data/data_reef_final.csv` — Combined transect-level dataset.
- `data/data_sites_final.csv` — Site-level metadata.
- `data/data_transects_final.csv` — Transect-level metadata and habitat data.
- `data/data_villages.csv` — Village locations for fishing pressure index.
- `data/meta_fish_species.csv` — Fish species metadata (e.g. names, families and total abundance).
- `data/maps/Selat_Dampier.shp` — Shapefile for the Selat Dampier marine protected area.

---

