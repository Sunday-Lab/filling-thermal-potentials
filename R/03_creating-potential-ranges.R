## script to create potential range shapefiles for each species based on their thermal tolerance limits
library(tidyverse)
library(sp)
library(sf)
library(raster)
library(ncdf4)
library(rnaturalearth)
library(smoothr)
library(lwgeom)
library(rmapshaper)
library(RColorBrewer)
select <- dplyr::select



####################################################################################
#####                       READ IN TEMPERATURE RASTERS                       ######
####################################################################################
raster_terr_low <- raster("./data-processed/raster_terr_low.nc")
raster_terr_high <- raster("./data-processed/raster_terr_high.nc")
raster_marine_low <- raster("./data-processed/raster_marine_low.nc")
raster_marine_high <- raster("./data-processed/raster_marine_high.nc")
raster_intertidal_low <- raster("./data-processed/raster_intertidal_low.nc")
raster_intertidal_high <- raster("./data-processed/raster_intertidal_high.nc")

####################################################################################
#####                   RASTERIZE REALIZED RANGES                             ######
####################################################################################
## read in realized ranges:
realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp") %>%
  mutate(range_id = paste(species, source, sep = "_"))
split_realized_ranges <- st_read("data-processed/realized-ranges_split.shp")


## constrain realized range rasters by habitat allowed to be in the potential range
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1)

## read in realm mask layers:
t_mask <- raster("./data-processed/raster_terr_mask.nc")
m_mask <- raster("./data-processed/raster_marine_mask.nc")
i_mask <- raster("./data-processed/raster_intertidal_mask.nc")

i = 1
while (i < nrow(realized_ranges)+1) {
  ## for ranges that are very complex and take forevvvvver to rasterize, simplify first:
  if (i %in% c(100, 173, 190, 191, 192, 196, 197, 198, 199, 205, 206, 209, 225, 229, 426, 524)) {
    range <- ms_simplify(realized_ranges[i, ], keep_shapes = TRUE)
  }
  else {
    range <- realized_ranges[i, ]
  }
  rr_raster <- rasterize(range, r, background = NA, getCover = TRUE)
  rr_raster[rr_raster > 0] <- 1 ## set all values that overlap realized range in raster to 1
  rr_raster[rr_raster != 1] <- NA ## set all other values to NA
  
  ## constrain by habitat:
  if(range$realm == "Terrestrial" | range$realm == "Freshwater") {
    rr_raster <- mask(rr_raster, t_mask) ## set cells in mask that do not overlap realized range to NA
  }
  else if(range$realm == "Marine") {
    rr_raster <- mask(rr_raster, m_mask) 
  }
  else {
    rr_raster <- mask(rr_raster, i_mask) 
  }
  
  ## add to list of rasters
  if (i == 1) {
    rasterized_rrs <- rr_raster
  }
  else {
    rasterized_rrs <- addLayer(rasterized_rrs, rr_raster, updatevalue = NA)
  }
  
  print(paste("Finished range number:", i))
  i = i + 1
}

names(rasterized_rrs) <- realized_ranges$range_id

##saveRDS(rasterized_rrs, "data-processed/rasterized_rrs.rds")
rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")


####################################################################################
#####                   CREATING POTENTIAL RANGE SHAPEFILES                   ######
####################################################################################
## read in thermal limits:
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = " "))

## read in species traits:
traits <- read.csv("./data-processed/wrangled-traits.csv")

#######################################################
#####   TERRESTRIAL AND FRESHWATER SPECIES:      ######
#######################################################
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]
only_upper <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]

## clump contiguous habitat together:
clumped_temps <- raster_terr_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## SPEICES WITH BOTH THERMAL LIMITS:
####################################=
combined <- filter_by_tolerance_both(both_upper = both_upper, 
                                    both_lower = both_lower, 
                                    raster_high = raster_terr_high, 
                                    raster_low = raster_terr_low)

combined_dormancy <- filter_by_tolerance_with_dormancy(both_upper = both_upper, 
                                                       both_lower = both_lower, 
                                                       raster_high = raster_terr_high, 
                                                       raster_low = raster_terr_low)
plot(combined)
plot(combined_dormancy)

## restrict range to contiguous habitat that begins at the species realized range:
prs_terrestrial <- create_potential_ranges(clumped_temps = clumped_temps, 
                                           realized_ranges = realized_ranges, combined = combined)
prs_terrestrial_dormancy <- create_potential_ranges(clumped_temps = clumped_temps, 
                                           realized_ranges = realized_ranges, combined = combined_dormancy)

##saveRDS(prs_terrestrial_dormancy, "./data-processed/prs_terrestrial_dormancy_6mo.rds")
prs_terrestrial_dormancy <- readRDS("./data-processed/prs_terrestrial_dormancy_6mo.rds")

## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################
# high_filtered <- filter_by_tolerance_upper(only_upper = only_upper, 
#                                            raster_high = raster_terr_high)
# plot(high_filtered)
# low_filtered <- filter_by_tolerance_lower(only_lower = only_lower, 
#                                           raster_low = raster_terr_low)
# plot(low_filtered)
# 
# ## restrict range to contiguous habitat that begins at the species realized range:
# prs_terrestrial_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
#                                                                realized_ranges = realized_ranges,
#                                                                high_filtered = high_filtered,  
#                                                                low_filtered = low_filtered)



#######################################################
#####                MARINE SPECIES:             ######
#######################################################
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Marine")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Marine")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]
only_upper <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]

## clump contiguous habitat together:
clumped_temps <- raster_marine_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## SPECIES WITH BOTH THERMAL LIMITS:
####################################
combined <- filter_by_tolerance_both(both_upper = both_upper,
                                     both_lower = both_lower,
                                     raster_high = raster_marine_high, 
                                     raster_low = raster_marine_low)
plot(combined)

## restrict range to contiguous habitat that begins at the species realized range:
prs_marine <- create_potential_ranges(clumped_temps = clumped_temps, 
                                      realized_ranges = realized_ranges, combined = combined)

##saveRDS(prs_marine, "./data-processed/prs_marine.rds")
prs_marine <- readRDS("./data-processed/prs_marine.rds")

## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################
# high_filtered <- filter_by_tolerance_upper(only_upper = only_upper,
#                                            raster_high = raster_marine_high)
# plot(high_filtered)
# low_filtered <- filter_by_tolerance_lower(only_lower = only_lower,
#                                           raster_low = raster_marine_low)
# plot(low_filtered)
# 
# ## restrict range to contiguous habitat that begins at the species realized range:
# prs_marine_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
#                                                           realized_ranges = realized_ranges, 
#                                                           high_filtered = high_filtered,  
#                                                           low_filtered = low_filtered)


######################################################
#####            INTERTIDAL SPECIES:            ######
######################################################
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Intertidal")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Intertidal")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]
only_upper <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]

## clump contiguous habitat together:
clumped_temps <- raster_intertidal_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps$geometry)

## SPECIES WITH BOTH THERMAL LIMITS:
####################################
combined <- filter_by_tolerance_both(both_upper = both_upper,
                                     both_lower = both_lower,
                                     raster_high = raster_intertidal_high,
                                     raster_low = raster_intertidal_low)
plot(combined)

## restrict range to contiguous habitat that begins at the species realized range:
prs_intertidal <- create_potential_ranges(clumped_temps = clumped_temps, 
                                          realized_ranges = realized_ranges, combined = combined)

##saveRDS(prs_intertidal, "./data-processed/prs_intertidal.rds")
prs_intertidal <- readRDS("./data-processed/prs_intertidal.rds")


## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################
# high_filtered <- filter_by_tolerance_upper(only_upper = only_upper,
#                                            raster_high = raster_intertidal_high)
# plot(high_filtered)
# low_filtered <- filter_by_tolerance_lower(only_lower = only_lower,
#                                           raster_low = raster_intertidal_low)
# plot(low_filtered)
# 
# ## restrict range to contiguous habitat that begins at the species realized range:
# prs_intertidal_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
#                                                               realized_ranges = realized_ranges, 
#                                                               high_filtered = high_filtered,  
#                                                               low_filtered = low_filtered)


########################################################
##      COLLATING POTENTIAL RANGES ACROSS REALMS     ###
########################################################
potential_ranges_both_limits <- stack(prs_terrestrial, prs_marine, prs_intertidal)
potential_ranges_both_limits_dormancy <- stack(prs_terrestrial_dormancy, prs_marine, prs_intertidal)
##saveRDS(potential_ranges_both_limits, "data-processed/potential_ranges.rds")
##saveRDS(potential_ranges_both_limits, "data-processed/potential_ranges_notcutatequator.rds")
##saveRDS(potential_ranges_both_limits_dormancy, "data-processed/potential_ranges_notcutatequator_dormancy.rds")


#potential_ranges_one_limit <- stack(prs_terrestrial_one_limit, prs_marine_one_limit, 
#                                    prs_intertidal_one_limit)
##saveRDS(potential_ranges_one_limit, "data-processed/potential_ranges_one_limit_coldwarm.rds")


## split ranges with both limits into a version with equatorward and poleward ranges
potential_ranges_split 
## combine with potential ranges with only one limit
potantial_ranges_split_all <- stack(potential_ranges_split, potential_ranges_one_limit)






######################################################
##                  FUNCTIONS:                      ##
######################################################

## functions to filter rasterized temperature data (raster_high, raster_low) by thermal tolerance limits (either upper, lower, or both) 
filter_by_tolerance_both <- function(both_upper, 
                                     both_lower, 
                                     raster_high, 
                                     raster_low) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  species = 1
  while (species < nrow(both_upper) + 1) {
    high <- addLayer(high, raster_high[[1]] - both_upper$thermal_limit[species]) 
    low <- addLayer(low, raster_low[[1]] - both_lower$thermal_limit[species]) 
    
    species = species + 1
  }
  
  ## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
  high[high > 0] <- NA
  low[low < 0] <- NA
  ##plot(high)
  ##plot(high)
  
  ## combine to find cells where seasonal high temp is less than CTmax and seasonal low temp is greater than CTmin 
  combined <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  i = 1  
  while (i < nrow(both_upper) + 1) {
    combined <- addLayer(combined, mask(high[[i]], low[[i]]), updatevalue = NA)
    
    i = i + 1
  }
  names(combined) <- paste(both_upper$Genus, both_upper$Species, sep = "_")
  ##plot(combined)
  
  return (combined)
}

filter_by_tolerance_with_dormancy <- function(both_upper, 
                                                   both_lower, 
                                                   raster_high, 
                                                   raster_low) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  species = 1
  high_names <- c()
  low_names <- c()
  tracker <- c()
  while (species < nrow(both_upper) + 1) {
    species_traits <- traits[which(traits$genus_species == both_upper$genus_species[species]),]
    
    if (species_traits$cold_season_dormancy_ == "Yes" & species_traits$hot_season_dormancy_ == "Yes") {
      high <- addLayer(high, raster_high[[1:2]] - both_upper$thermal_limit[species]) 
      low <- addLayer(low, raster_low[[1:2]] - both_lower$thermal_limit[species]) 
      high_names <- append(high_names, paste(paste(species_traits$Genus, species_traits$Species, sep = "_"),
                                             c(0,6), "months_both_dormancy", sep = "_"))
      low_names <- append(low_names, paste(paste(species_traits$Genus, species_traits$Species, sep = "_"),
                                           c(0,6), "months_both_dormancy", sep = "_"))
      tracker <- append(tracker, rep("both", 2))
    }
    else if (species_traits$cold_season_dormancy_ == "Yes") {
      high <- addLayer(high, raster_high[[1]] - both_upper$thermal_limit[species]) 
      high <- addLayer(high, raster_high[[1]] - both_upper$thermal_limit[species]) 
      low <- addLayer(low, raster_low[[1:2]] - both_lower$thermal_limit[species]) 
      high_names <- append(high_names, rep(paste(species_traits$Genus, species_traits$Species, 
                                             "no_hot_dormancy", sep = "_"), 2))
      low_names <- append(low_names, paste(paste(species_traits$Genus, species_traits$Species, sep = "_"),
                                            c(0,6), "months_cold_dormancy", sep = "_"))
      tracker <- append(tracker, rep("cold", 2))
    }
    else if (species_traits$hot_season_dormancy_ == "Yes") {
      high <- addLayer(high, raster_high[[1:2]] - both_upper$thermal_limit[species]) 
      low <- addLayer(low, raster_low[[1]] - both_lower$thermal_limit[species])
      low <- addLayer(low, raster_low[[1]] - both_lower$thermal_limit[species])
      high_names <- append(high_names, paste(paste(species_traits$Genus, species_traits$Species, sep = "_"),
                                            c(0,6), "months_hot_dormancy", sep = "_"))
      low_names <- append(low_names, rep(paste(species_traits$Genus, species_traits$Species, "no_cold_dormancy",
                                           sep = "_"),2))
      tracker <- append(tracker, rep("hot", 2))
    }
    else {
      high <- addLayer(high, raster_high[[1]] - both_upper$thermal_limit[species]) 
      low <- addLayer(low, raster_low[[1]] - both_lower$thermal_limit[species]) 
      high_names <- append(high_names, paste(species_traits$Genus, species_traits$Species, 
                                             "no_dormancy", sep = "_"))
      low_names <- append(low_names, paste(species_traits$Genus, species_traits$Species, "no_dormancy",
                                           sep = "_"))
      tracker <- append(tracker, "neither")
    }
    
    species = species + 1
  }
  
  ## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
  high[high > 0] <- NA
  low[low < 0] <- NA
  names(high) <- high_names
  names(low) <- low_names
  ##plot(high)
  ##plot(high)
  
  ## combine to find cells where seasonal high temp is less than CTmax and seasonal low temp is greater than CTmin 
  combined <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  c_names <- c()
  i = 1  
  while (i < nlayers(high) + 1) {
    combined <- addLayer(combined, mask(high[[i]], low[[i]]), updatevalue = NA)
    c_names <- append(c_names, ifelse(tracker[i] == "both", high_names[i],
                                      ifelse(tracker[i] == "cold", low_names[i], 
                                             ifelse(tracker[i] == "hot", high_names[i],
                                                    ifelse(tracker[i] == "neither", low_names[i],
                                                           "oop")))))
    i = i + 1
  }
  names(combined) <- c_names
  ##plot(combined)
  
  return (combined)
}

filter_by_tolerance_upper <- function(only_upper, 
                                      raster_high) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  ## for upper limits, only filter use seasonal high temps 
  species = 1
  while (species < nrow(only_upper) + 1) {
    high <- addLayer(high, raster_high[[1]] - only_upper$thermal_limit[species]) 
    
    species = species + 1
  }
  ## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
  high[high > 0] <- NA
  names(high) <- c(paste(only_upper$Genus, only_upper$Species, sep = "_"))
  ##plot(high)
  
  return(high)
}

filter_by_tolerance_lower <- function(only_lower, 
                                      raster_low) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## for lower limits, only filter use seasonal low temps 
  species = 1
  while (species < nrow(only_lower) + 1) {
    low <- addLayer(low, raster_low[[1]] - only_lower$thermal_limit[species]) 
    
    species = species + 1
  }
  
  ## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
  low[low < 0] <- NA
  if (!is.na(minValue(low)[1])) {
    names(low) <- c(paste(only_lower$Genus, only_lower$Species, sep = "_"))
  }
  ##plot(low)
  
  return(low)
}


## returns a potential range raster containing a layer for all species in combined that represents their potential range containing only contiguous habitat in clumps that their realized range touches 
create_potential_ranges <- function (clumped_temps, 
                                     combined, 
                                     realized_ranges) {
  
  pr_restricted_all <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  names <- c()
  i = 1
  while (i < length(names(combined)) + 1) {
    ## get unrestricted potential range 
    potential_raster <- combined[[i]] 
    potential_range <- potential_raster %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf()
    
    ## get realized range
    species <- names(potential_raster) %>%
      str_replace_all("_", " ") 
    
    if(str_detect(species, "dormancy") ==  TRUE) {
      split <- str_split_fixed(species, " ", n = 3)
      
      species <- paste(split[1,1], split[1,2], sep = " ")
    }
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
    
    ## if species has multiple realized ranges (IUCN and GBIF), loop through them:
    num = 1
    while (num < nrow(realized_range)+1) {
      ## overlay realized range with potential range
      ## restrict potential range to clumps that overlap with the realized range
      intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
      rr <- filter(clumped_temps, intersects == TRUE)
      
      intersects_potential <- st_intersects(rr, potential_range, sparse = FALSE)[,]
      
      if (nrow(rr) == 1) {
        pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                  arr.ind=TRUE),]
        pr_rasterized <- rasterize(pr_multipolygons, r, getCover = TRUE)
        pr_restricted <- mask(pr_rasterized, pr_multipolygons)
      }
      else {
        pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                  arr.ind=TRUE)[,2],]
        pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
        pr_restricted <- mask(pr_rasterized, pr_multipolygons)
      }
      
      ## if realized range does not cross equator, cut off potential range at the equator:
      ## DEC 14: commented out and rerun
      #lat_mp <- realized_range$lat_mp[num]
      
      # if (!realized_range$hemisphere[num] == "EQUATOR" & lat_mp > 0) {
      #   p <- Polygon(matrix(c(-180,0,-180,90,180,90,180,0,-180,0),
      #                       ncol=2, byrow=TRUE))
      #   rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
      # }
      # else if (!realized_range$hemisphere[num] == "EQUATOR" & lat_mp < 0) {
      #   p <- Polygon(matrix(c(-180,-90,-180,0,180,0,180,-90,-180,-90),
      #                       ncol=2, byrow=TRUE))
      #   rect <- SpatialPolygons(list(Polygons(list(p), "p1"))) 
      # }
      # else {
      #   p <- Polygon(matrix(c(-180, 90,-180,-90,180,-90,180, 90,-180, 90),
      #                       ncol=2, byrow=TRUE))
      #   rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
      # }
      # 
      # rect_raster <- rasterize(rect, r, getCover=TRUE)
      # rect_raster[rect_raster==0] <- NA
      
      #plot(pr_restricted, col = "orange")
      
      #pr_restricted <- mask(pr_restricted, rect_raster)
      
      ##plot(pr_restricted, col = "red", add = TRUE)
      
      ##plot(pr_restricted)
      ##plot(pr_multipolygons,  add=TRUE, col = "purple")
      ##plot(potential_range,  add=TRUE, col = "yellow")
      
      pr_restricted_all <- addLayer(pr_restricted_all, pr_restricted, updatevalue = NA)
      names <- append(names, paste(names(potential_raster), realized_range$source[num], sep = "_"))
      num = num + 1
    }
    
    i = i + 1
  }
  names(pr_restricted_all) <- names
  plot(pr_restricted_all)
  
  return(pr_restricted_all)
}


## returns a potential range raster containing a layer for all species in high_filtered and low_filtered that represents their potential equatorward and poleward range respectively containing only contiguous habitat in clumps that their realized range touches 
create_potential_ranges_one_limit <- function (clumped_temps, 
                                               realized_ranges, 
                                               high_filtered,  
                                               low_filtered) {
  ## WARM RANGES:
  pr_restricted_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                               crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  names <- c()
  i = 1
  if (nlayers(high_filtered) > 1) {
    while (i < length(names(high_filtered)) + 1) {
      ## get unrestricted potential range 
      potential_raster <- high_filtered[[i]] 
      potential_range <- potential_raster %>%
        clump(., directions = 8) %>%
        rasterToPolygons(., dissolve = TRUE) %>%
        st_as_sf()
      
      ## get realized range
      species <- names(potential_raster) %>%
        str_replace_all("_", " ")
      
      realized_range <- realized_ranges[which(as.character(realized_ranges$species)
                                                    %in% species),] 
      st_crs(realized_range) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"
      
      ## if species has multiple realized ranges (IUCN and GBIF), loop through them:
      num = 1
      while (num < nrow(realized_range)+1) {
        ## overlay realized range with potential range
        ## restrict potential range to clumps that overlap with the realized range
        intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
        rr <- filter(clumped_temps, intersects == TRUE)
        
        intersects_potential <- st_intersects(rr, potential_range, sparse = FALSE)[,]
        
        if (is.null(intersects_potential)) {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE),]
          pr_rasterized <- r
        }
        else if (nrow(rr) == 1) {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE),]
          pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
          pr_rasterized[pr_rasterized == 0] <- NA
          pr_restricted <- mask(pr_rasterized, pr_multipolygons)
        }
        else {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE)[,2],]
          pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
          pr_rasterized[pr_rasterized == 0] <- NA
          pr_restricted <- mask(pr_rasterized, pr_multipolygons)
        }
        
        ## if species range does not cross the equator, cut off potential range at the equator based on where realized range latitudinal midpoint is (below or above equator)
        if (realized_range$hemisphere[num] == "EQUATOR") {
          p <- Polygon(matrix(c(-180,-90,-180,90,180,90,180,-90,-180,-90),
                              ncol=2, byrow=TRUE))
          rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
        }
        else if (realized_range$hemisphere[num] == "N") {
          p <- Polygon(matrix(c(-180, 90,-180,0,180,0,180, 90,-180, 90),
                              ncol=2, byrow=TRUE))
          rect <- SpatialPolygons(list(Polygons(list(p), "p1"))) 
        }
        else {
          p <- Polygon(matrix(c(-180,-90,-180,0,180,0,180,-90,-180,-90),
                              ncol=2, byrow=TRUE))
          rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
        }
        
        rect_raster <- rasterize(rect, r, getCover=TRUE)
        rect_raster[rect_raster==0] <- NA
        
        plot(pr_restricted, col = "orange", main = paste(realized_range$species[num], "warm", 
                                                         realized_range$hemsphere[num]))
        
        pr_restricted <- mask(pr_restricted, rect_raster)
        
        plot(pr_restricted, add = TRUE)
        plot(rect, add = TRUE)
        ##plot(pr_multipolygons,  add=TRUE, col = "purple")
        ##plot(potential_range,  add=TRUE, col = "yellow")
        
        dev.copy(png, filename = paste("figures/one-limits/", realized_range$species[num], "_warm_",
                                       realized_range$source[num], ".png", sep = ""), 
                 width = 614, height = 472);
        dev.off()
        
        pr_restricted_high <- addLayer(pr_restricted_high, pr_restricted, updatevalue = NA)
        names <- append(names, paste(species, realized_range$source[num], "warm", sep = "_"))
        num = num + 1
      }
      
      i = i + 1
    }
    names(pr_restricted_high) <- names
    ##plot(pr_restricted_high)
  }
  
  
  ## COLD RANGES:
  pr_restricted_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  names <- c()
  i = 1
  if (!is.na(minValue(low_filtered)[1])) {
    while (i < length(names(low_filtered)) + 1) {
      ## get unrestricted potential range 
      potential_raster <- low_filtered[[i]] 
      potential_range <- potential_raster %>%
        clump(., directions = 8) %>%
        rasterToPolygons(., dissolve = TRUE) %>%
        st_as_sf()
      
      ## get realized range
      species <- names(potential_raster) %>%
        str_replace_all("_", " ")
      
      realized_range <- realized_ranges[which(as.character(realized_ranges$species)
                                                    %in% species),] 
      st_crs(realized_range) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"
      
      ## if species has multiple realized ranges (IUCN and GBIF), loop through them:
      num = 1
      while (num < nrow(realized_range)+1) {
        ## overlay realized range with potential range
        ## restrict potential range to clumps that overlap with the realized range
        intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
        rr <- filter(clumped_temps, intersects == TRUE)
        
        intersects_potential <- st_intersects(rr, potential_range, sparse = FALSE)[,]
        
        if (nrow(rr) == 1) {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE),]
          pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
          pr_rasterized[pr_rasterized == 0] <- NA
          pr_restricted <- mask(pr_rasterized, pr_multipolygons)
        }
        else {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE)[,2],]
          pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
          pr_rasterized[pr_rasterized == 0] <- NA
          pr_restricted <- mask(pr_rasterized, pr_multipolygons)
        }
        
        
        ## mark hemisphere since later equator species will need poleward range filling calculated twice
          plot(pr_restricted, col = "orange", main = paste(realized_range$species[num], "cold", 
                                                           realized_range$hemisphere[num]),
               legend = FALSE)
         
          ##plot(pr_multipolygons,  add=TRUE, col = "purple")
          ##plot(potential_range,  add=TRUE, col = "yellow")
          
          dev.copy(png, filename = paste("figures/one-limits/", realized_range$species[num], "_cold_",
                                         realized_range$source[num], ".png", sep = ""), 
                   width = 614, height = 472);
          dev.off()
          
          pr_restricted_low <- addLayer(pr_restricted_low, pr_restricted, updatevalue = NA)
          names <- append(names, paste(species, realized_range$source[num], "cold", 
                                       sep = "_"))

        
        num = num + 1
      }
      
      i = i + 1
    }
    names(pr_restricted_low) <- names
    plot(pr_restricted_low)
  }
  
  if(!is.na(minValue(high_filtered)[1]) & !is.na(minValue(low_filtered)[1])){
    pr_restricted_all <- addLayer(pr_restricted_high, pr_restricted_low)
  }
  else if(!is.na(minValue(low_filtered)[1]) & is.na(minValue(high_filtered)[1])) {
    pr_restricted_all <- pr_restricted_low
  }
  else {
    pr_restricted_all <- pr_restricted_high
  }
  
  return(pr_restricted_all)
}


## plots 
plot_ranges_overlap  <- function (potential_ranges) {
  potential_ranges <- readRDS("data-processed/potential_ranges_notcutatequator_dormancy.rds")
  rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")
  
  rrs <- as.data.frame(rasterized_rrs, xy=TRUE)
  colnames(rrs)[1:2] <- c("longitude", "latitude")
  
  prs <- as.data.frame(potential_ranges, xy=TRUE)
  colnames(prs)[1:2] <- c("longitude", "latitude")
  
  land <- as.data.frame(raster("./data-processed/raster_terr_mask.nc"), xy=TRUE)
  colnames(land)[1:3] <- c("longitude", "latitude", "mask")
  ocean <- as.data.frame(raster("./data-processed/raster_marine_mask.nc"), xy=TRUE)
  colnames(ocean)[1:3] <- c("longitude", "latitude", "mask")
  intertidal <- as.data.frame(raster("./data-processed/raster_intertidal_mask.nc"), xy=TRUE)
  colnames(intertidal)[1:3] <- c("longitude", "latitude", "mask")
  
  traits <- read.csv("./data-processed/wrangled-traits.csv")
  
  ## ggplots of all:
  i = 1
  while (i < ncol(prs) - 1) {
    range <- colnames(prs)[i+2]
    
    ## get realized range:
    species <- range %>%
      str_replace_all("_", ".")
    
    if(str_detect(species, "dormancy") ==  TRUE) {
      split <- str_split_fixed(range, "\\_", n = 7)
      source <- split[1,7]
    } 
    else {
      split <- str_split_fixed(range, "\\_", n = 3)
      source <- split[1,3]
    }
    
    ## get species name and realm
    species <- paste(split[1,1], split[1,2], sep = ".")
    realm <- traits$Realm[which(traits$genus_species == paste(split[1,1], split[1,2], sep = " "))]
    
    if(realm == "Terrestrial" | realm == "Freshwater") {
      backdrop = land
    }
    else if (realm == "Marine") {
      backdrop = ocean
    }
    else {
      backdrop = intertidal
    }
    
    rr_index <- which(str_detect(colnames(rrs), species) & str_detect(colnames(rrs), source))
    
    r <- rrs[,c(1:2, rr_index)] 
    colnames(r)[3] <- "rr"
    r <- left_join(r, prs[,c(1:2, which(colnames(prs) == range))]) 
    colnames(r)[4] <- "pr"
    
    r_gg <- r %>%
      ggplot(., aes(x = longitude, y = latitude)) +
      xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
      geom_raster(aes(fill=as.factor(rr))) + 
      scale_fill_manual(values = c("yellow"), aesthetics = 'fill', labels = ) +
      annotate(geom="raster", x=r$longitude, y=r$latitude, alpha=.5,
               fill = r$pr) +
      annotate(geom="raster", x=land$longitude, y=land$latitude, alpha=.1,
             fill = backdrop$mask) +
      labs(title = range,
           y = "Latitude",
           x = "Longitude") +
      scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
      scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) +
      theme(legend.position = "none")
    
    
    ## write to file:
    ggsave(r_gg, path = "figures/range-plots/", 
           filename = paste(range, ".png",sep = "_"), 
           height = 6, width = 10, units = "in", device = "png")
    
    i = i + 1
  }
  
}







### garbage ####
## writes plots of realized range and restricted potential range to files for species with both thermal limits 
plot_ranges <- function (clumped_temps, 
                         combined, 
                         realized_ranges) {
  i = 1
  while (i < length(names(combined)) + 1) {
    ## get unrestricted potential range 
    potential_range <- combined[[i]] 
    potential_raster <- potential_range %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf() 
    
    ## get realized range:
    species <- names(potential_range) %>%
      str_replace_all("_", " ")
    
    if(str_detect(species, "dormancy") ==  TRUE) {
      split <- str_split_fixed(species, " ", n = 3)
      
      species <- paste(split[1,1], split[1,2], sep = " ")
    }
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
    
    ##plot(st_geometry(realized_range), add = TRUE, col = "red")
    
    ## overlay realized range with potential range and restrict potential range to clumps that overlap with the realized range
    intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
    sub <- filter(clumped_temps, intersects == TRUE)
    
    intersects_potential <- st_intersects(sub, potential_raster, sparse = FALSE)[,]
    if (nrow(sub) == 1) {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE),]
    }
    else {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE)[,2],]
    }
    
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - realized range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(ms_simplify(realized_range, keep_shapes = TRUE)), add = TRUE, col = "red")
    
    dev.copy(png, filename = paste("figures/plot_ranges/", names(potential_range), "_realized-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - potential range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(sub_potential), add = TRUE, col = "yellow")
    
    dev.copy(png, filename = paste("figures/plot_ranges/", names(potential_range), "_potential-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    i = i + 1
  }
}


## writes plots of realized range and restricted potential poleward or equatorward range to files for species with only one thermal limit
plot_ranges_one_limit <- function (clumped_temps, 
                                   realized_ranges, 
                                   high_filtered, 
                                   low_filtered) {
  i = 1
  while (i < length(names(high_filtered)) + 1) {
    ## get unrestricted potential range 
    potential_range <- high_filtered[[i]] 
    potential_raster <- potential_range %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf() 
    
    ## get realized range:
    species <- names(potential_range) %>%
      str_replace_all("_", " ")
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
    
    ##plot(st_geometry(realized_range), add = TRUE, col = "red")
    
    ## overlay realized range with potential range and restrict potential range to clumps that overlap with the realized range
    intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
    sub <- filter(clumped_temps, intersects == TRUE)
    
    intersects_potential <- st_intersects(sub, potential_raster, sparse = FALSE)[,]
    if (nrow(sub) == 1) {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE),]
    }
    else {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE)[,2],]
    }
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - realized range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(ms_simplify(realized_range, keep_shapes = TRUE)), add = TRUE, col = "red")
    
    dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), "_realized-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - potential range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(sub_potential), add = TRUE, col = "yellow")
    
    dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), "_equatorward-potential-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    i = i + 1
  }
  
  i = 1
  while (i < length(names(low_filtered)) + 1 & length(names(low_filtered)) > 1) {
    ## get unrestricted potential range 
    potential_range <- low_filtered[[i]] 
    potential_raster <- potential_range %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf() 
    
    ## get realized range:
    species <- names(potential_range) %>%
      str_replace_all("_", " ")
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
    
    ##plot(st_geometry(realized_range), add = TRUE, col = "red")
    
    ## overlay realized range with potential range and restrict potential range to clumps that overlap with the realized range
    intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
    sub <- filter(clumped_temps, intersects == TRUE)
    plot(sub)
    
    intersects_potential <- st_intersects(sub, potential_raster, sparse = FALSE)[,]
    if (nrow(sub) == 1) {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE),]
    }
    else {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE)[,2],]
    }
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - realized range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(ms_simplify(realized_range, keep_shapes = TRUE)), add = TRUE, col = "red")
    
    dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), "_realized-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - potential range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(sub_potential), add = TRUE, col = "yellow")
    
    dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), "_poleward-potential-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    i = i + 1
  }
}





