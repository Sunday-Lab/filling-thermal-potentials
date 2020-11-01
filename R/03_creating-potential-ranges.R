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
#####                       SEASONAL TEMPERATURE RASTERS                      ######
####################################################################################
## read in seasonal high and low temp data:
terr_seasonal_high <- read.csv("data-processed/terrestrial_seasonal-max-temps.csv")
terr_seasonal_low <- read.csv("data-processed/terrestrial_seasonal-min-temps.csv")

## rasterize:
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

raster_terr_high <- rasterize(terr_seasonal_high[, 1:2], r, terr_seasonal_high[,3], fun=mean) 
raster_terr_high[is.infinite(raster_terr_high)] <- NA
names(raster_terr_high) <- "seasonal_high_temp"
plot(raster_terr_high, asp = 1)

raster_terr_low <- rasterize(terr_seasonal_low[, 1:2], r, terr_seasonal_low[,3], fun=mean)
raster_terr_low[is.infinite(raster_terr_low)] <- NA
names(raster_terr_low) <- "seasonal_low_temp"
plot(raster_terr_low, asp = 1)


## read in seasonal high and low temp data:
marine_seasonal_high <- read.csv("data-processed/marine_seasonal-max-temps.csv") 
marine_seasonal_low <- read.csv("data-processed/marine_seasonal-min-temps.csv")

## rasterize:
raster_marine_high <- rasterize(marine_seasonal_high[, 1:2], r, marine_seasonal_high[,3], fun=mean)
raster_marine_high[is.infinite(raster_marine_high)] <- NA
names(raster_marine_high) <- "seasonal_high_temp"
plot(raster_marine_high, asp = 1)

raster_marine_low <- rasterize(marine_seasonal_low[, 1:2], r, marine_seasonal_low[,3], fun=mean)
raster_marine_low[is.infinite(raster_marine_low)] <- NA
names(raster_marine_low) <- "seasonal_low_temp"
plot(raster_marine_low, asp = 1)


## create intertidal temperature data:
## create polygon representing the edge of land:
land <- raster_terr_high 
land[is.infinite(land)] = NA 
land[land > 0 | land < 0] <- 1

polygon <- land %>%
  rasterToPolygons(., dissolve = TRUE) %>% 
  st_as_sf()
smooth <- smooth(polygon, method = "ksmooth", smoothness = 10) 
st_crs(smooth) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")

## plot(st_geometry(smooth), axes = TRUE) 

## create buffer around intertidal area into sea and onto land
buffer_sea <- st_buffer(smooth, dist = 1)
buffer_land <- st_buffer(smooth, dist = -1)

##plot(st_geometry(buffer_sea))
##plot(st_geometry(buffer_land))

## subset temperature data to include only temperatures in buffer 
intertidal_sea_high <- raster_marine_high %>%
  mask(., buffer_sea) 

intertidal_land_high <- raster_terr_high %>%
  mask(., buffer_land, inverse = TRUE)

intertidal_sea_low <- raster_marine_low %>%
  mask(., buffer_sea) 

intertidal_land_low <- raster_terr_low %>%
  mask(., buffer_land, inverse = TRUE)

## combine land and sea temperatures, giving priority to sea temperatures 
raster_intertidal_high <- merge(intertidal_sea_high, intertidal_land_high)
raster_intertidal_low <- merge(intertidal_sea_low, intertidal_land_low)

plot(raster_intertidal_high)
plot(raster_intertidal_low)









####################################################################################
#####                   CREATING POTENTIAL RANGE SHAPEFILES                   ######
####################################################################################

## read in thermal limits:
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = " "))

## read in realized ranges:
realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp")

############################################
#####   TWO THERMAL LIMIT RANGES      ######
############################################
## TERRESTRIAL AND FRESHWATER SPECIES:
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]

## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
terr_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                    crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
terr_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
species = 1
while (species < nrow(both_upper) + 1) {
  terr_high <- addLayer(terr_high, raster_terr_high[[1]] - both_upper$thermal_limit[species]) 
  terr_low <- addLayer(terr_low, raster_terr_low[[1]] - both_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(terr_high) <- c(paste(both_upper$Genus, both_upper$Species, sep = "_"))
names(terr_low) <- c(paste(both_lower$Genus, both_lower$Species, sep = "_"))
plot(terr_high)
plot(terr_low)


## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
terr_high[terr_high > 0] <- NA
terr_low[terr_low < 0] <- NA
names(terr_high) <- c(paste(both_upper$Genus, both_upper$Species, sep = "_"))
names(terr_low) <- c(paste(both_lower$Genus, both_lower$Species, sep = "_"))
plot(terr_high)
plot(terr_low)


## combine to find cells where seasonal high temp is less than CTmax and seasonal low temp is greater than CTmin 
combined <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
i = 1  
while (i < nrow(both_upper) + 1) {
  ## "updatevalue = NA" sets cells where one thermal limit is exceeded to NA
  combined <- addLayer(combined, mask(terr_high[[i]], terr_low[[i]]), updatevalue = NA)
  
  i = i + 1
}
names(combined) <- paste(both_upper$Genus, both_upper$Species, sep = "_")
plot(combined)

## clump contiguous habitat together:
clumped_temps <- raster_terr_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## restrict range to contiguous habitat that begins at the species realized range:
plot_ranges(clumped_temps = clumped_temps, realized_ranges = realized_ranges, combined = combined)
prs_terrestrial <- create_potential_ranges(clumped_temps = clumped_temps, 
                                           realized_ranges = realized_ranges, combined = combined)



## MARINE SPECIES:
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Marine")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Marine")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]

## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
marine_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                      crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
marine_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
species = 1
while (species < nrow(both_upper) + 1) {
  marine_high <- addLayer(marine_high, raster_marine_high[[1]] - both_upper$thermal_limit[species]) 
  marine_low <- addLayer(marine_low, raster_marine_low[[1]] - both_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(marine_high) <- c(paste(both_upper$Genus, both_upper$Species, sep = "_"))
names(marine_low) <- c(paste(both_lower$Genus, both_lower$Species, sep = "_"))
plot(marine_high)
plot(marine_low)


## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
marine_high[marine_high > 0] <- NA
marine_low[marine_low < 0] <- NA
plot(marine_high)
plot(marine_low)


## combine to find cells where seasonal high temp is less than CTmax and seasonal low temp is greater than CTmin 
combined <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
i = 1  
while (i < nrow(both_upper) + 1) {
  ## "updatevalue = NA" sets cells where one thermal limit is exceeded to NA
  combined <- addLayer(combined, mask(marine_high[[i]], marine_low[[i]]), updatevalue = NA)
  
  i = i + 1
}
names(combined) <- paste(both_upper$Genus, both_upper$Species, sep = "_")
plot(combined)

## clump contiguous habitat together:
clumped_temps <- raster_marine_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## restrict range to contiguous habitat that begins at the species realized range:
plot_ranges(clumped_temps = clumped_temps, realized_ranges = realized_ranges, combined = combined)
prs_marine <- create_potential_ranges(clumped_temps = clumped_temps, 
                                      realized_ranges = realized_ranges, combined = combined)


## INTERTIDAL SPECIES:
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Intertidal")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Intertidal")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]

## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
intertidal_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                          crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
intertidal_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                         crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
species = 1
while (species < nrow(both_upper) + 1) {
  intertidal_high <- addLayer(intertidal_high, raster_intertidal_high[[1]] - both_upper$thermal_limit[species]) 
  intertidal_low <- addLayer(intertidal_low, raster_intertidal_low[[1]] - both_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(intertidal_high) <- c(paste(both_upper$Genus, both_upper$Species, sep = "_"))
names(intertidal_low) <- c(paste(both_lower$Genus, both_lower$Species, sep = "_"))
plot(intertidal_high)
plot(intertidal_low)


## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
intertidal_high[intertidal_high > 0] <- NA
intertidal_low[intertidal_low < 0] <- NA
plot(intertidal_high)
plot(intertidal_low)


## combine to find cells where seasonal high temp is less than CTmax and seasonal low temp is greater than CTmin 
combined <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
i = 1  
while (i < nrow(both_upper) + 1) {
  ## "updatevalue = NA" sets cells where one thermal limit is exceeded to NA
  combined <- addLayer(combined, mask(intertidal_high[[i]], intertidal_low[[i]]), updatevalue = NA)
  
  i = i + 1
}
names(combined) <- paste(both_upper$Genus, both_upper$Species, sep = "_")
plot(combined)

## clump contiguous habitat together:
clumped_temps <- raster_intertidal_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps$geometry[1], col = "red")

## restrict range to contiguous habitat that begins at the species realized range:
plot_ranges(clumped_temps = clumped_temps, realized_ranges = realized_ranges, combined = combined)
prs_intertidal <- create_potential_ranges(clumped_temps = clumped_temps, 
                                          realized_ranges = realized_ranges, combined = combined)



######################################################
## combine all potential ranges into single object:
potential_ranges <- stack(prs_terrestrial, prs_marine, prs_intertidal)
##saveRDS(potential_ranges, "data-processed/potential_ranges.rds")


## investigate
i=1
while (i < nlayers(potential_ranges) + 1) {
  plot(potential_ranges[[i]], main = names(potential_ranges)[i])
  
  i = i+1
}




############################################
#####   ONE THERMAL LIMIT RANGES      ######
############################################
## TERRESTRIAL AND FRESHWATER SPECIES:
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

only_upper <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]

## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
terr_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                    crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
terr_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
## for upper limits, only filter use seasonal high temps 
species = 1
while (species < nrow(only_upper) + 1) {
  terr_high <- addLayer(terr_high, raster_terr_high[[1]] - only_upper$thermal_limit[species]) 
  
  species = species + 1
}
## for lower limits, only filter use seasonal low temps 
species = 1
while (species < nrow(only_lower) + 1) {
  terr_low <- addLayer(terr_low, raster_terr_low[[1]] - only_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(terr_high) <- c(paste(only_upper$Genus, only_upper$Species, sep = "_"))
names(terr_low) <- c(paste(only_lower$Genus, only_lower$Species, sep = "_"))
plot(terr_high)
plot(terr_low)

## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
terr_high[terr_high > 0] <- NA
terr_low[terr_low < 0] <- NA
names(terr_high) <- c(paste(only_upper$Genus, only_upper$Species, sep = "_"))
names(terr_low) <- c(paste(only_lower$Genus, only_lower$Species, sep = "_"))
plot(terr_high)
plot(terr_low)

## clump contiguous habitat together:
clumped_temps <- raster_terr_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## restrict range to contiguous habitat that begins at the species realized range:
plot_ranges_one_limit(clumped_temps = clumped_temps, realized_ranges = realized_ranges, 
                      high_filtered = terr_high,  low_filtered = terr_low)
prs_terrestrial_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
                                                     realized_ranges = realized_ranges, 
                                                     high_filtered = terr_high,  
                                                     low_filtered = terr_low)

## MARINE SPECIES:
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Marine")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Marine")

only_upper <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]

## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
marine_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                    crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
marine_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
## for upper limits, only filter use seasonal high temps 
species = 1
while (species < nrow(only_upper) + 1) {
  marine_high <- addLayer(marine_high, raster_marine_high[[1]] - only_upper$thermal_limit[species]) 
  
  species = species + 1
}
## for lower limits, only filter use seasonal low temps 
species = 1
while (species < nrow(only_lower) + 1) {
  marine_low <- addLayer(marine_low, raster_marine_low[[1]] - only_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(marine_high) <- c(paste(only_upper$Genus, only_upper$Species, sep = "_"))
names(marine_low) <- c(paste(only_lower$Genus, only_lower$Species, sep = "_"))
plot(marine_high)
plot(marine_low)

## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
marine_high[marine_high > 0] <- NA
marine_low[marine_low < 0] <- NA
names(marine_high) <- c(paste(only_upper$Genus, only_upper$Species, sep = "_"))
names(marine_low) <- c(paste(only_lower$Genus, only_lower$Species, sep = "_"))
plot(marine_high)
plot(marine_low)

## clump contiguous habitat together:
clumped_temps <- raster_marine_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## restrict range to contiguous habitat that begins at the species realized range:
plot_ranges_one_limit(clumped_temps = clumped_temps, realized_ranges = realized_ranges, 
                      high_filtered = marine_high,  low_filtered = marine_low)
prs_marine_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
                                                               realized_ranges = realized_ranges, 
                                                               high_filtered = marine_high,  
                                                               low_filtered = marine_low)


## INTERTIDAL SPECIES:
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Intertidal")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Intertidal")

only_upper <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]

## create an individual raster layer of difference between thermal limit and seasonal temperature for
intertidal_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                      crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
intertidal_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
## for upper limits, only filter use seasonal high temps 
species = 1
while (species < nrow(only_upper) + 1) {
  intertidal_high <- addLayer(intertidal_high, raster_intertidal_high[[1]] - only_upper$thermal_limit[species]) 
  
  species = species + 1
}
## for lower limits, only filter use seasonal low temps 
species = 1
while (species < nrow(only_lower) + 1) {
  intertidal_low <- addLayer(intertidal_low, raster_intertidal_low[[1]] - only_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(intertidal_high) <- c(paste(only_upper$Genus, only_upper$Species, sep = "_"))
names(intertidal_low) <- c(paste(only_lower$Genus, only_lower$Species, sep = "_"))
plot(intertidal_high)
plot(intertidal_low)

## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
intertidal_high[intertidal_high > 0] <- NA
intertidal_low[intertidal_low < 0] <- NA
names(intertidal_high) <- c(paste(only_upper$Genus, only_upper$Species, sep = "_"))
names(intertidal_low) <- c(paste(only_lower$Genus, only_lower$Species, sep = "_"))
plot(intertidal_high)
plot(intertidal_low)

## clump contiguous habitat together:
clumped_temps <- raster_intertidal_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## restrict range to contiguous habitat that begins at the species realized range:
plot_ranges_one_limit(clumped_temps = clumped_temps, realized_ranges = realized_ranges, 
                      high_filtered = intertidal_high,  low_filtered = intertidal_low)
prs_intertidal_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
                                                          realized_ranges = realized_ranges, 
                                                          high_filtered = intertidal_high,  
                                                          low_filtered = intertidal_low)


potential_ranges_one_limit <- stack(prs_terrestrial_one_limit, prs_marine_one_limit, 
                                    prs_intertidal_one_limit)

##saveRDS(potential_ranges_one_limit, "data-processed/potential_ranges_one_limit.rds")



######################################################
##                  FUNCTIONS:                      ##
######################################################

## returns a potential range raster containing a layer for all species in combined that represents their potential range containing only contiguous habitat in clumps that their realized range touches 
## returns a potential range raster containing a layer for all species in combined that represents their potential range containing only contiguous habitat in clumps that their realized range touches 
create_potential_ranges <- function (clumped_temps, combined, realized_ranges) {
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
      
      ##plot(pr_restricted)
      ##plot(pr_multipolygons,  add=TRUE, col = "purple")
      ##plot(potential_range,  add=TRUE, col = "yellow")
      
      pr_restricted_all <- addLayer(pr_restricted_all, pr_restricted, updatevalue = NA)
      names <- append(names, paste(species, realized_range$source[num], sep = "_"))
      num = num + 1
    }
    
    i = i + 1
  }
  names(pr_restricted_all) <- names
  plot(pr_restricted_all)
  
  return(pr_restricted_all)
}


## returns a potential range raster containing a layer for all species in high_filtered and low_filtered that represents their potential equatorward and poleward range respectively containing only contiguous habitat in clumps that their realized range touches 
create_potential_ranges_one_limit <- function (clumped_temps = clumped_temps, 
                                               realized_ranges = realized_ranges, 
                                               high_filtered = terr_high,  
                                               low_filtered = terr_low) {
  ## EQUATORWARD RANGES:
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
        
        ##plot(pr_restricted)
        ##plot(pr_multipolygons,  add=TRUE, col = "purple")
        ##plot(potential_range,  add=TRUE, col = "yellow")
        
        ## restrict potential range to equator if species range does not cross the equator 
        
        pr_restricted_high <- addLayer(pr_restricted_high, pr_restricted, updatevalue = NA)
        names <- append(names, paste(species, realized_range$source[num], "ONLY_UPPER", sep = "_"))
        num = num + 1
      }
      
      i = i + 1
    }
    names(pr_restricted_high) <- names
    plot(pr_restricted_high)
  }
  
  
  ## POLEWARD RANGES:
  pr_restricted_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  names <- c()
  i = 1
  if (nlayers(low_filtered) > 1) {
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
        
        ##plot(pr_restricted)
        ##plot(pr_multipolygons,  add=TRUE, col = "purple")
        ##plot(potential_range,  add=TRUE, col = "yellow")
        
        pr_restricted_low <- addLayer(pr_restricted_low, pr_restricted, updatevalue = NA)
        names <- append(names, paste(species, realized_range$source[num], "ONLY_LOWER", sep = "_"))
        num = num + 1
      }
      
      i = i + 1
    }
    names(pr_restricted_low) <- names
    plot(pr_restricted_low)
  }
  
  if(nlayers(pr_restricted_low) > 1 & nlayers(pr_restricted_high) > 1){
    pr_restricted_all <- addLayer(pr_restricted_high, pr_restricted_low)
  }
  else if(nlayers(pr_restricted_low) > 1) {
    pr_restricted_all <- pr_restricted_low
  }
  else {
    pr_restricted_all <- pr_restricted_high
  }
  
  return(pr_restricted_all)
}


## writes plots of realized range and restricted potential range to files for species with both thermal limits 
plot_ranges <- function (clumped_temps, combined, realized_ranges) {
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
    
    dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), "_potential-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    i = i + 1
  }
}

## writes plots of realized range and restricted potential poleward or equatorward range to files for species with only one thermal limit
plot_ranges_one_limit <- function (clumped_temps, realized_ranges, high_filtered, low_filtered) {
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












## garbage:
###########
land_only <- nc_open("data-raw/HadCRUT.4.6.0.0.median.nc")

lat = ncvar_get(land_only, "latitude")
lon = ncvar_get(land_only, "longitude")
temps_land = ncvar_get(land_only, "temperature_anomaly")

nc_close(land_only)

## figure out which coordinates are on land vs the ocean:
map_df <- crossing(lat, lon) %>%
  arrange(-lat) %>%
  mutate(temp = as.vector(t(temps_land[ , , 1])))

x = 2
while(x < 2048) {
  layer <- as.vector(t(temps_land[ , , x]))
  map_df <- map_df %>%
    mutate(new_temp = layer) %>%
    mutate(temp = ifelse(!is.na(new_temp), new_temp, temp))
  x = x+1
}

## any cells that remain NA have no temperature data and therefore are ocean 
map_df <- map_df %>%
  mutate(is_land = !is.na(map_df$temp)) %>%
  select(lat, lon, is_land) 


## create new ocean raster (NA in land_only) and land raster(!NA in land_only)
land <- subset(map_df, is_land == TRUE) 
ocean <-  subset(map_df, is_land == FALSE) 




## raster to polygon code:
combinedtest <- combined[[-1]]
plot(combinedtest)

copy <- combinedtest
copy[copy < 0] = 1 

plot(copy[[11]], axes = TRUE) 

polygon <- copy[[11]] %>%
  rasterToPolygons(., dissolve = TRUE) %>% 
  st_as_sf()
smooth <- smooth(polygon, method = "ksmooth", smoothness = 2) 

st_crs(smooth) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")

plot(smooth, axes = TRUE) 

upper <- st_crop(smooth, xmin = -180, xmax = 180, ymin = 10, ymax = 85)

st_crop(smooth, xmin = -180, xmax = 180, ymax = 10, ymin = -90)

plot(upper, axes = TRUE) 
plot(lower, axes = TRUE) 

upper_raster <- land_high_tmp[[12]]
upper_raster[upper_raster < 0] = 1 
upper_raster <- crop(upper_raster, y = c(-180, 180, -30, 15))

lower_raster <- copy[[11]]
lower_raster <- crop(lower_raster, y = c(-180, 180, -90, -30))


plot(upper_raster)


## write code to restrict clumps to those crossing a set of lines (will be midpoints)
## and then clumps within xx km of those clumps 

potential <- combined[[2]]
plot(potential)

## create multilinestring representing latitudinal midpoints of range polygons 
multiline <- st_multilinestring(list(rbind(c(-170,0),c(-115,0)), 
                                     rbind(c(75,15),c(89,15)), 
                                     rbind(c(10,-85),c(-10,-85)),
                                     rbind(c(-170,-10),c(110,-10))))

polyraster <- clump(potential, directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()

plot(st_geometry(polyraster))
plot(multiline, add = TRUE) 

intersects <- st_intersects(polyraster, range, sparse = FALSE)[,]
sub <- filter(polyraster, intersects == TRUE)

plot(sub)
plot(st_geometry(multiline), add = TRUE, col = "RED")


## checking out how CTmax vs CTmin affects potential range:
i = 1
while (i < nlayers(terr_low)) {
  sp <- str_replace(names(terr_high)[i], "_", " ")
  rr <- realized_ranges[which(as.character(realized_ranges$species) %in% sp),] %>%
    st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
  
  par(mfrow=c(1,2))
  plot(terr_high[[i]], main = paste(names(terr_high)[i], "\n Seasonal high", sep = ""))
  plot(st_geometry(rr), add = TRUE, col = "red")
  plot(terr_low[[i]], main = "Seasonal low")
  plot(st_geometry(rr), add = TRUE, col = "red")
  
  dev.copy(png, filename = paste("figures/seasonal-temp-filtering/", names(terr_high)[i], sep = ""), width = 1000, height =600, res = 150);
  dev.off() 
  
  i = i + 1
}


## investigating clumping:
clumped_temps <- raster_intertidal_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps, col = colorRampPalette(brewer.pal(8, "Accent"))(48))

clumped_temps <- raster_marine_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps, col =(brewer.pal(4, "Accent")))

clumped_temps <- raster_terr_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf(bbox = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90)) 
plot(clumped_temps, col =  colorRampPalette(brewer.pal(8, "Accent"))(61))
