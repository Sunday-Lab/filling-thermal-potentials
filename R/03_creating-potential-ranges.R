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


raster_terr_high <- rasterize(terr_seasonal_high[, 1:2], r, terr_seasonal_high[,3], fun=mean) 
raster_terr_high[is.infinite(raster_terr_high)] <- NA
names(raster_terr_high) <- "seasonal_high_temp"
##plot(raster_terr_high, asp = 1)

raster_terr_low <- rasterize(terr_seasonal_low[, 1:2], r, terr_seasonal_low[,3], fun=mean)
raster_terr_low[is.infinite(raster_terr_low)] <- NA
names(raster_terr_low) <- "seasonal_low_temp"
##plot(raster_terr_low, asp = 1)


## read in seasonal high and low temp data:
marine_seasonal_high <- read.csv("data-processed/marine_seasonal-max-temps.csv") 
marine_seasonal_low <- read.csv("data-processed/marine_seasonal-min-temps.csv")

## rasterize:
raster_marine_high <- rasterize(marine_seasonal_high[, 1:2], r, marine_seasonal_high[,3], fun=mean)
raster_marine_high[is.infinite(raster_marine_high)] <- NA
names(raster_marine_high) <- "seasonal_high_temp"
##plot(raster_marine_high, asp = 1)

raster_marine_low <- rasterize(marine_seasonal_low[, 1:2], r, marine_seasonal_low[,3], fun=mean)
raster_marine_low[is.infinite(raster_marine_low)] <- NA
names(raster_marine_low) <- "seasonal_low_temp"
##plot(raster_marine_low, asp = 1)


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
##plot(raster_intertidal_high)
##plot(raster_intertidal_low)




####################################################################################
#####                   CREATING POTENTIAL RANGE SHAPEFILES                   ######
####################################################################################
## read in thermal limits:
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = " "))

## read in realized ranges:
realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp")
split_realized_ranges <- st_read("data-processed/realized-ranges_split.shp")

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
clumped_temps <- raster_terr_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## SPEICES WITH BOTH THERMAL LIMITS:
####################################=
combined = filter_by_tolerance_both(both_upper = both_upper, 
                                    both_lower = both_lower, 
                                    raster_high = raster_terr_high, 
                                    raster_low = raster_terr_low)
plot(combined)

## restrict range to contiguous habitat that begins at the species realized range:
##plot_ranges(clumped_temps = clumped_temps, realized_ranges = realized_ranges, combined = combined)
prs_terrestrial <- create_potential_ranges(clumped_temps = clumped_temps, 
                                           realized_ranges = realized_ranges, combined = combined)

## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################
high_filtered <- filter_by_tolerance_upper(only_upper = only_upper, 
                                           raster_high = raster_terr_high)
plot(high_filtered)
low_filtered <- filter_by_tolerance_lower(only_lower = only_lower, 
                                          raster_low = raster_terr_low)
plot(low_filtered)

## restrict range to contiguous habitat that begins at the species realized range:
## plot_ranges_one_limit(clumped_temps = clumped_temps, realized_ranges = realized_ranges,high_filtered = terr_high,  low_filtered = terr_low)
prs_terrestrial_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
                                                               realized_ranges = realized_ranges,
                                                               high_filtered = high_filtered,  
                                                               low_filtered = low_filtered)



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
clumped_temps <- raster_marine_high %>%
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
##plot_ranges(clumped_temps = clumped_temps, realized_ranges = realized_ranges, combined = combined)
prs_marine <- create_potential_ranges(clumped_temps = clumped_temps, 
                                      realized_ranges = realized_ranges, combined = combined)

## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################
high_filtered <- filter_by_tolerance_upper(only_upper = only_upper,
                                           raster_high = raster_marine_high)
plot(high_filtered)
low_filtered <- filter_by_tolerance_lower(only_lower = only_lower,
                                          raster_low = raster_marine_low)
plot(low_filtered)

## restrict range to contiguous habitat that begins at the species realized range:
##plot_ranges_one_limit(clumped_temps = clumped_temps, realized_ranges = realized_ranges, high_filtered = marine_high,  low_filtered = marine_low)
prs_marine_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
                                                          realized_ranges = realized_ranges, 
                                                          high_filtered = high_filtered,  
                                                          low_filtered = low_filtered)


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
clumped_temps <- raster_intertidal_high %>%
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
##plot_ranges(clumped_temps = clumped_temps, realized_ranges = realized_ranges, combined = combined)
prs_intertidal <- create_potential_ranges(clumped_temps = clumped_temps, 
                                          realized_ranges = realized_ranges, combined = combined)

## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################
high_filtered <- filter_by_tolerance_upper(only_upper = only_upper,
                                           raster_high = raster_intertidal_high)
plot(high_filtered)
low_filtered <- filter_by_tolerance_lower(only_lower = only_lower,
                                          raster_low = raster_intertidal_low)
plot(low_filtered)

## restrict range to contiguous habitat that begins at the species realized range:
##plot_ranges_one_limit(clumped_temps = clumped_temps, realized_ranges = realized_ranges, high_filtered = intertidal_high,  low_filtered = intertidal_low)
prs_intertidal_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
                                                              realized_ranges = realized_ranges, 
                                                              high_filtered = high_filtered,  
                                                              low_filtered = low_filtered)


########################################################
##      COLLATING POTENTIAL RANGES ACROSS REALMS     ###
########################################################
potential_ranges_both_limits <- stack(prs_terrestrial, prs_marine, prs_intertidal)
##saveRDS(potential_ranges_both_limits, "data-processed/potential_ranges.rds")
writeRaster(potential_ranges_both_limits, "data-processed/potential_ranges_both_limits.nc")

potential_ranges_one_limit <- stack(prs_terrestrial_one_limit, prs_marine_one_limit, 
                                    prs_intertidal_one_limit)
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
      lat_mp <- realized_range$lat_mp[num]
      
      if (!realized_range$hemisphere[num] == "EQUATOR" & lat_mp > 0) {
        p <- Polygon(matrix(c(-180,0,-180,90,180,90,180,0,-180,0),
                            ncol=2, byrow=TRUE))
        rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
      }
      else if (!realized_range$hemisphere[num] == "EQUATOR" & lat_mp < 0) {
        p <- Polygon(matrix(c(-180,-90,-180,0,180,0,180,-90,-180,-90),
                            ncol=2, byrow=TRUE))
        rect <- SpatialPolygons(list(Polygons(list(p), "p1"))) 
      }
      else {
        p <- Polygon(matrix(c(-180, 90,-180,-90,180,-90,180, 90,-180, 90),
                            ncol=2, byrow=TRUE))
        rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
      }
      
      rect_raster <- rasterize(rect, r, getCover=TRUE)
      rect_raster[rect_raster==0] <- NA
      
      #3plot(pr_restricted, col = "orange")
      
      pr_restricted <- mask(pr_restricted, rect_raster)
      
      ##plot(pr_restricted, col = "red", add = TRUE)
      
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


