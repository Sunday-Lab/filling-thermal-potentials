## script to create potential range shapefiles for each species based on their thermal tolerance limits
library(tidyverse)
library(sp)
library(sf)
library(raster)
library(ncdf4)
library(rnaturalearth)
library(smoothr)
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
buffer_sea <- st_buffer(smooth, dist = 2)
buffer_land <- st_buffer(smooth, dist = -2)

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
## realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp")


## TERRESTRIAL AND FRESHWATER SPECIES:
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]
only_upper<- upper_limits[!upper_limits$genus_species %in% both_upper$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% both_lower$genus_species,]


## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
species = 1
while (species < nrow(both_upper) + 1) {
  raster_terr_high <- addLayer(raster_terr_high, raster_terr_high[[1]] - both_upper$thermal_limit[species]) 
  raster_terr_low <- addLayer(raster_terr_low, raster_terr_low[[1]] - both_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(raster_terr_high) <- c("seasonal_high", paste(both_upper$Genus, both_upper$Species, sep = "_"))
names(raster_terr_low) <- c("seasonal_low", paste(both_lower$Genus, both_lower$Species, sep = "_"))
plot(raster_terr_high)
plot(raster_terr_low)


## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
raster_terr_high <- raster_terr_high[[-1]]
raster_terr_low <- raster_terr_low[[-1]]
raster_terr_high[raster_terr_high > 0] <- NA
raster_terr_low[raster_terr_low < 0] <- NA

plot(raster_terr_high)
plot(raster_terr_low)


## combine to find cells where seasonal high temp is less than CTmax and seasonal low temp is greater than CTmin 
combined <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
i = 1  
while (i < nrow(both_upper) + 1) {
  ## "updatevalue = NA" sets cells where one thermal limit is exceeded to NA
  combined <- addLayer(combined, mask(raster_terr_high[[i]], raster_terr_low[[i]]), updatevalue = NA)
  
  i = i + 1
}
names(combined) <- paste(both_upper$Genus, both_upper$Species, sep = "_")
plot(combined)



## restrict range to contiguous habitat that begins at the species realized range:
i = 1
while (i < nrow(both_upper) + 1) {
  ## get unrestricted potential range and clump contiguous habitat together:
  potential_range <- combined[[i]] 
  clumped_pr <- potential_range %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf()
  
  ##plot(potential_range)
  plot(st_geometry(clumped_pr), main = names(potential_range))
  
  ## get realized range:
  species <- names(potential_range) %>%
    str_replace_all("_", " ")
  
  realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),]
  
  plot(st_geometry(realized_range), add = TRUE, col = "red")
  
  ## overlay realized range with potential range and restrict potential range to clumps that overlap with the realized range
  intersects <- st_intersects(clumped_pr, realized_range, sparse = FALSE)[,]
  sub <- filter(clumped_pr, intersects == TRUE)
  
  plot(sub, add = TRUE, col = "blue")

  
  dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), ".png", sep = ""), width = 1000, height = 500, res = 200);
  dev.off()
  
  i = i + 1
}

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





####################################################################################
#####                            REALIZED RANGE POLYGONS                      ######
####################################################################################
library(tmap)
polygons <- st_read("data-raw/polygons/Filtered occurences ectotherm animals_020817.shp")


head(polygons)
qtm(polygons)
















## garbage:
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

<- st_crop(smooth, xmin = -180, xmax = 180, ymax = 10, ymin = -90)

plot(upper, axes = TRUE) 
plot(lower, axes = TRUE) 

upper_raster <- land_high_tmp[[12]]
upper_raster[upper_raster < 0] = 1 
upper_raster <- crop(upper_raster, y = c(-180, 180, -30, 15))

lower_raster <- copy[[11]]
lower_raster <- crop(lower_raster, y = c(-180, 180, -90, -30))


plot(upper_raster)
