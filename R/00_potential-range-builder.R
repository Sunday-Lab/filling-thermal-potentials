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

absolute <- nc_open("data-raw/absolute.nc")

lat = ncvar_get(absolute, "lat")
lon = ncvar_get(absolute, "lon")
time = ncvar_get(absolute, "time")
temps = ncvar_get(absolute, "tem")

nc_close(absolute)

## save rasters of each month
raster <- raster(xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
x = 1 
while (x  < 13) {
  slice <- temps[ , , x] 
  
  r <- raster(t(slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  raster <- addLayer(raster, r) 
  x = x+1
}

names(raster) <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")

## collapse into data frame where each row is a raster square 
tmp <- data.frame(rasterToPoints(raster)) 


## mean high temperature from the warmest month
high_tmp <- tmp %>%
  mutate(seasonal_high = pmax(.$January, .$February, .$March, .$April, .$May, 
                              .$June, .$July, .$August, .$September, .$October, .$November, .$December)) %>%
  select(x, y, seasonal_high) 


## mean low temperature from the coldest month 
low_tmp <- tmp %>%
  mutate(seasonal_low = pmin(.$January, .$February, .$March, .$April, .$May, 
                             .$June, .$July, .$August, .$September, .$October, .$November, .$December)) %>%
  select(x, y, seasonal_low)

## create raster layers of high and low seasonal temperature
r = raster(nrow = nrow(lat), ncol = nrow(lon),
           crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
raster_high <- rasterize(high_tmp[, 1:2], r, high_tmp[,3], fun=mean)
names(raster_high) <- "seasonal_high"
plot(raster_high, asp = 1)

raster_low <- rasterize(low_tmp[, 1:2], r, low_tmp[,3], fun=mean)
names(raster_low) <- "seasonal_low"
plot(raster_low, asp = 1)


## separate ocean data from land data
countries <- ne_countries(returnclass = "sp") 
countries <- spTransform(countries, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
land_high_tmp <- crop(raster_high, extent(countries)) %>%
  mask(., countries)
ocean_high_tmp <- crop(raster_high, extent(countries)) %>%
  mask(., countries, inverse = TRUE)

land_low_tmp <- crop(raster_low, extent(countries)) %>%
  mask(., countries)
ocean_low_tmp <- crop(raster_low, extent(countries)) %>%
  mask(., countries, inverse = TRUE)

plot(ocean_high_tmp)
plot(ocean_low_tmp)
plot(land_high_tmp)
plot(land_low_tmp)
  





####################################################################################
#####                   CREATING POTENTIAL RANGE SHAPEFILES                   ######
####################################################################################

## filter thermal limits to include only species we have realized ranges for 


##### TERRESTRIAL ##### 
## repeat for marine species after 
## read in thermal limit data for each species that has both thermal tolerance metrics 

thermal_limits <- read_csv("./data-raw/globtherm_full_dataset_2019.csv") %>%
  filter(thermy == "ectotherm")

upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Terrestrial")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Terrestrial")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]

names_high <- c("seasonal_high", paste(both_upper$Genus, both_upper$Species, sep = "_"))
names_low <- c("seasonal_low", paste(both_lower$Genus, both_lower$Species, sep = "_"))

## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
species = 1
while (species < nrow(both_upper) + 1) {
  land_high_tmp <- addLayer(land_high_tmp, land_high_tmp[[1]] - both_upper$thermal_limit[species]) 

  species = species + 1
}
names(land_high_tmp) <- names_high
plot(land_high_tmp)

species = 1
while (species < nrow(both_lower) + 1) {
  land_low_tmp <- addLayer(land_low_tmp, land_low_tmp[[1]] - both_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(land_low_tmp) <- names_low
plot(land_low_tmp)


## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
land_high_tmp[land_high_tmp > 0] <- NA
land_low_tmp[land_low_tmp < 0] <- NA

plot(land_high_tmp)
plot(land_low_tmp)


## combine: 
combined <- land_high_tmp[[1]]
i = 2  
while (i < nrow(both_upper) + 1) {
  combined <- addLayer(combined, mask(land_high_tmp[[i]], land_low_tmp[[i]]))
  
  i = i + 1
}

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

