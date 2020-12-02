## transitioning into the thermal dimension
library(tidyverse)
library(raster)
library(sf)
library(gridExtra)
select <- dplyr::select

## goal: to measure thermal range overfilling and underfilling in the thermal dimension

## this requires:
## 1. geting max and minimum temperatutes in the raster squares under each realized range (maybe need average values later?)
##    a. combine marine and terrestrial temperatute data into one layer 

r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1)

## read in seasonal high and low temp data:
terr_seasonal_high <- read.csv("data-processed/terrestrial_seasonal-max-temps.csv")
terr_seasonal_low <- read.csv("data-processed/terrestrial_seasonal-min-temps.csv")
marine_seasonal_high <- read.csv("data-processed/marine_seasonal-max-temps.csv") 
marine_seasonal_low <- read.csv("data-processed/marine_seasonal-min-temps.csv")

## rasterize:
raster_terr_high <- rasterize(terr_seasonal_high[, 1:2], r, terr_seasonal_high[,3], fun=mean) 
raster_terr_high[is.infinite(raster_terr_high)] <- NA
names(raster_terr_high) <- "seasonal_high_temp"
##plot(raster_terr_high, asp = 1)
raster_marine_high <- rasterize(marine_seasonal_high[, 1:2], r, marine_seasonal_high[,3], fun=mean)
raster_marine_high[is.infinite(raster_marine_high)] <- NA
names(raster_marine_high) <- "seasonal_high_temp"
##plot(raster_marine_high, asp = 1)

raster_terr_low <- rasterize(terr_seasonal_low[, 1:2], r, terr_seasonal_low[,3], fun=mean)
raster_terr_low[is.infinite(raster_terr_low)] <- NA
names(raster_terr_low) <- "seasonal_low_temp"
##plot(raster_terr_low, asp = 1)
raster_marine_low <- rasterize(marine_seasonal_low[, 1:2], r, marine_seasonal_low[,3], fun=mean)
raster_marine_low[is.infinite(raster_marine_low)] <- NA
names(raster_marine_low) <- "seasonal_low_temp"
##plot(raster_marine_low, asp = 1)

## merge marine high and low temps into a single raster layer:
## if cells overlap at the land-ocean boundary, keep the more extreme of the two values 
high_temps <- mosaic(raster_marine_high, raster_terr_high, fun = max)
##plot(high_temps, asp = 1)
low_temps <- mosaic(raster_marine_low, raster_terr_low, fun = min)
##plot(low_temps, asp = 1)



##    b. extract temperatures under each realized range for each species 
## read in rasterized realized ranges
rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")

## loop through ranges and extract temperatures under each 
i = 1
realized_temps <- c()
while (i < nlayers(rasterized_rrs) + 1) {
  range <- rasterized_rrs[[i]]
  
  ## use realized range raster as a mask to extract temperatures underneath
  highs <- mask(high_temps, range)
  if (length(which(is.na(values(highs)) == FALSE)) > 0) {
    high_vals <- data.frame(high_or_low = "high", temps = values(highs)[which(is.na(values(highs)) == FALSE)],
                            type = "realized")
    ## plot(highs)
  }
  else {
    high_vals <- data.frame(high_or_low = "high", temps = NA, type = "realized")
  }
  lows <- mask(low_temps, range)
  if (length(which(is.na(values(lows)) == FALSE)) > 0) {
    low_vals <- data.frame(high_or_low = "low", temps = values(lows)[which(is.na(values(lows)) == FALSE)],
                           type = "realized")
    ## plot(lows)
  }
  else {
    low_vals <- data.frame(high_or_low = "low", temps = NA,type = "realized")
  }
  
  if(i == 1) {
    realized_temps <- rbind(high_vals, low_vals) %>%
      mutate(range = names(rasterized_rrs)[i]) %>%
      select(range, everything())
  }
  else {
    realized_temps <- rbind(high_vals, low_vals) %>%
      mutate(range = names(rasterized_rrs)[i]) %>%
      select(range, everything()) %>%
      rbind(realized_temps, .)
  }
  i = i + 1
}

## save:
write.csv(realized_temps, "./data-processed/thermal-dimension_realized-temps.csv", row.names = FALSE)


## 2. comparing these temperatures to thermal tolerance limits in various ways 
##    a. for species with both thermal tolerance limits, get the temperatures 
##        in their potential range using the same method

## read in potential ranges:
potential_ranges <- readRDS("data-processed/potential_ranges.rds")

## loop through ranges and extract temperatures under each 
i = 1
potential_temps <- c()
while (i < nlayers(potential_ranges) + 1) {
  range <- potential_ranges[[i]]
  
  ## use realized range raster as a mask to extract temperatures underneath
  highs <- mask(high_temps, range)
  if (length(which(is.na(values(highs)) == FALSE)) > 0) {
    high_vals <- data.frame(high_or_low = "high", temps = values(highs)[which(is.na(values(highs)) == FALSE)],
                            type = "potential")
    ## plot(highs)
  }
  else {
    high_vals <- data.frame(high_or_low = "high", temps = NA, type = "potential")
  }
  lows <- mask(low_temps, range)
  if (length(which(is.na(values(lows)) == FALSE)) > 0) {
    low_vals <- data.frame(high_or_low = "low", temps = values(lows)[which(is.na(values(lows)) == FALSE)],
                           type = "potential")
    ## plot(lows)
  }
  else {
    low_vals <- data.frame(high_or_low = "low", temps = NA, type = "potential")
  }
  
  if(i == 1) {
    potential_temps <- rbind(high_vals, low_vals) %>%
      mutate(range = names(potential_ranges)[i]) %>%
      select(range, everything())
  }
  else {
    potential_temps <- rbind(high_vals, low_vals) %>%
      mutate(range = names(potential_ranges)[i]) %>%
      select(range, everything()) %>%
      rbind(potential_temps, .)
  }
  i = i + 1
}

## save:
write.csv(potential_temps, "./data-processed/thermal-dimension_potential-temps.csv", row.names = FALSE)



## plot them all: 
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = "."))


all_temps <- rbind(realized_temps, potential_temps)
num = 1 
while (num < length(unique(potential_temps$range)) + 1) {
  h <- all_temps %>%
    filter(range == unique(potential_temps$range)[num], high_or_low == "high")
  
  l <- all_temps %>%
    filter(range == unique(potential_temps$range)[num], high_or_low == "low")
  
  ## get thermal limits:
  lims <- thermal_limits[which(thermal_limits$genus_species == str_split_fixed(h$range[1], "_", n = 2)[1,1]),]
  ctmax <- lims$thermal_limit[which(lims$type == "max")]
  ctmin <- lims$thermal_limit[which(lims$type == "min")]
  
  hplot <- h %>%
    ggplot(., aes(x = temps, y =..density.., fill = type)) +  geom_density(alpha = 0.3) +
    labs(fill = "Range:", y = "Density", x = "Seasonal high temperature", title = h$range[1]) + 
    geom_vline(xintercept = ctmax, linetype="dotted", size = 0.5)
  
  lplot <- l %>%
    ggplot(., aes(x = temps, y =..density.., fill = type)) +  geom_density(alpha = 0.3) +
    labs(fill = "Range:", y = "Density", x = "Seasonal low temperature", title = l$range[1]) + 
    geom_vline(xintercept = ctmin, linetype="dotted", size = 0.5)
  
  hlplot <- grid.arrange(lplot, hplot, ncol = 2)
  
  ggsave(hlplot, path = "figures/thermal-dimension/", 
        filename = paste(l$range[1], ".png",sep = "_"), 
        height = 5, width = 11, units = "in", device = "png")
  
  num = num + 1
}
 



 



