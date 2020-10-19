library(ncdf4)
library(tidyverse)
library(raster)
select <- dplyr::select


####################################################################################
#####               TERRESTRIAL SEASONAL TEMPERATURE HIGH AND LOW             ######
####################################################################################

## max and min terrestrial temps:
filename <- paste("/Volumes/ADATA HV620/temperature-data/terrestrial/Complete_TMAX_Daily_LatLong1_1950.nc", sep = "")
ncfile <- nc_open(filename)

## create variables for things needed to use data
lat <- ncvar_get(ncfile, "latitude")
long <- ncvar_get(ncfile, "longitude")

## close the file
nc_close(ncfile)

## create arrays to hold mean temperatures from each 10 year dataset
mean_max <- array(dim = c(360, 180, 3650))
mean_max[,,] <- NaN
mean_min <- array(dim = c(360, 180, 3650))
mean_min[,,] <- NaN

rep = 1950
while (rep < 2020) {
  ## open max and minimum files 
  filename_max <- paste("/Volumes/ADATA HV620/temperature-data/terrestrial/Complete_TMAX_Daily_LatLong1_", rep, ".nc", sep = "")
  ncfile_max <- nc_open(filename_max)
  filename_min <- paste("/Volumes/ADATA HV620/temperature-data/terrestrial/Complete_TMIN_Daily_LatLong1_", rep, ".nc", sep = "")
  ncfile_min <- nc_open(filename_min)
  
  ## create variables for data
  date <- ncvar_get(ncfile_max, "date_number")
  arr.anom_max <- ncvar_get(ncfile_max, "temperature")
  arr.clim_max <- ncvar_get(ncfile_max, "climatology")
  arr.anom_min <- ncvar_get(ncfile_min, "temperature")
  arr.clim_min <- ncvar_get(ncfile_min, "climatology")
  
  ## close files
  nc_close(ncfile_max)
  nc_close(ncfile_min)
  
  ## figure out which years in this time frame are leap years
  leap_years <- seq(from = rep, to = (rep+9), by = 1) %% 4 == 0
  
  ## figure out which index represents the end of each year 
  x <- c(1)
  i = 1
  while(i < length(leap_years)) {
    if(leap_years[i] == FALSE) {
      x <- append(x, x[i] + 365)
    }
    else {
      x <- append(x, x[i] + 366)
    }
    i = i+1
  }
  
  leap_years <- data.frame(leap_year = leap_years, index = x)
  
  ## create arrays to store temperatures in 
  arr.temps_max <- array(dim = c(nrow(arr.clim_max), ncol(arr.clim_max), 3650))
  arr.temps_min <- array(dim = c(nrow(arr.clim_min), ncol(arr.clim_min), 3650))
  
  ## loop through each element (unique pairs of row x column)
  row =  1
  while(row < nrow(arr.anom_max) + 1) {
    col = 1
    while(col < ncol(arr.anom_max) + 1) {
      ## retrieve climatology and anomaly in the cell 
      anom_max <- arr.anom_max[row, col, ]
      anom_min <- arr.anom_min[row, col, ]
      
      year <- 1
      this_index <- leap_years
      ## for leap years, remove element 60 within that year in anomaly data 
      while(year < nrow(leap_years) + 1) {
        if (this_index$leap_year[year] == TRUE) {
          index <- this_index$index[year]
          anom_max <- anom_max[-index+60]
          anom_min <- anom_min[-index+60]
          
          this_index$index <- this_index$index - 1
        }
        year = year + 1
      }
      ## extend climatology and add anomaly to climatology to get temperature 
      clim_max <- rep(arr.clim_max[row, col, ], times = 10)
      temps_max <- clim_max + anom_max 
      clim_min <- rep(arr.clim_min[row, col, ], times = 10)
      temps_min <- clim_min + anom_min 
      
      ## store temperatures for that cell 
      arr.temps_max[row, col, ] <- temps_max
      arr.temps_min[row, col, ] <- temps_min
      
      ## calculate average temp in each cell for each day over the period and store 
      x = 1
      iteration = (rep-1950)/10
      while (x < 366) {
        mean_max[row, col, x+365*iteration] <- mean(arr.temps_max[row, col, seq(x, 3650, 365)], 
                                                    na.rm = TRUE)
        mean_min[row, col, x+365*iteration] <- mean(arr.temps_min[row, col, seq(x, 3650, 365)], 
                                                    na.rm = TRUE)
        x = x +1
      }
      col = col + 1
    }
    row = row + 1
  }
  rep = rep + 10
}


mean_max_new <- array(dim = c(360, 180, 365))
mean_min_new <- array(dim = c(360, 180, 365))
final_max <- data.frame(matrix(ncol = 3))
colnames(final_max) = c("latitude", "longitude", "seasonal_high_temp")
final_min <- data.frame(matrix(ncol = 3))
colnames(final_min) = c("latitude", "longitude", "seasonal_low_temp")

x = 1
row = 1
while(row < nrow(mean_max) + 1) {
  col = 1
  while (col < ncol(mean_max) +1) {
    x = 1
    while (x < 366) {
      ## get mean min and max temp for each cell over all 70 years:
      mean_max_new[row, col, x] <- mean(mean_max[row, col, seq(x, 3650, 365)], 
                                        na.rm = TRUE)
      mean_min_new[row, col, x] <- mean(mean_min[row, col, seq(x, 3650, 365)], 
                                        na.rm = TRUE)
      x = x+1
    }
    ## get maximum mean daily temperature and minimum mean daily temperature for cell:
    final_max <- rbind(final_max, c(lat[col], long[row], max(mean_max_new[row, col,], na.rm=TRUE)))
    final_min <- rbind(final_min, c(lat[col], long[row],  min(mean_min_new[row, col,], na.rm=TRUE)))
    
    col = col + 1
  }
  row = row + 1
}

final_max <- final_max[-1,] %>%
  select(longitude, latitude, seasonal_high_temp)
final_min <- final_min[-1,] %>%
  select(longitude, latitude, seasonal_low_temp)

## save data:
write.csv(final_max, "data-processed/terrestrial_seasonal-max-temps.csv", row.names = FALSE)
write.csv(final_min, "data-processed/terrestrial_seasonal-min-temps.csv", row.names = FALSE)



####################################################################################
#####                 MARINE SEASONAL TEMPERATURE HIGH AND LOW                ######
####################################################################################
firstset <- "/Volumes/ADATA HV620/temperature-data/marine/sst.wkmean.1981-1989.nc"
secondset <- "/Volumes/ADATA HV620/temperature-data/marine/sst.wkmean.1990-present.nc"
landmask <- "/Volumes/ADATA HV620/temperature-data/marine/lsmask.nc" 

ncfile_first <- nc_open(firstset)
ncfile_second <- nc_open(secondset)
ncfile_mask <- nc_open(landmask)

## create variables for things needed to use data
lat <- ncvar_get(ncfile_first, "lat")
long <- ncvar_get(ncfile_first, "lon")
mask <- ncvar_get(ncfile_mask, "mask")

## close the files
nc_close(ncfile_first)
nc_close(ncfile_second)
nc_close(ncfile_mask)

lt_weekly_means <- array(dim = c(360, 180, 52))

firstset <-
  paste("/Volumes/ADATA HV620/temperature-data/marine/sst.wkmean.1981-1989.nc",
        sep = "")
secondset <-
  paste("/Volumes/ADATA HV620/temperature-data/marine/sst.wkmean.1990-present.nc",
        sep = "")
ncfile_first <- nc_open(firstset)
ncfile_second <- nc_open(secondset)

## week centred on Sunday
## first day is October 29, 1981
## last day is December 28, 1989
time_first <- ncvar_get(ncfile_first, "time")

## week centred on Wednesday
## first day is December 31, 1989
## last day is October 4, 2020
time_second <- ncvar_get(ncfile_second, "time")

## create variables for data
sst_first <- ncvar_get(ncfile_first, "sst")
sst_second <- ncvar_get(ncfile_second, "sst")

## close files
nc_close(ncfile_first)
nc_close(ncfile_second)

max_weekly_mean <- data.frame(matrix(ncol = 3))
colnames(max_weekly_mean) = c("latitude", "longitude", "seasonal_high_temp")
min_weekly_mean <- data.frame(matrix(ncol = 3))
colnames(min_weekly_mean) = c("latitude", "longitude", "seasonal_low_temp")

## loop through each element (unique pairs of row x column)
row =  1
while (row < nrow(sst_first) + 1) {
  col = 1
  while (col < ncol(sst_first) + 1) {
    ## if it is a cell on land, skip it
    if (mask[row, col] == 0) {
      col = col + 1
    }
    else {
      ## retrieve temps in cell
      temps <- sst_first[row, col,]
      temps <- append(temps, sst_second[row, col,])
      
      ## calculate average temp in each cell for each week over the period and store
      x = 1
      while (x < 53) {
        lt_weekly_means[row, col, x] <- mean(temps[seq(x, 2033, 52)],
                                             na.rm = TRUE)
        
        x = x + 1
      }
      ## get maximum weekly mean temperature and minimum weekly mean temperature for cell
      max_weekly_mean <- rbind(max_weekly_mean, c(lat[col], long[row], max(lt_weekly_means[row, col,], na.rm=TRUE)))
      min_weekly_mean <- rbind(min_weekly_mean, c(lat[col], long[row],  min(lt_weekly_means[row, col,], na.rm=TRUE)))
      
      col = col + 1
    }
  }
  row = row + 1
}

max_weekly_mean <- max_weekly_mean[-1,] %>%
  select(longitude, latitude, seasonal_high_temp) %>%
  mutate(longitude = ifelse(longitude < 180, longitude, longitude - 360))

min_weekly_mean <- min_weekly_mean[-1,] %>%
  select(longitude, latitude, seasonal_low_temp) %>%
  mutate(longitude = ifelse(longitude < 180, longitude, longitude - 360))


## save these datasets as seasonal high and low temps
write.csv(max_weekly_mean, "data-processed/marine_seasonal-max-temps.csv", row.names = FALSE)
write.csv(min_weekly_mean, "data-processed/marine_seasonal-min-temps.csv", row.names = FALSE)
