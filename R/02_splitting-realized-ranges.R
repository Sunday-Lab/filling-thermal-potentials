## script to split realized ranges into equatorward and poleward components and then rasterize them
library(tidyverse)
library(sp)
library(sf)
library(raster)
library(ncdf4)
library(rnaturalearth)
library(smoothr)
library(rmapshaper)
select <- dplyr::select
## devtools::install_github("r-spatial/lwgeom")
library(lwgeom)


####################################################################################
#####                       SPLITTING RANGE SHAPEFILES                        ######
####################################################################################
## read in thermal limits:
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = " "))

## read in shape files
IUCN <- st_read("/Volumes/ADATA HV620/IUCN/FILTERED/IUCN-ectotherms.shp") %>%
  select(binomial, geometry) %>%
  rename(species = binomial) 
GBIF <- st_read("/Volumes/ADATA HV620/polygons/Filtered occurences ectotherm animals_020817.shp")
GBIF <- GBIF[GBIF$species %in% thermal_limits$genus_species, ] ## get rid of species not in thermal ectotherm data

realized_ranges <- rbind(IUCN, GBIF)

## split into equator-crossers, northern hemisphere and southern hemisphere ranges:
equator <- st_linestring(rbind(c(-180, 0), c(180, 0)))
n_hemi <- st_polygon(list(matrix(c(-180,0,-180,90,180,90,180,0,-180,0),ncol=2, byrow=TRUE)))
s_hemi <- st_polygon(list(matrix(c(-180,-90,-180,0,180,0,180,-90,-180,-90),ncol=2, byrow=TRUE)))

equator_check <- st_intersects(realized_ranges, equator, sparse = FALSE)[,]
does_not <- filter(realized_ranges, equator_check == FALSE)

in_north <- st_intersects(does_not, n_hemi, sparse = FALSE)[,]
in_south <- st_intersects(does_not, s_hemi, sparse = FALSE)[,]

in_both <- filter(does_not, in_north == TRUE & in_south == TRUE)
## only one in the north and south, mostly in the north so consider to be in_north for now

## create sf object for ranges with information about hemisphere:
northern_ranges <- filter(does_not, in_north == TRUE) %>%
  mutate(hemisphere = "N")

southern_ranges <- filter(does_not, in_north == FALSE & in_south == TRUE) %>%
  mutate(hemisphere = "S")

crosses_equator <- filter(realized_ranges, equator_check == TRUE) %>%
  mutate(hemisphere = "EQUATOR")

realized_ranges <- rbind(northern_ranges, southern_ranges, crosses_equator)

st_write(realized_ranges, "data-processed/realized-ranges_unsplit.shp")

## loop through all realized ranges
i = 1
while (i < length(st_geometry(realized_ranges)) + 1) {
  
  print(paste("On species number:", i, sep = " "))
  
  ## get range and latitudinal midpoint
  ## for ranges that are very complex and take forevvvvver to split, simplify first:
  if (i %in% c(100, 173, 190, 191, 192, 196, 197, 198, 199, 205, 206, 209, 225, 229, 426, 524)) {
    range <- ms_simplify(realized_ranges[i, ])
  }
  else {
    range <- realized_ranges[i, ]
  }
  
  extent <- extent(range)
  latitudinal_midpoint <- (as.numeric(ymin(extent)) + as.numeric(ymax(extent))) / 2
  lmp_ls <- st_linestring(rbind(c(-180, latitudinal_midpoint), c(180, latitudinal_midpoint)))
  
  ## for northern and southern species, chop the range at the latitudinal midpoint:
  if (range$hemisphere != "EQUATOR") {
    cut_line <- lmp_ls
    lat_of_cut <- latitudinal_midpoint
    print(paste("Species range is in the", range$hemisphere, "hemisphere.", sep = " "))
  }
  ## for equator crossers, chop the range at the equator:
  else if (range$hemisphere == "EQUATOR") {
    cut_line <- equator
    lat_of_cut <- 0
    print("Species range crosses the equator.")
  }
  
  ## split the range: 
  split <- st_split(range, cut_line)
  
  ## separate into polygons that are below cut vs polygons that are above cut line:
  polys <- list(st_geometry(split))[1][[1]][[1]]
  above <- which(unlist(lapply(polys, 
                               function (x) st_coordinates(st_centroid(x))[2]))
                 > lat_of_cut) %>%
    polys[.]
  
  below <- which(unlist(lapply(polys, function (x) st_coordinates(st_centroid(x))[2])) 
                 < lat_of_cut) %>%
    polys[.]
  
  ## create simple features to house the split pieces of the ranges:
  sf_below <- st_sfc(lapply(below, function (x) st_polygon(x))) %>%
    st_as_sf(.) %>%
    mutate(species = range$species) %>%
    mutate(latitudinal_midpoint = latitudinal_midpoint) %>%
    mutate(poleward_or_equatorward = "equatorward") %>%
    mutate(hemisphere = range$hemisphere)
  
  sf_above <- st_sfc(lapply(above, function (x) st_polygon(x))) %>%
    st_as_sf(.) %>%
    mutate(species = range$species) %>%
    mutate(latitudinal_midpoint = latitudinal_midpoint) %>%
    mutate(poleward_or_equatorward = "poleward") %>%
    mutate(hemisphere = range$hemisphere)
  
  if (i == 1) {
    sf_cumulative <- rbind(sf_above, sf_below)
  }
  else {
    sf_cumulative <- rbind(sf_cumulative, sf_above, sf_below)
  }
  
  i = i + 1
}


##saveRDS(sf_cumulative, "data-processed/sf_cumulative.rds")
rds <- readRDS("data-processed/sf_cumulative.rds")
st_write(sf_cumulative, "data-processed/realized-ranges_split.shp")

## code to visualize each range-splitting action:
countries <- ne_countries(returnclass = "sf") %>%
  st_transform(., "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")

plot(st_geometry(range), col = "lightpink") ## range boundary
plot(st_geometry(polys[[1]]), add = TRUE, col = "lightblue") ## plot polygons separately
plot(st_geometry(polys[[2]]), add = TRUE, col = "lightyellow")
plot(st_geometry(sf_above), add = TRUE, col = "orange2") ## plots all polygons above midpoint
plot(st_geometry(sf_below), add = TRUE, col = "yellow1") ## plots all polygons below midpoint
plot(st_geometry(countries), add = TRUE) 
plot(st_geometry(equator), add = TRUE, col = "red") ## equator
plot(st_geometry(lmp_ls), add = TRUE, col = "orange") ## latitudinal midpoint