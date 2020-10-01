## collating IUCN range maps
library(tidyverse)
library(rgdal)
library(sf)
library(rnaturalearth)
library(latticeExtra)

thermal_limits <- read_csv("./data-raw/globtherm_full_dataset_2019.csv") %>%
  filter(thermy == "ectotherm")


## figure out which fish polygon groups to download
families <- c("Pomacanthidae", "Blenniidae"	, "Albulidae", "Elopidae", "Megalopidae", "Chaetodontidae", "Epinephelidae", "Tetraodontidae", "Sparidae", "Centracanthidae", "Acanthuridae", "Syngnathidae", "Aulostomidae", "Centriscidae", "Fistulariidae", "Solenostomidae", "Istiophoridae", "Scombridae" , "Xiphiidae", "Labridae", "Scaridae")
classes <- c("Chondrichthyes", "Myxini")
orders <- c("Clupeiformes")

families_overlap <- families[which(families %in% thermal_limits$Family)]
classes_overlap <- classes[which(classes %in% thermal_limits$Class)]
orders_overlap <- orders[which(orders %in% thermal_limits$Order)]

## download sets that include: "Blenniidae", "Tetraodontidae", "Sparidae", "Labridae","Chondrichthyes", "Clupeiformes"



## sort through downloaded range maps, filter out species we do not have in thermal tolerance data
thermal_species <- unique(paste(thermal_limits$Genus, thermal_limits$Species, sep = " "))

amphibs_all <- st_read("/Volumes/ADATA HV620/IUCN/AMPHIBIANS/AMPHIBIANS.shp")
amphibs <- thermal_species[which(thermal_species %in% amphibs_all$binomial)]

amphibs_overlap <- amphibs_all %>%
  filter(binomial %in% amphibs)
rm(amphibs_all,amphibs)

blennies_all <- st_read("/Volumes/ADATA HV620/IUCN/BLENNIES/BLENNIES.shp")
blennies <- thermal_species[which(thermal_species %in% blennies_all$binomial)]

blennies_overlap <- blennies_all %>%
  filter(binomial %in% blennies)
rm(blennies_all,blennies)

clups_all <- st_read("/Volumes/ADATA HV620/IUCN/CLUPEIFORMES/CLUPEIFORMES.shp")
clups <- thermal_species[which(thermal_species %in% clups_all$binomial)]

clups_overlap <- clups_all %>%
  filter(binomial %in% clups)
rm(clups_all,clups)

puff_all <- st_read("/Volumes/ADATA HV620/IUCN/PUFFERFISH/PUFFERFISH.shp")
puff <- thermal_species[which(thermal_species %in% puff_all$binomial)]

puff_overlap <- puff_all %>%
  filter(binomial %in% puff)
rm(puff_all,puff)

reptiles_all <- st_read("/Volumes/ADATA HV620/IUCN/REPTILES/REPTILES.shp")
reptiles <- thermal_species[which(thermal_species %in% reptiles_all$binomial)]

reptiles_overlap <- reptiles_all %>%
  filter(binomial %in% reptiles)
rm(reptiles_all,reptiles)

seabream_all <- st_read("/Volumes/ADATA HV620/IUCN/SEABREAMS_PORGIES_PICARELS/SEABREAMS_PORGIES_PICARELS.shp")
seabream <- thermal_species[which(thermal_species %in% seabream_all$binomial)]

seabream_overlap <- seabream_all %>%
  filter(binomial %in% seabream)
rm(seabream_all,seabream)

sharks_all <- st_read("/Volumes/ADATA HV620/IUCN/SHARKS_RAYS_CHIMAERAS/SHARKS_RAYS_CHIMAERAS.shp")
sharks <- thermal_species[which(thermal_species %in% sharks_all$binomial)]

sharks_overlap <- sharks_all %>%
  filter(binomial %in% sharks)
rm(sharks_all,sharks)

wrasses_all <- st_read("/Volumes/ADATA HV620/IUCN/WRASSES_PARROTFISHES/WRASSES_PARROTFISHES.shp")
wrasses <- thermal_species[which(thermal_species %in% wrasses_all$binomial)]

wrasses_overlap <- wrasses_all %>%
  filter(binomial %in% wrasses)
rm(wrasses_all,wrasses)


combined <- rbind(amphibs_overlap, blennies_overlap, clups_overlap, puff_overlap, reptiles_overlap, seabream_overlap, sharks_overlap, wrasses_overlap) 


## write out to file:
st_write(combined, "/Volumes/ADATA HV620/IUCN/FILTERED/IUCN-ectotherms.shp", driver = "ESRI Shapefile")







##############################################################
## plot one:
one_amphib <- filter(combine, binomial == "Hylarana erythraea")
plot(st_geometry(one_amphib))


crs(one_amphib)

## transform it:

one_amphib <- st_transform(one_amphib, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
plot(st_geometry(one_amphib))
crs(one_amphib)


countries <- ne_countries(returnclass = "sf") 
countries <- st_transform(countries, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")

plot(st_geometry(countries))
plot(st_geometry(one_amphib), add = TRUE, col = "red")

ggplot(data = world) +
  geom_sf()
