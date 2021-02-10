## collating IUCN range maps
library(tidyverse)
library(rgdal)
library(sf)
library(rnaturalearth)
library(latticeExtra)

thermal_limits <- read_csv("data-raw/globtherm_full_dataset_2019.csv") %>%
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

seabream_all <- st_read("/Volumes/ADATA HV620/IUCN/SEABREAMS_PORGIES_PICARELS/SEABREAMS_fcrsPORGIES_PICARELS.shp")
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


## combine and filter out non-resident ranges:
combined <- rbind(amphibs_overlap, blennies_overlap, clups_overlap, puff_overlap, 
                   reptiles_overlap, seabream_overlap, sharks_overlap, wrasses_overlap) 
#%>%
  # filter(legend == "Extant (resident)")

## collect same ID number (species) into one MULTIPOLYGON:
combined <- aggregate(combined, list(combined$id_no), function(x) x[1])


## write out to file:
st_write(combined, "/Volumes/ADATA HV620/IUCN/FILTERED/IUCN-ectotherms.shp", driver = "ESRI Shapefile")

##############################################################
## plot one:
one_amphib <- filter(combined, binomial == "Eurycea bislineata")
plot(st_geometry(one_amphib))

crs(one_amphib)

## transform it to the correct projection:
one_amphib <- st_transform(one_amphib, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
plot(st_geometry(one_amphib))
crs(one_amphib)

countries <- ne_countries(returnclass = "sf") 
countries <- st_transform(countries, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")

plot(st_geometry(countries))
plot(st_geometry(one_amphib), add = TRUE, col = "red")





#########################################################################################
## figuring out how many range maps we have (IUCN + GBIF) and how many cross the equator 
IUCN <- st_read("/Volumes/ADATA HV620/IUCN/FILTERED/IUCN-ectotherms.shp") %>%
  select(binomial, geometry) %>%
  rename(species = binomial) 
GBIF <- st_read("/Volumes/ADATA HV620/polygons/Filtered occurences ectotherm animals_020817.shp")
GBIF <- GBIF[GBIF$species %in% thermal_species, ] ## get rid of species not in thermal ectotherm data

st_crs(IUCN)
st_crs(GBIF)

## combine:
realized_ranges <- rbind(IUCN, GBIF) ## okay, we have 524 ectotherm ranges 
length(unique(realized_ranges$species)) ## 439 unique species 

length(unique(IUCN$species)) ## 214
length(unique(GBIF$species)) ## 310

## create sf that represents the equator (a line)
equator <- st_linestring(rbind(c(-180, 0), c(180, 0)))
plot(st_geometry(equator))

## test equator crossing check:
one_species <- filter(IUCN, species == "Hylarana erythraea")
equator_check <- st_intersects(one_species, equator, sparse = FALSE)[,] ## should return TRUE
another_species <- filter(IUCN, species == "Eurycea bislineata")
equator_check <- st_intersects(another_species, equator, sparse = FALSE)[,] ## should return FALSE

both_species <- rbind(one_species, another_species) 
equator_check <- st_intersects(both_species, equator, sparse = FALSE)[,] ## should return TRUE, FALSE

## yay! it works



## see which ranges intersect the equator:
equator_check <- st_intersects(realized_ranges, equator, sparse = FALSE)[,]
crosses_equator <- filter(realized_ranges, equator_check == TRUE)

## see which remaining ranges have Northern and Southern parts
does_not <- filter(realized_ranges, equator_check == FALSE)
pts <- matrix(c(-180,-90,-180,90,180,90,180,-90,-180,-90),ncol=2, byrow=TRUE)
pts_n <- matrix(c(-180,0,-180,90,180,90,180,0,-180,0),ncol=2, byrow=TRUE)
pts_s <- matrix(c(-180,-90,-180,0,180,0,180,-90,-180,-90),ncol=2, byrow=TRUE)
both <- st_polygon(list(pts))
n_hemi <- st_polygon(list(pts_n))
s_hemi <- st_polygon(list(pts_s))
plot(both)
plot(n_hemi, add = TRUE)
plot(s_hemi, add = TRUE, col = "red")

in_north <- st_intersects(does_not, n_hemi, sparse = FALSE)[,]
in_south <- st_intersects(does_not, s_hemi, sparse = FALSE)[,]

in_both <- filter(does_not, in_north == TRUE & in_south == TRUE)
## only one in the north and south 

plot(st_geometry(countries))
plot(st_geometry(equator), add = TRUE)
plot(st_geometry(in_both)[1], add = TRUE, col = "red") 


## write out pngs of the ranges that cross the equator
i = 1
while (i < length(crosses_equator$species) + 1) {
  
  graphics.off()
  plot(st_geometry(countries), main = crosses_equator$species[i])
  plot(st_geometry(crosses_equator)[i], add = TRUE, col = "red")
  plot(st_geometry(equator), add = TRUE, col = "blue")
  
  dev.copy(png, filename=paste("data-processed/equator-crossers/", crosses_equator$species[i], ".png", sep = ""), width = 1000, height = 600);
  dev.off ();
  
  i = i + 1
}




## write new thermal limits database with only species that we have ranges for 
thermal_limits_new <- thermal_limits[paste(thermal_limits$Genus, thermal_limits$Species, sep = " ") %in% realized_ranges$species,]

length(unique(thermal_limits_new$genus_species))

write.csv(thermal_limits_new, "data-processed/thermal-limits_ectotherms-with-ranges.csv", row.names = FALSE)






## look at species with ranges in IUCN and made by Greta
## all 85 are squamata 
thermal_limits_new <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv")
duplicated <- IUCN$species[which(IUCN$species %in% GBIF$species)]

duplicate <- filter(realized_ranges, species == as.character(duplicated[1]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[1]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
plot(st_geometry(countries), add = TRUE)
legend(st_bbox(duplicate)$xmin-5, st_bbox(duplicate)$ymin+2, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)

duplicate <- filter(realized_ranges, species == as.character(duplicated[2]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[2]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
plot(st_geometry(countries), add = TRUE)
legend(st_bbox(duplicate)$xmin-5, st_bbox(duplicate)$ymin+10, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)

duplicate <- filter(realized_ranges, species == as.character(duplicated[3]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[3]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
plot(st_geometry(countries), add = TRUE)
legend(st_bbox(duplicate)$xmin+8, st_bbox(duplicate)$ymin+9, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)

duplicate <- filter(realized_ranges, species == as.character(duplicated[4]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[4]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
plot(st_geometry(countries), add = TRUE)
legend(st_bbox(duplicate)$xmin+28, st_bbox(duplicate)$ymin+22, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)


## see how different IUCN ranges are from Greta's polygons:
i = 1 
props <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(props) <- c("species", "IUCN_filled", "overlap", "GBIF_filled")

while(i < length(duplicated) + 1) {
  duplicate <- filter(realized_ranges, species == as.character(duplicated[i]))
  IUCNrange <- as_Spatial(st_geometry(duplicate)[1])
  GBIFrange <- as_Spatial(st_geometry(duplicate)[2])
  
  overlap <- intersect(GBIFrange, IUCNrange)
  overlap <- area(overlap) / 10000000
  
  leftover_IUCN <- (area(IUCNrange)/10000000) - overlap 
  leftover_GBIF <-  (area(GBIFrange)/10000000) - overlap 
  row <- data.frame(species = as.character(duplicated[i]), 
                    IUCN_filled = leftover_IUCN, overlap = overlap, GBIF_filled = leftover_GBIF)
  props <- rbind(props, row)
  
  i = i + 1
}

props <- gather(props, "leftover_type", "area", -species)
props$leftover_type <- factor(props$leftover_type, levels = c("GBIF_filled", "overlap", "IUCN_filled"))

ggplot(props, aes(fill = leftover_type, y = area, x=reorder(species, -area))) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  labs(y = "Area (km^2)", x = "Species", fill = "") +
  scale_fill_manual(values =  c("darkgoldenrod1", "azure4", "darkorange3"), 
                    labels=c("GBIF range not in IUCN range",
                             "Overlap between ranges",
                             "IUCN range not in GBIF range")) +
  theme(axis.text.y = element_text(size = 6))

ggsave(device = "png", filename = "figures/IUCN-GBIF-range-overlap.png", height = 6, width = 10)

## make plot to explain what overlaps mean:
plot(IUCNrange, col = "darkorange3")
plot(GBIFrange, add = TRUE, col = "darkgoldenrod1")
plot(overlap, add = TRUE, col = "azure4")
plot(st_geometry(countries), add = TRUE)

dev.copy(png, filename = "figures/IUCN-GBIF-range-overlap_example.png", width = 1000, height = 600);
dev.off()
