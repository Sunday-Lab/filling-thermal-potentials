## for  cleaning and manipulating the traits database 
library(tidyverse)

## read in all traits
traits_all <- read_csv("data-raw/globtherm_traits_collated_180617.csv") %>%
  mutate(genus_species = paste(.$Genus, .$Species, sep = "_"))

thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv")

## check that no species we have ranges for are missing from the trait database:
length(which(!thermal_limits$genus_species %in% traits_all$genus_species))
## 7 rows species missing! Which ones? 
new <- unique(thermal_limits$genus_species[which(!thermal_limits$genus_species %in% traits_all$genus_species)]) ## 4 species
new_infos <- thermal_limits[which(thermal_limits$genus_species %in% new),] 

types <- new_infos %>%
  select(genus_species, type) %>%
  unique() 

types <- aggregate(types$type, list(types$genus_species), paste, collapse = ", ") %>%
  rename("genus_species" = Group.1, "limit_type" = x)
  
new_infos <- new_infos %>%
  select(-type) %>%
  filter(!duplicated(genus_species)) %>%
  left_join(., types, by = "genus_species") 
  
## subset to only ectotherm species for which we have thermal limits and range
traits_sub <- traits_all[traits_all$genus_species %in% thermal_limits$genus_species,]

## add missing species:
new_species <- traits_sub[1:length(new),1:40] 
new_species[1:length(new),] <- NA 
new_species <- new_species %>%
  mutate(genus_species = new_infos$genus_species) %>%
  mutate(Genus = new_infos$Genus, Species = new_infos$Species, Family = new_infos$Family,
         Phylum = new_infos$Phylum, Class = new_infos$Class, 
         Order = new_infos$Order, Realm = new_infos$realm) %>%
  mutate('data gatherer' = "Nikki") 

traits_sub <- rbind(traits_sub, new_species)

## make new column saying whether species has one or both thermal limits 
lim_types <- thermal_limits %>%
  select(genus_species, type) %>%
  unique() 

lim_types <- aggregate(lim_types$type, list(lim_types$genus_species), paste, collapse = ", ") %>%
  rename("genus_species" = Group.1, "limit_type" = x)

traits_sub <- left_join(traits_sub, lim_types, by = "genus_species")


## write out and start filling in the missing ones!!
write.csv(traits_sub, "data-processed/globtherm_traits_collated_180617_ectotherms-with-limits.csv", row.names = FALSE)




###### inspecting dormancy trait completeness #####
## fix wonky column names 
oldnames <- colnames(traits_sub)
colnames(traits_sub) <- str_replace_all(colnames(traits_sub), pattern = " ", replacement = "_") %>%
  str_replace_all(., pattern = "\\(", replacement = "") %>%
  str_replace_all(., pattern = "\\)", replacement = "") %>%
  str_replace_all(., pattern = "\\;", replacement = "") %>%
  str_replace_all(., pattern = "\\/", replacement = "_") %>%
  str_replace_all(., pattern = "\\,", replacement = "") %>%
  str_replace_all(., pattern = "\\.", replacement = "") %>%
  str_replace_all(., pattern = "\\?", replacement = "")
colnames(traits_sub)

dormancy_sub <- traits_sub %>%
  filter(limit_type == "max, min") %>%
  filter(!is.na(cold_season_dormancy) & !is.na(hot_season_dormancy)) %>%
  mutate(cold_season_dormancy = ifelse(str_detect(cold_season_dormancy, "N") | 
                                         str_detect(cold_season_dormancy, "no"), 
                                       "No", ifelse(str_detect(cold_season_dormancy, "Y"), 
                                                    "Yes", cold_season_dormancy)
                                       
                                       )) %>%
  mutate(hot_season_dormancy = ifelse(str_detect(hot_season_dormancy, "N") | 
                                        str_detect(hot_season_dormancy, "no"), 
                                      "No", hot_season_dormancy))
  

unique(dormancy_sub$cold_season_dormancy)
unique(dormancy_sub$hot_season_dormancy)

## check how many are yes and no
length(which(dormancy_sub$cold_season_dormancy == "Yes")) ## 39/141 
length(which(dormancy_sub$hot_season_dormancy == "Yes")) ## 0/141
