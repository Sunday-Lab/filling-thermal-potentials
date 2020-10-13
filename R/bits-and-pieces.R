## code for various purposes that I cannot precisely describe in a succint sentence
library(tidyverse)



###### investigating range duplicates
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




## looking at how complete traits are:
traits <- read_csv("data-raw/globtherm_traits_collated_180617.csv") %>%
  mutate(genus_species = paste(.$Genus, .$Species, sep = "_"))

thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv")

traits <- traits[traits$genus_species %in% thermal_limits$genus_species,]

## fix wonky column names 
colnames(traits) <- str_replace_all(colnames(traits), pattern = " ", replacement = "_") %>%
  str_replace_all(., pattern = "\\(", replacement = "") %>%
  str_replace_all(., pattern = "\\)", replacement = "") %>%
  str_replace_all(., pattern = "\\;", replacement = "") %>%
  str_replace_all(., pattern = "\\/", replacement = "_") %>%
  str_replace_all(., pattern = "\\,", replacement = "") %>%
  str_replace_all(., pattern = "\\.", replacement = "") %>%
  str_replace_all(., pattern = "\\?", replacement = "")
colnames(traits)


## count how many are NA in each column 
na_count <- sapply(traits, function(y) sum(length(which(is.na(y)))))
na <- as.data.frame(na_count) %>%
  mutate(column = row.names(.)) %>%
  filter(!str_detect(column, "source")) %>%
  filter(!str_detect(column, "Source")) %>%
  filter(!str_detect(column, "notes"))

summary <- ggplot(na, aes(x = reorder(column, -na_count), y = na_count, fill = "orange")) + geom_col() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  coord_flip() + 
  labs(y = "Trait", x = "Number of values missing")

ggsave(summary, filename = "figures/traits_completeness-summary.png", device = "png")



## plotting expectations of model results:
df <- data.frame()
ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "|Latitudinal midpoint|", y = "Range filling") +
  scale_x_continuous(breaks = c(0,30,60,90), labels = c("0°","30°", "60°", "90°"), limits = c(0,91))

ggsave(filename = "figures/expectations_latitudinal-midpoint.png", device = "png", height = 2, width = 3, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Habitat type", y = "Range filling") +
scale_x_continuous(breaks = c(0,30,60,90, 120), labels = c("Marine", "Coastal", "Intertidal", "Terrestrial", "Freshwater"), limits = c(0,121))

ggsave(filename = "figures/expectations_habitat-type.png", device = "png", height = 2, width = 3, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Range size (km2)", y = "Range filling") +
  scale_x_continuous(breaks = waiver(), limits = c(0,1000000))

ggsave(filename = "figures/expectations_range-size.png", device = "png", height = 2, width = 3, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Trophic position", y = "Range filling") +
  scale_x_continuous(breaks = c(0,30,60,90, 120), labels = c("Primary \n producer", "Herbivore", "Insectivore", "Omnivore", "Carnivore"), limits = c(0,121))

ggsave(filename = "figures/expectations_trophic-position.png", device = "png", height = 2, width = 3.5, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Body size (cm)", y = "Range filling") +
  scale_x_continuous(breaks = waiver(), limits = c(0,580))

ggsave(filename = "figures/expectations_body-size.png", device = "png", height = 2, width = 3, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Dispersal distance category", y = "Range filling") +
  scale_x_continuous(breaks = c(0,30,60), labels = c("0-1 km", "1-10 km", "10+ km"), limits = c(0,61))

ggsave(filename = "figures/expectations_dispersal-distance.png", device = "png", height = 2, width = 3, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Dispersal type category", y = "Range filling") +
  scale_x_continuous(breaks = c(0,40,80,120,160,200), labels = c("walking", "flying", "non-pelagic \n development \n and sessile adults", "non-pelagic \n development \n and crawling adults", "non-pelagic \n development \n and swimming adults", "pelagic \n  development"), limits = c(0,210))

ggsave(filename = "figures/expectations_dispersal-ability.png", device = "png", height = 3, width = 7, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Warm season dormancy", y = "Range filling") +
  scale_x_continuous(breaks = c(20, 55), labels = c("Y", "N"), limits = c(0,75))

ggsave(filename = "figures/expectations_warm-season-dormancy.png", device = "png", height = 2, width = 3, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Cold season dormancy", y = "Range filling") +
  scale_x_continuous(breaks = c(20, 55), labels = c("Y", "N"), limits = c(0,75))

ggsave(filename = "figures/expectations_cold-season-dormancy.png", device = "png", height = 2, width = 3, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Migratory", y = "Range filling") +
  scale_x_continuous(breaks = c(20, 55), labels = c("Y", "N"), limits = c(0,75))

ggsave(filename = "figures/expectations_migratory.png", device = "png", height = 2, width = 3, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Acclimation response ratio", y = "Range filling") +
  scale_x_continuous(breaks = waiver(), limits = c(0,1.5))

ggsave(filename = "figures/expectations_range-size.png", device = "png", height = 2, width = 3, units = "in")

