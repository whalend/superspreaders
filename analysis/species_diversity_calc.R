#' # Species Diversity Metrics
#' Now I also want to calculate some plot-level community variables from the vegetation data: diversity and evenness metrics from the vegetation transects
#' 
#'  1. over and understory richness; count of number of species
#'  2. Shannon-Weiner Evenness & Diversity Indices
#'  3. Simpson's Diversity Index
#' 
#' The Shannon and Simpson's indices can be calculated using functions in the `vegan` package. I have done some exploring of this package in the 'vegan_diversity.Rmd' script. 
#' 
#+ load plot sampling data

## Ecological data: leaf counts, infection, at stem and plot levels
load("~/GitHub/superspreaders/stems_plots.RData")

## Vegetation transect data
veg_us_2005 <- read.csv("analysis/data/understory_veg_2005.csv")
veg_us_2011 <- read.csv("analysis/data/understory_veg_2011.csv")
canopy_cover <- read.csv("analysis/data/plot_canopy_cover_dens.csv")

## Untagged dbh data
untagged_dbh_2005 <- read.csv("analysis/data/untagged_dbh2005_summary.csv")
untagged_dbh_2014 <- read.csv("analysis/data/untagged_dbh2014_summary.csv")
alt05live <- read.csv("analysis/data/altsp05abund_live.csv")
alt05dead <- read.csv("analysis/data/altsp05abund_dead.csv")
alt14live <- read.csv("analysis/data/altsp14abund_live.csv")
alt14dead <- read.csv("analysis/data/altsp14abund_dead.csv")
#'
#'
#' Calculate understory richness
#+ understory species coverage data from transects

summary(veg_us_2005)# data from vegetation transects 2003, 2004, 2005
veg_us_2005 <- as.tbl(select(veg_us_2005, -X))
summary(veg_us_2011)# data from vegetation transects 2011
veg_us_2011 <- as.tbl(select(veg_us_2011, -X))
veg_us_2011# 1300 records
veg_us_2005# 1249 records

# Each record is a species in a plot, so it was observed in at least one transect sample point.

#' Transpose vegetation transect data from long to wide format. The wide format preps the data for the `diversity` function in the `vegan` package.
#' 
#+ transpose data long to wide
library(tidyr)
# 2005 vegetation data proportional species coverage
veg_us_2005_wide1 <- select(veg_us_2005, PlotID, UnderstorySp, pcov)
veg_us_2005_wide1 <- veg_us_2005_wide1 %>% spread(UnderstorySp, pcov)
summary(veg_us_2005_wide1)
# 2005 vegetation data point count species coverage
veg_us_2005_wide2 <- select(veg_us_2005, PlotID, UnderstorySp, ptct)
veg_us_2005_wide2 <- veg_us_2005_wide2 %>% spread(UnderstorySp, ptct)
summary(veg_us_2005_wide2)

# 2011 vegetation data proportianal species coverage
veg_us_2011_wide1 <- select(veg_us_2011, PlotID, UnderstorySp, pcov)
veg_us_2011_wide1 <- veg_us_2011_wide1 %>% spread(UnderstorySp, pcov)
summary(veg_us_2011_wide1)
# 2011 vegetation data point count species coverage
veg_us_2011_wide2 <- select(veg_us_2011, PlotID, UnderstorySp, ptct)
veg_us_2011_wide2 <- veg_us_2011_wide2 %>% spread(UnderstorySp, ptct)
summary(veg_us_2011_wide2)


veg_us <- rbind(veg_us_2005, veg_us_2011)
unique(veg_us$year)
veg_us$year[veg_us$year < 2011] <- 2005 # Creates the two most complete veg sampling
veg_us <- veg_us %>% rename(us_species = UnderstorySp, plotid = PlotID)
unique(veg_us$us_species)

# Change species ID to be continuous (no spaces) for plants identified to genus
veg_us$us_species <- as.character(veg_us$us_species)
veg_us[veg_us=="quercus sp."] <- "qusp"
veg_us[veg_us=="unknown sp. 1"] <- "unk1"
veg_us[veg_us=="vaccinium sp."] <- "vacc"
veg_us[veg_us=="vinca sp."] <- "vinca"
veg_us[veg_us=="prunus sp."] <- "prunus"
veg_us[veg_us=="rubus sp."] <- "rubus"
veg_us[veg_us=="ribes sp."] <- "ribes"
veg_us[veg_us=="ceanothus sp."] <- "ceanothus"
veg_us[veg_us=="lupine sp."] <- "lupine"
veg_us[veg_us=="mimulus sp."] <- "mimulus"
veg_us[veg_us=="mimulus sp"] <- "mimulus"
veg_us[veg_us=="unknown sp. 2"] <- "unk2"
veg_us[veg_us=="garrya sp."] <- "garrya"
veg_us[veg_us=="rose sp."] <- "rosa"
veg_us[veg_us=="quag x quke"] <- "quagXquke"
veg_us[veg_us=="quag x quga"] <- "quagXquke"
veg_us[veg_us=="rhamnus sp."] <- "rhamnus"
veg_us[veg_us=="salix sp."] <- "salix"
veg_us[veg_us=="cotoneaster sp"] <- "cotoneaster"
veg_us[veg_us=="quke x quwi"] <- "qukeXquwi"
veg_us[veg_us=="quke x quch"] <- "qukeXquch"
veg_us[veg_us=="rosa sp."] <- "rosa"
sort(unique(veg_us$us_species))

# I want to match the species to those used by Sarah in her second chapter
local_woody <- read.csv("analysis/data/presence_local_woody_n197.csv")
summary(local_woody)
unique(local_woody$spp.local)
woody_spp <- as.character(local_woody$spp.local)

sort(woody_spp) # only 29 species
sort(unique(veg_us$us_species))# 79 species
# It looks like I will want to add a few more to her list

# Add other selected woody species to the list for inclusion in species richness/diversity/evenness calculations
woody_spp <- c(woody_spp, "ceanothus", "cebe", "cecu", "frla", "garrya", "gemo", "juca", "mafu", "mane", "mimulus", "pisa", "quagXquke", "qube", "qusp", "quwi", "rhamnus", "rhcr", "ribes", "rica", "risa", "roca", "rosa", "rubus", "rudi", "salix", "syal", "vacc")
veg_sub <- filter(veg_us, veg_us$us_species %in% woody_spp)# 56 species
veg_sub <- droplevels(veg_sub)# drop the unused levels from veg species

# Calculate basic richness as number of different species in each plot each year
richness <- veg_sub %>%
      group_by(plotid, year) %>% 
      summarise(us_rich = length(unique(us_species)))
richness

#' # Calculating Un/tagged & Diversity
#' Next, overstory richness from the DBH data for tagged & untagged stems.

summary(untagged_dbh)
untagged_dbh <- as.tbl(untagged_dbh %>% select(-X) %>% rename(plotid = PlotID, species = Species))
untagged_dbh$plotid <- tolower(untagged_dbh$plotid)
untagged_dbh
unique(untagged_dbh$year)
filter(untagged_dbh, plotid == "ganay02")

untagged_dbh$year[untagged_dbh$year < 2014] <- 2005 # Create two most complete sample groups
unique(untagged_dbh$plotid)
untagged_dbh <- filter(untagged_dbh, plotid != "bush01")
untagged_dbh
unique(untagged_dbh$year)

# Make untagged species abundance data frames for 2005 & 2014
untagged_abund_2005 <- untagged_dbh %>% 
      filter(year == 2005) %>% 
      select(plotid, species, live_count)
untagged_abund_2014 <- untagged_dbh %>% 
      filter(year == 2014) %>% 
      select(plotid, species, live_count)

# This code provides the number of unique species from the untagged tree data
untagged_dbh %>% group_by(plotid, year) %>% summarise(untagged_rich = length(unique(species)))


# Calculate the richness of the tagged stems
stems <- as.tbl(stems)
summary(stems)
unique(stems$status)
stems <- stems %>% filter(status == "Alive" | status == "Dead")# filter out missing stems
stems <- droplevels(stems)
unique(stems$species)
unique(stems$plotid)
stems$species <- tolower(stems$species)
stems$species[stems$species == "quke x quag"] <- "quag x quke"
stems$plotid <- tolower(stems$plotid)
stems <- droplevels(stems)


live_tagged <- stems %>% filter(status == "Alive") %>% group_by(plotid, year, species) %>% summarise(live_count = length(species))
unique(live_tagged$species)
live_tagged

filter(live_tagged, year == 2005, plotid == "ganay02")
filter(untagged_dbh, plotid == "ganay02")
filter(untagged_dbh, year == 2014, plotid == "ganay02")

tagged_abund_2005 <- live_tagged %>% filter(year == 2005) %>% 
      tagged_abund_2005 <- ungroup(tagged_abund_2005) %>% select(-year)
tagged_abund_2014 <- live_tagged %>% filter(year == 2014)
tagged_abund_2014 <- ungroup(tagged_abund_2014) %>% select(-year)

tagged_abund_2005; unique(tagged_abund_2005$species)

untagged_abund_2005; unique(untagged_abund_2005$species)


dbh_abund_2005 <- rbind(tagged_abund_2005, untagged_abund_2005)
dbh_abund_2005[525:530,]
dbh_abund_2005_wide <- 
      dbh_abund_2014 <- rbind(tagged_abund_2014, untagged_abund_2014)
dbh_abund_2014

# Join the tagged and untagged species data
os_abundance <- rbind(live_tagged, select(untagged_dbh, plotid, species, year, live_count))
os_abundance
# This is the number of each species in each plot during each year.
unique(os_abundance$species)# 28 overstory species

# Summarize into overstory species richness - number of each species in each plot for each survey
os_richness <- left_join(
      stems %>% filter(year == 2005 | year == 2014) %>% 
            group_by(plotid, year) %>% 
            summarise(tagged_rich = length(unique(species))),
      untagged_dbh %>% 
            group_by(plotid, year) %>% 
            summarise(untagged_rich = length(unique(species))), by = c("plotid", "year")) %>% 
      mutate(os_rich = rowSums(cbind(tagged_rich,untagged_rich), na.rm = T))
os_richness
# broken down by tagged & untagged richness, and then total richness 'os_rich'

# Join overstory richness and understory richness (counts) into single data frame
richness <- left_join(os_richness, richness, by = c("plotid","year"))
richness
summary(richness)



