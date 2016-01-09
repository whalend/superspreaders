#' # Species Diversity Metrics
#' Now I also want to calculate some plot-level community variables from the vegetation data: diversity and evenness metrics from the vegetation transects
#' 
#'  1. over and understory richness; count of number of species
#'  2. Shannon-Weiner Evenness & Diversity Indices
#'  3. Simpson's Diversity Index
#' 
#' The Shannon and Simpson's indices can be calculated using functions in the `vegan` package. I have done some exploring of this package in the 'vegan_diversity.Rmd' script. 
#' 
#+ load plot sampling data ####

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
#+ understory species coverage data from transects ####

summary(veg_us_2005)# data from vegetation transects 2003, 2004, 2005
veg_us_2005 <- as.tbl(select(veg_us_2005, -X))
summary(veg_us_2011)# data from vegetation transects 2011
veg_us_2011 <- as.tbl(select(veg_us_2011, -X))
veg_us_2011# 1300 records
veg_us_2005# 1249 records

# Each record is a species in a plot, so it was observed in at least one transect sample point.

library(tidyr)

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
us_richness <- left_join(
      veg_sub %>% filter(year == 2005) %>% group_by(plotid) %>% 
            summarise(us_rich05 = length(unique(us_species))),
      veg_sub %>% filter(year == 2011) %>% group_by(plotid) %>% 
            summarise(us_rich14 = length(unique(us_species))))
us_richness


# 2005 vegetation data proportional species coverage
veg_us_2005_pcov <- veg_sub %>% filter(year == 2005) %>%
      ungroup(.) %>% 
      select(plotid, us_species, pcov)
veg_us_2005_pcov <- veg_us_2005_pcov %>% spread(us_species, pcov)
summary(veg_us_2005_pcov)
veg_us_2005_pcov[is.na(veg_us_2005_pcov)] <- 0
# 2005 vegetation data point count species coverage, proxy for abundance
veg_us_2005_ptct <- veg_sub %>% filter(year == 2005) %>% 
      ungroup(.) %>% 
      select(plotid, us_species, ptct)
veg_us_2005_ptct <- veg_us_2005_ptct %>% spread(us_species, ptct)
summary(veg_us_2005_ptct)
veg_us_2005_ptct[is.na(veg_us_2005_ptct)] <- 0

# 2011 vegetation data proportianal species coverage
veg_us_2011_pcov <- veg_sub %>% filter(year == 2011) %>% 
      ungroup(.) %>% 
      select(plotid, us_species, pcov)
veg_us_2011_pcov <- veg_us_2011_pcov %>% spread(us_species, pcov)
summary(veg_us_2011_pcov)
veg_us_2011_pcov[is.na(veg_us_2011_pcov)] <- 0
# 2011 vegetation data point count species coverage, proxy for abundance
veg_us_2011_ptct <-  veg_sub %>% filter(year == 2011) %>% 
      ungroup(.) %>% 
      select(plotid, us_species, ptct)
veg_us_2011_ptct <- veg_us_2011_ptct %>% spread(us_species, ptct)
summary(veg_us_2011_ptct)
veg_us_2011_ptct[is.na(veg_us_2011_ptct)] <- 0



#' ## Calculating Un/tagged & Diversity
#' Next, overstory richness from the DBH data for tagged & untagged stems.
#+ calculate overstory richness ####
summary(untagged_dbh_2005)
summary(untagged_dbh_2014)

untagged_dbh_2005 <- as.tbl(untagged_dbh_2005 %>% rename(plotid = PlotID, species = Species))
untagged_dbh_2014 <- as.tbl(untagged_dbh_2014 %>% rename(plotid = PlotID, species = Species))
untagged_dbh_2005$plotid <- tolower(untagged_dbh_2005$plotid)
untagged_dbh_2014$plotid <- tolower(untagged_dbh_2014$plotid)

unique(untagged_dbh_2005$species)
unique(untagged_dbh_2014$species)
untagged_dbh_2005$species <- as.character(untagged_dbh_2005$species)
untagged_dbh_2014$species <- as.character(untagged_dbh_2014$species)
untagged_dbh_2005[untagged_dbh_2005=="quercus sp."] <- "qusp"
untagged_dbh_2005[untagged_dbh_2005=="arcto sp."] <- "arcto"
untagged_dbh_2005[untagged_dbh_2005=="prunus sp."] <- "prunus"

untagged_dbh_2014[untagged_dbh_2014=="arcto sp."] <- "arcto"
untagged_dbh_2014[untagged_dbh_2014=="quercus sp."] <- "qusp"
untagged_dbh_2014[untagged_dbh_2014=="unknown sp. 1"] <- "unk1"
untagged_dbh_2014[untagged_dbh_2014=="prunus sp."] <- "prunus"

unique(untagged_dbh_2005$plotid)
unique(untagged_dbh_2014$plotid)
untagged_dbh_2005 <- filter(untagged_dbh_2005, plotid != "bush01")
untagged_dbh_2014 <- filter(untagged_dbh_2014, plotid != "bush01")

#' ## Examine live and dead species abundances for 2005 & 2014
#+ untagged abundances ####
untagged_dbh_2005 %>% select(plotid, species, live_count)
alt05live <- select(alt05live, -X)
dim(alt05live)
alt05dead <- select(alt05dead, -X)

untagged_dbh_2014 %>% select(plotid, species, live_count)
alt14live <- select(alt14live, -X)
dim(alt14live)
alt14dead <- select(alt14dead, -X)


# This code provides the number of unique species from the untagged tree data, i.e. the counted untagged species richness in each plot
untagged_dbh_2005
untagged_dbh_2014
untagged_dbh_2005 %>% group_by(plotid, year) %>% summarise(untagged_rich = length(unique(species)))
untagged_dbh_2014 %>% group_by(plotid, year) %>% summarise(untagged_rich = length(unique(species)))


# Calculate the richness of the tagged stems
stems <- as.tbl(stems)
summary(stems)
unique(stems$status)
stems <- stems %>% filter(status == "Alive" | status == "Dead")# filter out missing/removed stems
stems <- droplevels(stems)
stems
stems <- select(stems, -X) %>% rename(plotid = plot)
stems$plotid <- tolower(stems$plotid)
unique(stems$plotid)
unique(stems$species)
stems$species <- tolower(stems$species)
stems$species[stems$species == "quke x quag"] <- "quagXquke"
stems$species[stems$species == "unknown sp."] <- "unk1"
stems$species[stems$species == "quag x quwi"] <- "quagXquwi"
stems <- droplevels(stems)
stems


live_tagged <- stems %>% filter(status == "Alive") %>% 
      group_by(plotid, year, species) %>% 
      summarise(live_count = length(species))
unique(live_tagged$species)
live_tagged

filter(live_tagged, year == 2005, plotid == "ganay02")
filter(untagged_dbh_2005, plotid == "ganay02")
filter(live_tagged, year == 2014, plotid == "ganay02")
filter(untagged_dbh_2014, plotid == "ganay02")

tagged_abund_2005 <- live_tagged %>% filter(year == 2005)
tagged_abund_2005 <- ungroup(tagged_abund_2005) %>% select(-year)
length(unique(tagged_abund_2005$plotid))

tagged_abund_2014 <- live_tagged %>% filter(year == 2014)
tagged_abund_2014 <- ungroup(tagged_abund_2014) %>% select(-year)
length(unique(tagged_abund_2014$plotid))

# Join tagged and untagged stems for 2005
tagged_abund_2005; unique(tagged_abund_2005$species)
untagged_dbh_2005; unique(untagged_dbh_2005$species)

## 2005 data
abund_2005 <- rbind(tagged_abund_2005, select(untagged_dbh_2005, plotid, year, species, live_count))
unique(abund_2005$species)# 26 species
alt05live
abund_2005_wide <- ungroup(abund_2005) %>% 
      select(-year) %>% spread(species, live_count)

## 2014 data
abund_2014 <- rbind(tagged_abund_2014, select(untagged_dbh_2014, plotid, species, live_count))
unique(abund_2014$species)# 24 species
as.tbl(alt14live)
abund_2014_wide <- abund_2014 %>% spread(species, live_count)

# Summarize into overstory species richness - number of each species in each plot for each survey
os_rich2005 <- abund_2005 %>% summarise(os_rich05 = length(unique(species)))
os_rich2014 <- abund_2014 %>% group_by(plotid) %>% 
      summarise(os_rich14 = length(unique(species)))
os_richness <- left_join(ungroup(os_rich2005) %>% select(-year), os_rich2014)

# Join overstory richness and understory richness (counts) into single data frame
richness <- left_join(os_richness, us_richness, by = c("plotid"))
richness
summary(richness)


#' ## Calculate Diversity Indices
#' Calculate diversity metrics for overstory and understory species using tools from the `vegan` package.
#+ vegan diversity ####
library(vegan)

# Overstory Diversity
summary(abund_2005_wide)
abund_2005_wide[is.na(abund_2005_wide)] <- 0
summary(abund_2014_wide)
abund_2014_wide[is.na(abund_2014_wide)] <- 0

H.2005 <- diversity(select(abund_2005_wide, -plotid), index = "shannon", MARGIN = 1, base = exp(1))
D.2005 <- diversity(select(abund_2005_wide, -plotid), index = "simpson")
inv.D.2005 <- diversity(select(abund_2005_wide, -plotid), index = "invsimpson")

H.2014 <- diversity(select(abund_2014_wide, -plotid), index = "shannon", MARGIN = 1, base = exp(1))
D.2014 <- diversity(select(abund_2014_wide, -plotid), index = "simpson")
inv.D.2014 <- diversity(select(abund_2014_wide, -plotid), index = "invsimpson")

diversity.2005 <- select(abund_2005_wide, plotid) %>% 
      mutate(H.2005 = H.2005, D.2005 = D.2005, inv.D.2005 = inv.D.2005)
diversity.2014 <- select(abund_2014_wide, plotid) %>% 
      mutate(H.2014 = H.2014, D.2014 = D.2014, inv.D.2014 = inv.D.2014)
os.diversity <- left_join(diversity.2005, diversity.2014, by = "plotid")
summary(os.diversity)

# Add Pielou's evenness (J) to the overstory diversity dataframe
# J = H/log(species number)
os.diversity$J.2005 <- os.diversity$H.2005/log1p(specnumber(select(abund_2005_wide, -plotid)))
os.diversity$J.2014 <- os.diversity$H.2014/log1p(specnumber(select(abund_2014_wide, -plotid)))


# Understory Diversity - based on transect point counts
H.2005 <- diversity(select(veg_us_2005_ptct, -plotid), index = "shannon", MARGIN = 1, base = exp(1))
D.2005 <- diversity(select(veg_us_2005_ptct, -plotid), index = "simpson")
inv.D.2005 <- diversity(select(veg_us_2005_ptct, -plotid), index = "invsimpson")
H.2011 <- diversity(select(veg_us_2011_ptct, -plotid), index = "shannon", MARGIN = 1, base = exp(1))
D.2011 <- diversity(select(veg_us_2011_ptct, -plotid), index = "simpson")
inv.D.2011 <- diversity(select(veg_us_2011_ptct, -plotid), index = "invsimpson")

diversity.2005 <- select(veg_us_2005_ptct, plotid) %>% 
      mutate(us.H.2005 = H.2005, us.D.2005 = D.2005, us.inv.D.2005 = inv.D.2005)
diversity.2011 <- select(veg_us_2011_ptct, plotid) %>% 
      mutate(us.H.2011 = H.2011, us.D.2011 = D.2011, us.inv.D.2011 = inv.D.2011)
us.diversity <- left_join(diversity.2005, diversity.2011, by = "plotid")
summary(us.diversity)

# Add Pielou's evenness to the understory diversity data frame
us.diversity$us.J.2005 <- us.diversity$us.H.2005/log1p(specnumber(select(veg_us_2005_ptct, -plotid)))
us.diversity$us.J.2011 <- us.diversity$us.H.2011/log1p(specnumber(select(veg_us_2011_ptct, -plotid)))

os.us.diversity <- left_join(os.diversity, us.diversity, by = "plotid")
os.us.diversity <- left_join(os.us.diversity, richness, by = "plotid")

write.csv(os.us.diversity, "analysis/data/os_us_diversity.csv", row.names = F)

#'  - H is the Shannon-Weiner Index
#'  - D is Simpson's Index
#'  - inv.D is the inverse Simpson's Index
#'  - J is Pielou's Evenness
#' 
#' 
#' ## Calculate Species Abundance Models
#' Fisher's log-series estimates the expected number of species with *n* individuals, plotting species by frequencies.
#+ Fishers log series ####
k.2005 <- sample(nrow(select(abund_2005_wide, -plotid)), 1)
fish.2005 <- fisherfit(as.data.frame(select(abund_2005_wide, -plotid))[k.2005,])
plot(fish.2005)
fish.2005

k.2014 <- sample(nrow(select(abund_2014_wide, -plotid)), 1)
fish.2014 <- fisherfit(as.data.frame(select(abund_2014_wide, -plotid))[k.2014,])
plot(fish.2014)
fish.2014
#'
#'
#' ## Preston's Log-normal Model
#' Preston's log-normal model bins species into frequency classes of increasing sizes instead of plotting species by frequencies. There are two alternative fucntions in the `vegan` pacakge for fitting the lognormal model:
#' 
#'  1. `prestonfit` uses binning approach with arbitrary choices of limits and ties
#'    - for ties between adjacent octaves half were in 1st octave and half transferred to next octave
#'    - can either split ties or keep all limit cases in lower octave
#'  
#'  2. `prestondistr` maximizes truncated log-normal likelihood without binning data
#'    - recommended in `vegan` vignette paper
#'    
#+ Prestons log-normal model ####
prest.2005 <- prestondistr(as.data.frame(select(abund_2005_wide, -plotid))[k.2005,])
plot(prest.2005)
prest.2014 <- prestondistr(as.data.frame(select(abund_2014_wide, -plotid))[k.2014,])
plot(prest.2014)
#' This uses maximum-likelihood method and, I think, log base 2 (?)
#' 
#' 
#' ## Whittaker Plots - Ranked Abundance Distribution
#' The `vegan` function `radfit` uses maximum likelihood to fit popular models:
#' 
#'       - brokenstick: null model with no estimated parameters in `vegan`
#'       - preemption: one estimated parameter ($\alpha$)
#'       - log-normal: two estimated parameters ($\mu, sigma$)
#'       - Zipf: two estimated parameters (p1, $\gamma$)
#'       - Zipf-Mandlebrot: three estimated parameters (*c*, $\beta, gamma$)
#'       
#' The models are customarily defined for proportions, but `radfit` works directly with abundance data.
#+ Whittaker plots

# rad <- radfit(BCI[k,], family = poisson)
# plot(rad)
# radlattice(rad)
# rad

rad.2005 <- radfit(as.data.frame(select(abund_2005_wide, -plotid))[k.2005,], family = poisson)
plot(rad.2005)
radlattice(rad.2005)

rad.2014 <- radfit(as.data.frame(select(abund_2014_wide, -plotid))[k.2014,], family = poisson)
plot(rad.2014)
radlattice(rad.2014)

#' This function compares models using AIC or BIC (based on log-likelihood, penalized by number of estimated parameters).

