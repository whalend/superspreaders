#' # DBH Data
#' 
#' 

## Read in all tagged stem data, examine data structure ####
# This is only the biologically related field data

tagged_stems <- read.csv("analysis/data/TaggedStemQuery.csv")
head(tagged_stems)
summary(tagged_stems)
str(tagged_stems)

untagged_stems <- read.csv("analysis/data/Veg_TblVegDBH_Summary.csv")
head(untagged_stems)
summary(untagged_stems)
str(untagged_stems)


## Create 'Date' formatted date variable and year variable ####
tagged_stems$date <- as.Date(as.character(tagged_stems$Date), format = "%m/%d/%Y")
class(tagged_stems$date)
untagged_stems$date <- as.Date(as.character(untagged_stems$Date), format = "%m/%d/%Y")

library(lubridate)
tagged_stems$year <- year(tagged_stems$date)
summary(tagged_stems)
untagged_stems$year <- year(untagged_stems$date)
summary(untagged_stems)
detach("package:lubridate", unload=TRUE)

## Replace all `-9999` codes with NA ####
tagged_stems[tagged_stems == -9999] <- NA
summary(tagged_stems)


## Rename columns using `rename` funcion from `dplyr` ####
library(plyr); library(dplyr)
tagged_stems <- select(tagged_stems, -Date)
names(tagged_stems) <- tolower(names(tagged_stems))

tagged_stems <- rename(tagged_stems, plot_status = status,
                       cluster = clusterid, species = speciesid,
                       tag = tagnumber, initial_dbh = dbhcm,  
                       stem_status = stemstatus, sod_killed = sod_dead, 
                       umca_slc = sympleafcount, canker = cankerpresent,
                       lide_slc = lideslc, lide_leaf_symp = lidefoliarsymp, 
                       notes = stemsymptomnotes)
str(tagged_stems)

unique(tagged_stems$stem_status)
tagged_stems$plotid <- tolower(tagged_stems$plotid)
tagged_stems$species <- tolower(tagged_stems$species)
tagged_stems$sod_killed <- ifelse(tagged_stems$sod_killed == 0, "No", "Yes")
tagged_stems <- droplevels(filter(tagged_stems, stem_status == "Alive" | stem_status == "Dead"))
summary(tagged_stems)


untagged_stems <- select(untagged_stems, -Date, -TreeID, -StemID)
names(untagged_stems) <- tolower(names(untagged_stems))
untagged_stems <- rename(untagged_stems, dbh = measure, stem_status = dead)
untagged_stems$stem_status[untagged_stems$stem_status==0] <- "Alive"
untagged_stems$stem_status[untagged_stems$stem_status==1] <- "Dead"
untagged_stems$plotid <- tolower(untagged_stems$plotid)
untagged_stems$species <- tolower(untagged_stems$species)
summary(untagged_stems)

## Write out raw revised data sets ####
write.csv(tagged_stems, "analysis/data/tagged_stems.csv", row.names = F)
write.csv(untagged_stems, "analysis/data/untagged_stems.csv", row.names = F)
write.csv(rbind(
      select(tagged_stems, plotid, species, dbh, stem_status, date, year),
      untagged_stems),
      "analysis/data/all_stems.csv", row.names = F) 

## Further exploration of revised stem data ####
# tagged_stems <- read.csv("analysis/data/tagged_stems.csv") # Read in if necessary
str(tagged_stems)
summary(tagged_stems)
unique(tagged_stems$species)# 8 host species listed, two hybrids
# The hybrid identification and QUWI are rare in the study area.
length(unique(filter(tagged_stems, species == "quwi")$tag))# 36 stems, 1 canker observation during the study period
length(unique(filter(tagged_stems, species == "quag x quwi")$tag))# 11 stems
length(unique(filter(tagged_stems, species == "quke x quag")$tag))# 2 stems

## Commented out because this script starts with corrected data.
# tagged_stems <- read.csv("analysis/data/stem_summary_qry2_04-14.csv")
# summary(tagged_stems)
# 
# tagged_stems$Date <- as.Date(tagged_stems$Date, format = "%m/%d/%Y")
# class(tagged_stems$Date)
# library(lubridate)
# tagged_stems$year <- year(as.Date(tagged_stems$Date, format = "%m/%d/%Y"))
# summary(tagged_stems)
# tagged_stems$year[is.na(tagged_stems$year)] <- 2012
# tagged_stems[tagged_stems == -9999] <- NA
# summary(tagged_stems)


## Subset DBH data to create a data frame of remeasured tagged_stems ####

dbh2003 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2003)
# 11 observations


dbh2004 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2004)
# 2991 observations
anti_join(dbh2003, dbh2004, by = "tag")
# Tagged stem DBH measurements are unique between 2004 and 2005
dbh_a <- union(dbh2004, dbh2003)


dbh2005 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2005)
# 681 observations
anti_join(dbh2005, dbh_a, by = "tag")# this indicates that 2 stems were measured again in 2005
dbh_a <- union(dbh_a, anti_join(dbh2005, dbh_a, by = "tag"))
# I think that this is the complete set of stems that makes up the "plot-establishment" data set. 
dbh_establishment <- dbh_a
write.csv(dbh_establishment, "analysis/data/tagged-stems-2003_2005.csv")


dbh2006 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2006)
anti_join(dbh2006, dbh_a, by = "tag")# 2003, 2004, 2005, 2006 are independent sets of observations, 69 new stems tagged in 2006
dbh_a <- union(dbh_a, anti_join(dbh2006, dbh_a, by = "tag"))


dbh2007 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2007)
anti_join(dbh2007, dbh_a, by = "tag")# there is one tag with remeasured dbh in 2007, 34 new stems tagged in 2007
dbh_a <- union(dbh_a, anti_join(dbh2007, dbh_a, by = "tag"))
anyDuplicated(dbh_a$tag) # check for repeated tags


dbh2008 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2008)
anti_join(dbh2008, dbh_a, by = "tag")# 12 new stems in 2008
dbh_a <- union(dbh_a, anti_join(dbh2008, dbh_a, by = "tag"))


dbh2009 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2009)
anti_join(dbh2009, dbh_a, by = "tag")# 20 new tagged stems in 2009
dbh_a <- union(dbh_a, anti_join(dbh2009, dbh_a, by = "tag"))


dbh2010 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2010)
anti_join(dbh2010, dbh_a, by = "tag")# 101 new stem entering the study
dbh_a <- union(dbh_a, anti_join(dbh2010, dbh_a, by = "tag"))


dbh2011 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2011)
anti_join(dbh2011, dbh_a, by = "tag")# 74 new stems entering the study
dbh_a <- union(dbh_a, anti_join(dbh2011, dbh_a, by = "tag"))


dbh2012 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2012)
anti_join(dbh2012, dbh_a, by = "tag")# remeasured 3839 of 3991 stems that were in the study at some point (may have been in abandoned/decommissioned plots) and tagged 284 new stems due to recruitment or plot quality control.
dbh_a <- union(dbh_a, anti_join(dbh2012, dbh_a, by = "tag"))


dbh2014 <- tagged_stems %>% 
      select(plotid, plot_status, species:dbh, stem_status, year) %>% 
      filter(dbh >= 2.0, year == 2014)
anti_join(dbh2014, dbh_a, by = "tag") # this inidcates that there are 126 stems that were measured and tagged for the first time in 2014.
dbh_a <- union(dbh_a, anti_join(dbh2014, dbh_a, by = "tag"))
anyDuplicated(dbh_a)
# At this point the 'dbh_a' dataframe should be a record of the first season each tagged stem was measured during the study period. In theory, these things happened at the same time (tagging and measuring DBH).
summary(dbh_a)

dbh_a <- rename(dbh_a, dbh_entry = dbh)
tagged_stems <- as.tbl(left_join(tagged_stems,
                                 select(dbh_a, plotid, tag, dbh_entry),
                                 by = c("plotid", "tag")))
summary(tagged_stems)
summary(filter(tagged_stems, is.na(dbh_entry)))# NA are due to stems not meeting the 2.0 cm DBH requirement.

# dropping stems that don't meet DBH threshold requirement
tagged_stems <- filter(tagged_stems, dbh >= 2.0)
write.csv(tagged_stems, "analysis/data/tagged-stems-corrected.csv")


## Explore DBH data: negative change in DBH ####
library(ggplot2)
# Scatterplot of first dbh measurement (y-axis) against year
qplot(year, dbh_entry, data = dbh_a, geom="jitter", color = stem_status)
# The majority of all stems had been measured by 2005

# Dotplot of most recent dbh measurement
qplot(stem_status, dbh, data = dbh2014, geom = "jitter", color = stem_status)

# Data frame of tagged stems that were measured or remeasured in 2014
unique(dbh_establishment$year)
dbh_establishment$year <- 2005
unique(dbh2014$year)

# Match plots between establishment and 2014 remeasurement
length(unique(dbh_establishment$plotid))# 203 plots
length(unique(dbh2014$plotid))# 195 plots
tmp1 <- droplevels(unique(dbh_establishment$plotid))
tmp2 <- droplevels(unique(dbh2014$plotid))

dbh_changes <- rbind(dbh_establishment, dbh2014)
summary(dbh_changes)
qplot(factor(year), dbh, data = dbh_changes, geom = "jitter", color = stem_status, facets = .~species)

