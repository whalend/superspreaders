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
write.csv(tagged_stems, "analysis/data/tagged_stems_revised.csv", row.names = F)
write.csv(untagged_stems, "analysis/data/untagged_stems_revised.csv", row.names = F)
write.csv(all_stems <- rbind(
      select(tagged_stems, plotid, species, dbh, stem_status, date, year),
      untagged_stems),
      "analysis/data/all_stems_revised.csv", row.names = F) 

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
write.csv(dbh_establishment, "analysis/data/tagged-stems-2003_2005.csv", row.names = F)


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
summary(too_small_stems <- filter(tagged_stems, is.na(dbh_entry)))# NA are due to stems not meeting the 2.0 cm DBH requirement, except at least one curiosity for the initial_dbh
filter(tagged_stems, is.na(dbh_entry), initial_dbh>2)
# this is a plot that was abandoned very early in the study, may have been visited only one time
# Take a closer look at these too small stems
# fix(too_small_stems)
# many are from abandoned plots, or are missing stems, and mostly did not meet DBH of 2.0cm upon first measurement

# dropping 194 stem records that don't meet DBH threshold requirement on their entry DBH measurement
tagged_stems_corrected <- filter(tagged_stems, is.na(dbh_entry)==F)


## Examine untagged stems data ####
summary(untagged_stems)
length(unique(untagged_stems$plotid))# 185 plots with untagged stems recorded
length(unique(tagged_stems_corrected$plotid))# 203 plots
unique(untagged_stems$plotid) %in% unique(tagged_stems_corrected$plotid)
# there are three plots in the untagged stems data that aren't in the corrected tagged stems data
unique(untagged_stems$plotid)[c(28,163,170)]
# badge01 does not have any tagged hosts; ann28 and halow01 were abandoned in 2003/2004
# remove the abandoned plots, retain badge01
untagged_stems <- filter(untagged_stems, plotid != "ann28", plotid != "halow01")
write.csv(untagged_stems, "analysis/data/untagged_stems_corrected.csv", row.names = F)


## Explore DBH data ####
library(ggplot2)
# Scatterplot of first dbh measurement (y-axis) against year
qplot(year, dbh_entry, data = dbh_a, geom="jitter", color = stem_status)
# The majority of all stems had been measured by 2005

# Dotplot of most recent dbh measurement
qplot(stem_status, dbh, data = dbh2014, geom = "jitter", color = stem_status)

## Match plots between establishment and 2014 remeasurement ####
unique(dbh_establishment$year)
dbh_establishment$year <- 2005
unique(dbh2014$year)

length(unique(dbh_establishment$plotid))# 203 plots
length(unique(dbh2014$plotid))# 195 plots
dbh_establishment$plotid <- as.factor(dbh_establishment$plotid)
tmp1 <- droplevels(unique(dbh_establishment$plotid))
dbh2014$plotid <- as.factor(dbh2014$plotid)
tmp2 <- droplevels(unique(dbh2014$plotid))
setdiff(tmp1,tmp2)# identify plots not in 2014 data
# These all make sense except for mroth03. For some reason, the 2014 plot visit is not being included in the Access database query.

# Add the plot and stems, and update the DBH measurements according to the values in the database
dbh2014 <- rbind(dbh2014, dbh2012 %>% filter(plotid == "mroth03"))
# enter 2014 DBH values from database
dbh2014$dbh[dbh2014$tag==1668] <- 68
dbh2014$dbh[dbh2014$tag==1669] <- 49.2
dbh2014$dbh[dbh2014$tag==1670] <- 56.1
dbh2014$dbh[dbh2014$tag==1671] <- 51.2
dbh2014$dbh[dbh2014$tag==1672] <- 66.9
dbh2014$year <- 2014

# These measurements for 2014 are then also missing from the tagged_stems dataframe
filter(dbh2014, plotid == "mroth03")
summary(filter(tagged_stems, plotid == "mroth03") %>% 
              select(plotid, tag, date))
names(tagged_stems)
names(dbh2014)
mroth_tmp1 <- filter(tagged_stems, plotid == "mroth03", year == 2011) %>% 
      select(-plot_status, -species, -location, -dbh, -stem_status, -year)
mroth_tmp <- filter(dbh2014, plotid == "mroth03")
mroth_tmp <- left_join(mroth_tmp, mroth_tmp1)
names(mroth_tmp)
mroth_tmp <- select(mroth_tmp, plotid, plot_status, cluster, species, tag, location, dbh, initial_dbh, stem_status, sod_killed:date, year, dbh_entry)
tagged_stems <- rbind(tagged_stems, mroth_tmp)
tagged_stems_corrected <- rbind(tagged_stems_corrected, mroth_tmp)
rm(mroth_tmp); rm(mroth_tmp1)
write.csv(tagged_stems_corrected, "analysis/data/tagged-stems-corrected.csv", row.names = F)

# revise all_stems dataframe with corrected tagged and untagged stems
all_stems_corrected <- rbind(
      select(tagged_stems_corrected, plotid, species, dbh, stem_status, date, year),
      untagged_stems)
summary(all_stems_corrected)
write.csv(all_stems_corrected, "analysis/data/all_stems_corrected.csv", row.names = F)

# Reassess the plots not in 2014 data
tmp1 <- droplevels(unique(dbh_establishment$plotid))
tmp2 <- droplevels(unique(dbh2014$plotid))
setdiff(tmp1,tmp2)
# These make sense as plots that were not resampled in 2014 due to being lost for various resons earlier in the study. For assessing changes in basal area and abundance at the two time steps I want only the plots that were measured during both.

# Remove plots from the establishment data set that were not revisited in 2014
dbh_establishment <- filter(dbh_establishment, plotid != "bush01", 
                            plotid != "mcneil01", plotid != "mcneil03", 
                            plotid != "ponti01", plotid != "sumtv02", 
                            plotid != "sweet01", plotid != "votru01")

## Combine establishment and 2014 remeasurement tagged stems dataframes ####
tagged_dbh_changes <- rbind(dbh_establishment, dbh2014)
summary(tagged_dbh_changes)
tagged_dbh_changes <- select(tagged_dbh_changes, -plot_status)
tagged_dbh_changes$species <- as.factor(tagged_dbh_changes$species)
tagged_dbh_changes$year <- as.factor(tagged_dbh_changes$year)

tagged_dbh_changes <- left_join(tagged_dbh_changes,
      as.tbl(tagged_dbh_changes %>% filter(year == 2005) %>% 
                   rename(dbh1 = dbh) %>% 
                   select(plotid, tag, dbh1)))
tagged_dbh_changes <- left_join(tagged_dbh_changes,
      as.tbl(tagged_dbh_changes %>% filter(year == 2014) %>% 
                   rename(dbh2 = dbh) %>% 
                   select(plotid, tag, dbh2)))
tagged_dbh_changes <- tagged_dbh_changes %>% 
      mutate(delta_dbh = dbh2 - dbh1, abs_incr_grwth = dbh2/dbh1, 
             inst_grwth_rate = log(dbh2/dbh1) / (2014 - 2005))
summary(tagged_dbh_changes)

# export to data file
write.csv(tagged_dbh_changes, "analysis/data/tag-dbh_0514-corrected.csv", row.names = F)

## Exploratory plots of tagged stem DBH changes ####
qplot(delta_dbh, dbh1, data = tagged_dbh_changes %>% 
            group_by(plotid, tag, species, year) %>% 
            filter(year == 2014),
      color = stem_status, facets = species ~.)
filter(tagged_dbh_changes, delta_dbh < -15, year == 2014, stem_status == "Alive", location == "In")

qplot(year, abs_incr_grwth, data = tagged_dbh_changes, geom = "boxplot")

## Compare Basal & Species Abundances ####
par(mfrow=c(2,1))
boxplot(dbh ~ year, 
        data = tagged_dbh_changes %>% 
              filter(species == "UMCA", stem_status == "Alive"), 
        main = "Live UMCA DBHs")
boxplot(dbh ~ year, 
        data = tagged_dbh_changes %>% 
              filter(species == "UMCA", stem_status == "Dead"), 
        main = "Dead UMCA DBHs")

qplot(factor(year), dbh, data = tagged_dbh_changes, geom = "jitter", color = stem_status, facets = .~species)

# basal_abund <- dbhs %>% 
#       #filter(status == "Alive") %>% 
#       group_by(plotid, species, year, status) %>% 
#       summarise(avg_ba_m2 = mean(pi*dbh^2/40000, na.rm = T), 
#                 tot_ba_m2 = sum(pi*dbh^2/40000, na.rm = T),
#                 abundance = length(dbh))
# summary(basal_abund)
# write.csv(basal_abund, "analysis/data/tag-dbh-plot.csv", row.names = F)


## Assess untagged stems data ####
summary(untagged_stems)
unique(untagged_stems$year)# 5 unique years
# I need to figure out if things might be duplicated between the years
tmp03 <- filter(untagged_stems, year == 2003)
tmp04 <- filter(untagged_stems, year == 2004)
tmp05 <- filter(untagged_stems, year == 2005)

unique(tmp03$plotid) %in% unique(tmp04$plotid)
unique(tmp03$plotid)[c(28:30)]
filter(tmp03, plotid == "ganay02"); filter(tmp04, plotid == "ganay02")
filter(untagged_stems, year == 2012, plotid == "ganay02")
# this would support keeping 2004 data, removing 2003 data
untagged_stems <- filter(untagged_stems, year != 2003 | plotid != "ganay02")

filter(tmp03, plotid == "lins01"); filter(tmp04, plotid == "lins01")
filter(untagged_stems, year == 2012, plotid == "lins01")
# the 2003 and 2004 data appears to be unique, so leave alone for now

filter(tmp03, plotid == "mitsu02"); filter(tmp04, plotid == "mitsu02")
filter(untagged_stems, year == 2012, plotid == "mitsu02")
# little difference between 2003 and 2004; keep 2004 data
untagged_stems <- filter(untagged_stems, year != 2003 | plotid != "mitsu02")

rm(tmp03)

# set all 2003 year to 2004
untagged_stems$year[untagged_stems$year == 2003] <- 2004
tmp04 <- filter(untagged_stems, year == 2004)

# compare 2005 plot coverage to 2004/2003 plot coverage
unique(tmp05$plotid) %in% unique(tmp04$plotid)# two duplicates
unique(tmp05$plotid)[c(25,26)]
filter(tmp04, plotid == "cook02"); filter(tmp05, plotid == "cook02")
filter(untagged_stems, year == 2014, plotid == "cook02")
# assume the 2004 data is correct for this plot
untagged_stems <- filter(untagged_stems, year != 2005 | plotid != "cook02")

filter(tmp04, plotid == "jlsp08"); filter(tmp05, plotid == "jlsp08")
filter(untagged_stems, year == 2014, plotid == "jlsp08")
# did we miss a bunch of SESE in 2012 & 2014??? I'm dropping those from the 2004 data
untagged_stems <- filter(untagged_stems, year != 2005 | plotid != "jlsp08")
untagged_stems <- filter(untagged_stems, year != 2004 | plotid != "jlsp08"|
                         species != "sese")
rm(tmp04); rm(tmp05)

#set all 2004 year to 2005
untagged_stems$year[untagged_stems$year == 2004] <- 2005

# compare 2012 plot coverage to 2014 plot coverage
tmp12 <- filter(untagged_stems, year == 2012)
tmp14 <- filter(untagged_stems, year == 2014)
unique(tmp12$plotid) %in% unique(tmp14$plotid)
unique(tmp12$plotid)[c(4,12,13,20,21,52,53,58,111,114)]
# 10 plots not in the 2014 data that should be, i.e. none of these were abandoned/decommissioned and should have been visited in 2014
tmp12 <- filter(tmp12, plotid == "arbit01" | plotid == "lins01" |
                      plotid == "lupin01" | plotid == "mitsu02" |
                      plotid == "skile01" | plotid == "spaul01" |
                      plotid == "ann09" | plotid == "ganay02" |
                      plotid == "halow02")
untagged_stems <- filter(untagged_stems, year != 2012)
untagged_stems <- rbind(untagged_stems, tmp12)
untagged_stems$year[untagged_stems$year == 2012] <- 2014
write.csv(untagged_stems, "analysis/data/untagged-stems_0514-corrected.csv", row.names = F)



