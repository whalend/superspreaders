--------
# Exploring UMCA stem infection through time #
# Superspreaders chapter of dissertation      
--------

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
dbh2014 <- tagged_stems %>% filter(dbh>0, year == 2014)
length(unique(dbh2014$tag))
# 4158, indicating there are 294 tagged_stems with no dbh measured in 2014

dbh2003 <- tagged_stems %>% select(plot, tag, dbh, year) %>% 
      filter(dbh>0, year == 2003) # these 2 are the same as in 2005
# 11 observations

dbh2004 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2004)
# 3109 observations
anti_join(dbh2003, dbh2004, by = "tag")
dbh_a <- union(dbh2004, dbh2003)


dbh2005 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2005)
# 712 observations
anti_join(dbh2005, dbh_a, by = "tag") # this indicates that the dbh data for 2004 and 2005 are unique
dbh_a <- union(dbh_a, anti_join(dbh2005, dbh_a, by = "tag"))

dbh2006 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2006)
anti_join(dbh2006, dbh_a, by = "tag") # 2004, 2005, 2006 are independent sets
dbh_a <- union(dbh_a, anti_join(dbh2006, dbh_a, by = "tag"))

dbh2007 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2007)
anti_join(dbh2007, dbh_a, by = "tag") # there is one tag with remeasured dbh in 2007
dbh_a <- union(dbh_a, anti_join(dbh2007, dbh_a, by = "tag"))
length(unique(dbh_a$tag)) # check for repeated tags, if obs > value, then there is a duplicate

dbh2008 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2008)
anti_join(dbh2008, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2008, dbh_a, by = "tag"))

dbh2009 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2009)
anti_join(dbh2009, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2009, dbh_a, by = "tag"))

dbh2010 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2010)
anti_join(dbh2009, dbh_a, by = "tag") # no unique dbh values for 2010, so no new recruitment, or tagged_stems entering the study

dbh2011 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2011)
anti_join(dbh2011, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2011, dbh_a, by = "tag"))

dbh2012 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2012)
anti_join(dbh2012, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2012, dbh_a, by = "tag"))

dbh2014 <- tagged_stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2014)
anti_join(dbh2014, dbh_a, by = "tag") # this inidcates that there are 113 tagged_stems that were measured for the first time in 2014

#### Explore DBH data: negative change in DBH ####
library(ggplot2)
# Scatterplot of first dbh measurement (y-axis) against year
qplot(year, dbh, data = dbh_a, geom="jitter")
# Scatterplot of most recent dbh measurement
qplot(year, dbh, data = dbh2014, geom = "jitter")

# Data frame of tagged_stems that were measured or remeasured in 2014
dbh <- right_join(dbh_a, dbh2014, by = c("plot","tag"))
summary(dbh)
sum(dbh$dbh.x, na.rm = TRUE)
sum(dbh$dbh.y, na.rm = TRUE)
sum(dbh$dbh.y, na.rm = TRUE) - sum(dbh$dbh.x, na.rm = TRUE)
# [1] 6345.03 centimeters of total growth in diameter for all tagged_stems remeasured in 2014

dbh <- rename(dbh, dbh1 = dbh.x, dbh2 = dbh.y, year_dbh1 = year.x, year_dbh2 = year.y) # This data frame has all the tagged_stems for which we have an initial measurement AND a remeasurement in 2014
dbh$delta_dbh <- dbh$dbh2 - dbh$dbh1
write.csv(dbh, "analysis/data/dbh_remeasures.csv")
summary(dbh)
# filter(dbh, delta_dbh < -20)

library(ggplot2)
qplot(dbh1, dbh2, data = dbh)
qplot(delta_dbh, dbh2, data = dbh, color = year_dbh1)

#Join dbh measure, remeasure and change data to `tagged_stems` data frame ####
dbh <- read.csv("analysis/data/dbh_2014_remeasures.csv")
tagged_stems <- left_join(tagged_stems, dbh, by = c("tag","plot"))
str(tagged_stems)
summary(tagged_stems)
qplot(delta_dbh, dbh2, data = tagged_stems %>% select(plot, cluster, tag, species, delta_dbh, dbh2, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead"), 
      color = status)

qplot(year, dbh, data = tagged_stems %>% filter(status == "Alive" | status == "Dead"), color = status, geom = "jitter")


tagged_stems %>% select(plot, cluster, tag, species, delta_dbh, year, year_dbh1, status) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(delta_dbh < -5, year == 2014, status == "Alive")
# I changed values for 9 of these 14 tagged_stems

tagged_stems %>% select(plot, cluster, tag, species, delta_dbh, year, year_dbh1, status) %>% group_by(plot, cluster, tag) %>% 
      filter(between(delta_dbh, -15, -10), year == 2014, status == "Alive")
# I changed values for 7 of 9 tagged_stems

tagged_stems %>% select(plot, cluster, tag, species, delta_dbh, year, year_dbh1, status) %>% group_by(plot, cluster, tag) %>% 
      filter(between(delta_dbh, -10, -7), year == 2014, status == "Alive")
# I changed values for 11 of 20 tagged_stems

tagged_stems %>% select(plot, cluster, tag, species, delta_dbh, year, year_dbh1, status) %>% group_by(plot, cluster, tag) %>% 
      filter(between(delta_dbh, -6.9, -5), year == 2014, status == "Alive")
# I changed values for 6 of 18 tagged_stems


#### Create a plot-level data frame of infection status ####
tagged_stems <- read.csv("analysis/data/all_tagged_stems_corrected.csv")
library(dplyr)

plots_umca <- tagged_stems %>%
      select(plot, Date, species, slc, year, status) %>%
      group_by(plot, year, species) %>%
      filter(species == "UMCA", status == "Alive") %>%
      summarise(uninfected_bay_ct = length(which(slc==0)), infected_bay_ct = length(which(slc > 0)), tot_bay = length(species)) 

# plots_umca$infected <- ifelse(plots_umca$infected_bay_ct==0, 0, 1)
summary(plots_umca)
plots_umca <- droplevels(plots_umca)


length(which(plots_umca$ct_bay_NA > 0))
# filter(plots_umca, infected == 0)
summary(filter(plots_umca, infected_bay_ct == 0))# 49 plots, max year = 2010

write.csv(plots_umca, "analysis/data/plots_umca_infection.csv")
# I sent this file to Francesco for informing the spread model


plot_qusp <- tagged_stems %>%
      select(plot, date, species, canker, sod_dead, year, status) %>%
      group_by(plot, year, species) %>%
      filter(species == "QUAG" | species == "QUKE" | 
                   species == "QUWI" | species == "QUCH", 
             status == "Alive" | status == "Dead") %>%
      summarise(uninfected_oak_ct = length(which(canker==0)), infected_oak_ct = length(which(canker > 0)), tot_oak = length(species))
summary(plot_qusp)
plot_qusp <- droplevels(plot_qusp)

library(tidyr)
spread(plot_qusp, key = species, value = infected_oak_ct, drop = TRUE)

gather(plot_qusp, inf_status, count, ends_with("ct")) %>% 
      unite(species, inf_status) %>% 
      spread(inf_status, )

library(reshape2)
tmp <- plot_qusp %>%
      dcast(year + plot ~ species, value.var = "infected_oak_ct")
head(tmp)

# Checking large positive DBH changes ####
## So at this point I think I have corrected all the negative dbh change values that I legitimately could in the database. I have gone back and exported the query again and use this file below, and then recycling code from the beginning of this document.
tagged_stems <- read.csv("analysis/data/Stem_Summary_negDBH_corrected_201402.csv")
summary(tagged_stems)
str(tagged_stems)

tagged_stems$Date <- as.character(tagged_stems$Date)
class(tagged_stems$Date)
tagged_stems$date <- as.Date(tagged_stems$Date, format = "%m/%d/%Y")
library(lubridate)
tagged_stems$year <- year(tagged_stems$date)
summary(tagged_stems)
library(plyr)
library(dplyr)
filter(tagged_stems, is.na(year))
filter(tagged_stems, PlotID == "SECRET02", SpeciesID == "QUAG")

tagged_stems$year[is.na(tagged_stems$year)] <- 2012
filter(tagged_stems, is.na(date))
tagged_stems$date[is.na(tagged_stems$date)] <- "2012-05-22"
summary(tagged_stems)
tagged_stems[tagged_stems == -9999] <- NA
summary(tagged_stems)

tagged_stems <- rename(tagged_stems, plotid = PlotID)
tagged_stems <- rename(tagged_stems, cluster = HostID, tag = TagNumber, 
                species = SpeciesID, status = tagged_stemstatus, 
                alive_class = AliveClass, dead_class= DeadClass,
                slc = SympLeafCount, canker = CankerPresent, dbh = DBH,
                lide_lf_symp = LideFoliarSymp, sod_dead = SOD_Dead, 
                location = Location)
str(tagged_stems)
tagged_stems <- rename(tagged_stems, notes = tagged_stemsymptomNotes)
summary(tagged_stems$status)
tagged_stems$status[tagged_stems$status=="Nearly dead"] <- "Alive"
tagged_stems <- droplevels(tagged_stems)
class(tagged_stems)
tagged_stems <- as.tbl(tagged_stems)

# Subset DBH data to create a data frame of remeasured tagged_stems: Round 2 ####

detach("package:lubridate", unload=TRUE)# b/c it has `union` function masking the one from `dplyr` that apparently behaves a little differently, creating a list instead of a data frame


dbh2003 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>% 
      filter(dbh>0, year == 2003) # these 2 are the same as in 2005
# 11 observations

dbh2004 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2004)
# 3109 observations
anti_join(dbh2003, dbh2004, by = c("plotid", "tag"))# 11 DBH observations different between 2003 and 2004
dbh_a <- union(dbh2004, dbh2003)

dbh2005 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2005)
# 702 observations
anti_join(dbh2005, dbh_a, by = c("plotid","tag")) # this indicates that the dbh data for 2003/04 and 2005 are nearly unique, 702 new DBH values at 109 plots in 2005
# The 'union' with the 'anti_join' selects only the new records from 2005
dbh_a <- union(dbh_a, anti_join(dbh2005, dbh_a, by = c("plotid","tag")))
unique(dbh_a$plotid)# 204 plots on record at this point
dbh_t1 <- dbh_a# data frame for DBH at time 1

dbh2006 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2006)
anti_join(dbh2006, dbh_a, by = c("plotid","tag")) # 2004, 2005, 2006 are independent sets, 86 newly tagged DBH values in 2006
dbh_a <- union(dbh_a, anti_join(dbh2006, dbh_a, by = c("plotid","tag")))

dbh2007 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2007)
anti_join(dbh2007, dbh_a, by = c("plotid","tag")) # there is one tag with remeasured DBH in 2007 and 36 new DBH observations
dbh_a <- union(dbh_a, anti_join(dbh2007, dbh_a, by = c("plotid","tag")))
anyDuplicated(dbh_a$tag)# check for any repeated tags, stops at first encounter

dbh2008 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2008)
anti_join(dbh2008, dbh_a, by = c("plotid","tag"))# 12 new DBH observations
dbh_a <- union(dbh_a, anti_join(dbh2008, dbh_a, by = c("plotid","tag")))

dbh2009 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2009)
anti_join(dbh2009, dbh_a, by = c("plotid","tag"))# 21 new DBH observations
dbh_a <- union(dbh_a, anti_join(dbh2009, dbh_a, by = c("plotid","tag")))

dbh2010 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2010)
anti_join(dbh2009, dbh_a, by = c("plotid","tag")) # no unique DBH values for 2010, so no new recruitment or tagged_stems entering the study, just remeasurement of 109 tagged_stems

dbh2011 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2011)
anti_join(dbh2011, dbh_a, by = c("plotid","tag"))# 74 new DBH observations
dbh_a <- union(dbh_a, anti_join(dbh2011, dbh_a, by = c("plotid","tag")))

dbh2012 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2012)# remeasured all tagged_stems
anti_join(dbh2012, dbh_a, by = c("plotid","tag"))# 289 new DBH observations (maybe...)
dbh_a <- union(dbh_a, anti_join(dbh2012, dbh_a, by = c("plotid","tag")))
anyDuplicated(dbh_a$tag) # check for repeated tags
dbh_a[1091,]
filter(dbh_a, tag == 4042)
filter(dbh_a, plotid == "SUGAR02")
filter(dbh_a, plotid == "SUGAR30")
filter(dbh_a, plotid == "SUGAR02")$tag %in% filter(dbh_a, plotid == "SUGAR30")$tag
# correct based on database, where SUGAR02 was abandoned in 2004
dbh_a$plotid[dbh_a$plotid=="SUGAR02" & dbh_a$year == 2004] <- "SUGAR30"

tmp <- filter(dbh_a, plotid == "SUGAR30", year == 2012, tag != 7302, tag != 7303, tag != 7301)
dbh_a <- anti_join(dbh_a, tmp); rm(tmp)
anyDuplicated(dbh_a$tag)
# Final count of 4329 tagged_stems after 2012, up from 3822 after 2005.
length(unique(dbh_a$plotid)); length(unique(dbh_t1$plotid))


dbh2014 <- tagged_stems %>% select(plotid, tag, species, status, dbh, year, date, slc, canker, sod_dead, lide_lf_symp, location) %>%
      filter(dbh>0, year == 2014)
anyDuplicated(dbh2014$tag)# 1852 duplicated tags, WTF? 
filter(dbh2014, duplicated(tag))# nope, that's the row number
filter(dbh2014, tag == 4194)
dbh2014 <- filter(dbh2014, tag != 4194 | date != "2012-05-14")
summary(dbh2014)

anti_join(dbh2014, dbh_a, by = "tag") # this inidcates that there are 113 tagged_stems that were measured for the first time in 2014

no2014 <- anti_join(dbh_a, dbh2014, by = "tag")# this indicates that there were 284 tagged_stems previously measured that were not part of the remeasurement in 2014. This looks in part to be due to plot abandonments/decomissions, but may also be due to other reasons. It looks like a few of these had measurements < 2.0 cm, and the McNeil plots were removed by the landowner before 2014.

# 'dbh_a' has all measurements through 2012
# 'dbh2014' has all remeasurements in 2014, some plots lost at this point
# As of 2004 there were 3120 tagged_stems, by 2005 there were 3824 tagged_stems, and then after 2012 there were 4329 tagged_stems, but 284 of these were not remeasured in 2014.

summary(no2014); unique(no2014$plotid)# 73 plots; for proper comparison, whether using 2012 or 2014 data, I need to limit the set to plots that were in the establishment and final sampling.

# Compare Basal & Species Abundances ####
head(dbh_t1)
dbh_t1 <- select(dbh_t1, -date) %>% arrange(plotid)
dbh_t1$year <- 2005# change all to last sampling season included
dbh_t1 <- droplevels(dbh_t1)
head(dbh2014)
dbh2014 <- select(dbh2014, -date) %>% arrange(plotid)
dbh2014 <- droplevels(dbh2014)

# Match plots between establishment and 2014 remeasurement
length(unique(dbh_t1$plotid))# 204 plots
length(unique(dbh2014$plotid))# 195 plots
tmp1 <- droplevels(unique(dbh_t1$plotid))
tmp2 <- droplevels(unique(dbh2014$plotid))
setdiff(tmp1,tmp2)# identify plots not in 2014 data, note that one is SUGAR02
filter(dbh_t1, plotid == "SUGAR02" | plotid == "SUGAR30")
# change SUGAR02 to SUGAR30 for early data
dbh_t1$plotid[dbh_t1$plotid == "SUGAR02"] <- "SUGAR30"

tmp1 <- droplevels(unique(dbh_t1$plotid))
tmp2 <- droplevels(unique(dbh2014$plotid))
setdiff(tmp1,tmp2)
dbh2014 <- rbind(dbh2014, select(dbh2012, -date) %>% 
                       filter(plotid == "MROTH03"))
# enter 2014 DBH values from database
dbh2014$dbh[dbh2014$tag==1668] <- 68
dbh2014$dbh[dbh2014$tag==1669] <- 49.2
dbh2014$dbh[dbh2014$tag==1670] <- 56.1
dbh2014$dbh[dbh2014$tag==1671] <- 51.2
dbh2014$dbh[dbh2014$tag==1672] <- 66.9
dbh2014$year <- 2014

tmp1 <- droplevels(unique(dbh_t1$plotid))
tmp2 <- droplevels(unique(dbh2014$plotid))
setdiff(tmp1,tmp2)# reduced to 7 missing, 3 lost "late" in study
# filter out plots from 2005 data that aren't in 2014 data
dbh_t1 <- filter(dbh_t1, plotid != "BUSH01", plotid != "MCNEIL01",
                 plotid != "MCNEIL03", plotid != "PONTI01",
                 plotid != "SUMTV02", plotid != "SWEET01",
                 plotid != "VOTRU01")

dbhs <- rbind(dbh_t1, dbh2014)
str(dbhs)
unique(dbhs$plotid)
dbhs <- droplevels(dbhs)
summary(dbhs)
dbhs <- filter(dbhs, status == "Alive" | status == "Dead")
write.csv(dbhs, "analysis/data/tag-dbh_0514-corrected.csv", row.names = F)

par(mfrow=c(2,1))
boxplot(dbh ~ year, 
        data = dbhs %>% filter(species == "UMCA", status == "Alive"), 
        main = "Live UMCA DBHs")
boxplot(dbh ~ year, 
        data = dbhs %>% filter(species == "UMCA", status == "Dead"), 
        main = "Dead UMCA DBHs")


basal_abund <- dbhs %>% 
      #filter(status == "Alive") %>% 
      group_by(plotid, species, year, status) %>% 
      summarise(avg_ba_m2 = mean(pi*dbh^2/40000, na.rm = T), 
                tot_ba_m2 = sum(pi*dbh^2/40000, na.rm = T),
                abundance = length(dbh))
summary(basal_abund)
write.csv(basal_abund, "analysis/data/tag-dbh-plot.csv", row.names = F)


## Some plots of DBH measurements ####
library(ggplot2)
# Scatterplot of first dbh measurement (y-axis) against year
qplot(tag, dbh, data = dbh_a, geom="jitter", facets = year ~.)

# Scatterplot of most recent dbh measurement
qplot(tag, dbh, data = dbh2014, geom = "jitter", facets = year ~.)
qplot(dbh, data = dbh2014)

# Data frame of tagged_stems that were measured or remeasured in 2014
dbh <- right_join(dbh_a, dbh2014, by = c("plotid","tag"))
summary(dbh)

dbh <- rename(dbh, dbh1 = dbh.x, dbh2 = dbh.y, year_dbh1 = year.x, year_dbh2 = year.y) # This data frame has all the tagged_stems for which we have an initial measurement AND a remeasurement in 2014
dbh$delta_dbh <- dbh$dbh2 - dbh$dbh1
write.csv(dbh, "analysis/data/dbh_2014_remeasures.csv")
summary(dbh)
sum(dbh$delta_dbh, na.rm=T)
# filter(dbh, delta_dbh < -20)


#### Check the positive DBH change outliers of tagged_stems remeasured in 2014 ####
dbh <- read.csv("analysis/data/dbh_2014_remeasures.csv")
summary(dbh)
library(ggplot2)
qplot(dbh1, dbh2, data = dbh)
qplot(delta_dbh, dbh2, data = dbh, color = year_dbh1)
summary(dbh)
dbh <- as.tbl(dbh)

# Join dbh measure, remeasure, and change data to `tagged_stems` data frame
dbh <- select(dbh, -X) %>%
      rename(date_dbh1 = date.x, date_dbh2 = date.y)
summary(dbh)
summary(tagged_stems)
tagged_stems <- left_join(tagged_stems, dbh, by = c("tag","plot"))
str(tagged_stems)
summary(tagged_stems)

# Exploratory plot of dbh data
qplot(delta_dbh, dbh1, data = tagged_stems %>% select(plotid, cluster, tag, species, delta_dbh, dbh1, year, status) %>% 
            group_by(plotid, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead", species == "UMCA"), 
      color = status)
# A visual check shows that there are a handful of tagged_stems with growth > 20cm

#### Suggestions from Margaret about exploring DBH changes ####
# I emailed Richard and Margaret about thresholds for where I should really look more closely at the measurements. Richard noted that there isn't much out there on growth rates for these species. I did find something done for tree species of the Northeastern U.S. <http://www.fs.fed.us/ne/newtown_square/publications/research_papers/pdfs/scanned/OCR/ne_rp649.pdf>

# Margaret suggested a few things, including other ways of looking at the data:
# 1. "Have you tried looking at relative growth instead of the absolute increment growth?  Either DBH2/DBH1, or the instantaneous rate ln(DBH2/DBH1) / (time2-time1).  20 cm growth means a lot more to a 5 cm tree than a 50 cm tree… either way, it's still big.  I'd also consider this species by species, because I bet some grow much faster than others."

tagged_stems$dbh2_1_ratio <- tagged_stems$dbh2/tagged_stems$dbh1# absolute increment growth
tagged_stems$inst_grwth_rate <- log(tagged_stems$dbh2/tagged_stems$dbh1) / (tagged_stems$year_dbh2 - tagged_stems$year_dbh1)
#tagged_stems$time_diff <- (tagged_stems$year_dbh2 - tagged_stems$year_dbh1)
summary(tagged_stems)
library(lubridate)
tagged_stems$doy <- yday(tagged_stems$date)# create a day of year variable
write.csv(tagged_stems, "analysis/data/tagged_stems_growth.csv")
summary(tagged_stems$status)
summary(tagged_stems$tag)

qplot(dbh2_1_ratio, data = tagged_stems %>% select(plotid, tag, species, dbh2_1_ratio, dbh1, year, status, doy) %>% 
            group_by(plotid, tag) %>% 
            filter(status == "Alive" | status == "Dead", dbh1 >= 2.0),
      color = status)

qplot(tagged_stems$dbh2_1_ratio, main = "Relative Growth")

qplot(dbh1, dbh2_1_ratio, data = tagged_stems %>% select(plotid, cluster, tag, species, dbh2_1_ratio, dbh1, year, status, dbh) %>% 
            group_by(plotid, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead", dbh1 >= 2.0), 
      color = status, main = "Change in DBH vs. Original DBH")

qplot(dbh1, inst_grwth_rate, data = tagged_stems %>% select(plotid, cluster, tag, species, inst_grwth_rate, dbh1, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead", dbh1 >= 2.0), 
      color = status, main = "Instantaneous Growth Rate vs. Original DBH")

qplot(dbh1, inst_grwth_rate, data = tagged_stems %>% select(plotid, cluster, tag, species, inst_grwth_rate, dbh1, year_dbh1, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(species == "QUAG", year == 2014, status == "Dead" | status == "Alive", dbh1 >= 2.0), 
      color = status, main = "QUAG Instantaneous Growth Rate vs Original DBH")

qplot(dbh1, inst_grwth_rate, data = tagged_stems %>% select(plotid, cluster, tag, species, inst_grwth_rate, dbh1, year_dbh1, year, status) %>% 
            group_by(plotid, cluster, tag) %>%
            filter(species == "UMCA", year == 2014, status == "Dead" | status == "Alive", dbh1 >= 2.0), 
      color = status, main = "UMCA Instantaneous Growth Rate vs Original DBH")

boxplot(dbh ~ year, data = tagged_stems %>% filter(species == "UMCA", status == "Alive", year == 2005|year==2012|year == 2014))

tagged_stems %>% filter(species == "UMCA", status == "Alive", year == 2005|year==2012|year == 2014)

u.2005 <- tagged_stems %>%
      filter(species == "UMCA", status == "Alive", year == 2005|year==2014)

t.test(dbh~year, data = tagged_stems %>% filter(species == "UMCA", status == "Alive", year_dbh1 == 2005|year_dbh2 == 2014))
t.test(dbh~year, data = tagged_stems %>% filter(species == "QUAG", status == "Alive", year == 2005|year == 2014))
t.test(dbh~year, data = tagged_stems %>% filter(species == "QUAG", status == "Dead", year == 2005|year == 2014))
t.test(dbh~year, data = tagged_stems %>% filter(species == "QUKE", status == "Alive", year == 2005|year == 2014))
t.test(dbh~year, data = tagged_stems %>% filter(species=="QUAG"|species == "QUKE", status == "Alive", year == 2005|year == 2014))



tagged_stems %>% select(plot, cluster, tag, species, year, year_dbh1, status, delta_dbh, dbh1, dbh2, notes) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(delta_dbh > 20, year == 2014, status == "Alive")
# 8 tagged_stems

tagged_stems %>% select(plot, tag, cluster, species, year, year_dbh1, status, dbh, delta_dbh, dbh1, dbh2, notes) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(plot == "LAZY05", species == "QUKE", status == "Alive")

tagged_stems %>% select(plot, tag, cluster, species, year, year_dbh1, status, dbh, delta_dbh, dbh1, dbh2, notes) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(plot == "BOWES02", species == "UMCA", status == "Alive")

tagged_stems %>% arrange(desc(delta_dbh)) %>% 
      select(plot, cluster, tag, species, delta_dbh, dbh1, dbh2, year, year_dbh1, status) %>% group_by(plot, cluster, tag) %>% 
      filter(between(delta_dbh, 10, 20), year == 2014, status == "Alive")
# 39 tagged_stems
summary(tagged_stems)
rm(list = ls(pattern = "dbh20*"))

## Exploring Symptomatic Leafcount Data ####
tagged_stems <- read.csv("analysis/data/tagged_stems_growth.csv")
str(tagged_stems)
summary(tagged_stems)
tagged_stems$date <- as.Date(tagged_stems$date)
library(plyr);library(dplyr)
tagged_stems <- as.tbl(tagged_stems)
summary(tagged_stems)
tagged_stems <- select(tagged_stems, -X, -Date)
summary(tagged_stems$date)
summary(tagged_stems$tag)
filter(tagged_stems, tag == 4194)

bay_laurel <- tagged_stems %>%
      select(plot, date, cluster, tag, slc, species, status, dbh1, dbh2, year_dbh1, year_dbh2, delta_dbh, dbh2_1_ratio, inst_grwth_rate, doy, year, location) %>% 
      filter(species == "UMCA")
str(bay_laurel)
summary(bay_laurel)
bay_laurel <- droplevels(bay_laurel)

## I found that one of the tag numbers has some kind of duplicate entries, where `dbh2` is either 2.0 or 2.1
filter(bay_laurel, tag == 4194)
#### Assign the entries for tag 4194 with slc == 60 the 2012 sampling date
bay_laurel$date[bay_laurel$tag == 4194 & bay_laurel$date == "2014-04-18" & bay_laurel$slc == 60] <- "2012-05-14"
#### Remove the records where dbh2 == 2.0
bay_laurel <- bay_laurel[!(bay_laurel$tag==4194 & bay_laurel$dbh2==2.0),]
summary(bay_laurel)
#### Make the same corrections to the `tagged_stems` data frame
tagged_stems$date[tagged_stems$tag == 4194 & tagged_stems$date == "2014-04-18" & tagged_stems$slc == 60] <- "2012-05-14"
tagged_stems <- tagged_stems[!(tagged_stems$tag==4194 & tagged_stems$dbh2==2.0),]# remove rows based on conditions of tag equaling 4194 and dbh2 equaling 2.0
write.csv(tagged_stems, "analysis/data/tagged_stems_growth.csv")# overwrite the incorrect data

#### Read corrected data set back in
tagged_stems <- read.csv("analysis/data/tagged_stems_growth.csv")
summary(tagged_stems)
str(tagged_stems)
tagged_stems$date <- as.Date(tagged_stems$date)

summary(bay_laurel)# there are 1075 `NA` values for `slc`
summary(filter(bay_laurel, is.na(slc)))# it appears that most are either dead or missing tagged_stems at the time of observation
summary(filter(bay_laurel, status == "Alive" & is.na(slc)))# 105 are alive with NA for slc
alive.0.slc <- (filter(tagged_stems, species == "UMCA" & status == "Alive" & is.na(slc)))
summary(alive.0.slc)
unique(alive.0.slc$tag)
alive.0.slc$notes# Most of these were due to a status change during a subsequent year during quality control efforts. A stem may have been "Dead" one year, but has a leaf count the next, indicating that it was not in fact dead, just really unhealthy.
rm(alive.0.slc)


summary(bay_laurel)
library(ggplot2)
qplot(doy, slc, data = bay_laurel %>% filter(status == "Alive", year != 2003), facets = year ~ .)
qplot(dbh2_1_ratio, slc, data = bay_laurel %>% filter(status == "Alive", year != 2003), facets = year ~ .)
qplot(dbh2_1_ratio, log(slc+1), data = bay_laurel %>% filter(status == "Alive", year != 2003), facets = year ~ .)
cor(bay_laurel$dbh2_1_ratio, log(bay_laurel$slc+1), use = "complete")
cor(bay_laurel$dbh2_1_ratio, bay_laurel$slc, use = "na.or.complete")
par(mfrow = c(1,2))
hist(bay_laurel$slc)
hist(log(bay_laurel$slc + 1))

## Create `pairs` plot of variables
#### Histogram & correlation functions from help documentation
panel.hist <- function(x, ...)
{
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(usr[1:2], 0, 1.5) )
      h <- hist(x, plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- abs(cor(x, y, use = "na.or.complete"))
      txt <- format(c(r, 0.123456789), digits = digits)[1]
      txt <- paste0(prefix, txt)
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(bay_laurel %>% filter(status == "Alive") %>% 
            select(slc, dbh1, dbh2, delta_dbh, dbh2_1_ratio, inst_grwth_rate, year, plot),
      upper.panel = panel.cor,
      lower.panel = panel.smooth,
      diag.panel = panel.hist,
      cex.labels = 2, font.labels = 2)

summarise(bay_laurel %>% filter(year == 2005), mdoy = max(doy))
filter(bay_laurel, year == 2005, doy == 146)


###### These may be useful functions for reference #### 
matplot(cbind(avg_rain_year,avg_slc_year),type="b")
matplot(cbind(rainy_days,avg_slc_year),type="b")

library(plotrix)
twoord.stackplot(lx=years, rx=years, ldata=avg_slc_year, rdata=avg_rain_year, lcol="red", rcol="blue", ltype="b", rtype="b", rylab="Mean SLC", lylab="Mean Rainfall", xlab="Year")

library(ggplot2)
qplot(years, avg_slc_year, geom = c("point", "line"))

ggplot() +
      geom_line(data = bay_dat1, aes(x = years, y = colMeans(cbind(amou2004, amou2005, amou2006, amou2007, amou2008, amou2009, amou2010, amou2011, amou2012)), colour="blue")) +
      geom_line(data=bay_dat1,aes(x=years,y=colMeans(cbind(cbind(avgslc04, avgslc05, avgslc06, avgslc07, avgslc08, avgslc09, avgslc10, avgslc11, avgslc12))), colour="red"))

