--------
# Exploring UMCA stem infection through time #
# Superspreaders chapter of dissertation      
--------

#### Read in all stem data, examine data structure ####
#This is only the biologically related field data

stems <- read.csv("analysis/data/stem_summary_query_2004-2014.txt")
head(stems)
summary(stems)
str(stems)


### Need to do some data cleaning and prep ####

## Create date variable as posix and integer/factor year variable ###
stems$Date <- as.character(stems$Date)
class(stems$Date)
library(lubridate)
stems$year <- year(as.Date(stems$Date, format = "%m/%d/%Y"))
summary(stems)
## Why are there NA's in the `year` variable?
library(plyr)
library(dplyr)
stems_sub <- filter(stems, is.na(year))
stems_sub2 <- filter(stems, PlotID == "SECRET02", SpeciesID == "QUAG", year == 2012)
# So the 19 NA's for in the table are from the SECRET02 plot and belong to 2012. I am not sure why these did not query correctly in the database because the data is entered with the proper date...
# I replaced the NA's for the `year` variable with 2012
stems$year[is.na(stems$year)] <- 2012
summary(stems)
rm(stems_sub)
rm(stems_sub2)

## Currently unused code for string substitution and date calculation
# stems$date <- gsub("/", "-",stems$Date)
# stems$date <- dmy_hms(stems$date)
# stems$year <- year(stems$date)

## Replace all `-9999` codes with NA ###
stems[stems == -9999] <- NA
summary(stems)


#### Rename columns using `rename` funcion from `dplyr` ####
library(dplyr)
stems <- rename(stems, plot = PlotID)
stems <- rename(stems, cluster = ClusterID, tag = TagNumber, 
                species = SpeciesID, status = StemStatus, 
                alive_class = AliveClass, dead_class= DeadClass,
                slc = SympLeafCount, canker = CankerPresent, dbh = DBH,
                lide_lf_symp = LideFoliarSymp, sod_dead = SOD_Dead, 
                location = Location)
str(stems)

write.csv(stems,"analysis/data/all_stems.csv") # Write the revised stem data to CSV

#### Further exploration of revised stem data ####
stems <- read.csv("analysis/data/all_stems.csv") # Read in if necessary
str(stems)
summary(stems)
summary(stems$species) # 5 host species listed, so this appears legit
summary(stems$status) # There are a lot of different `X` statuses, which I am not particularly concerned about. These will all just be treated as missing. However, there are 26 stems that appear to have a `blank` status type. 
class(stems$status)
View(filter(stems, status == ""))
## A good number of these stems needed an `X` status, but there were at least 8 that needed correcting in the database. The primary cause for the blank status appeared to be a duplication of the visit date for a single stem. Once one of these records was deleted the data appeared. So, I am going to export the query from the data base again and reassess the data with these corrections.
stems <- read.csv("analysis/data/stem_summary_qry2_04-14.csv")
summary(stems)

stems$Date <- as.Date(stems$Date, format = "%m/%d/%Y")
class(stems$Date)
library(lubridate)
stems$year <- year(as.Date(stems$Date, format = "%m/%d/%Y"))
summary(stems)
stems$year[is.na(stems$year)] <- 2012
stems[stems == -9999] <- NA
summary(stems)


library(dplyr)
stems <- rename(stems, plot = PlotID)
stems <- rename(stems, cluster = ClusterID, tag = TagNumber, 
                species = SpeciesID, status = StemStatus, 
                alive_class = AliveClass, dead_class= DeadClass,
                slc = SympLeafCount, canker = CankerPresent, dbh = DBH,
                lide_lf_symp = LideFoliarSymp, sod_dead = SOD_Dead, 
                location = Location, notes = StemSymptomNotes)
str(stems)
summary(stems)
summary(stems$species) # 5 host species listed, so this appears legit
summary(stems$status) # Now there are no blanks, 36864 Alive records, 2482 dead records, 1 "Nearly dead" so changing that to "Alive"
stems$status[stems$status=="Nearly dead"] <- "Alive"
stems <- droplevels(stems)# drop unused factor levels

stems <- as.tbl(stems) # coerce to dplyr style data table
class(stems)
length(unique(stems$tag))
# 4452 unique tag numbers in the data
write.csv(stems, "analysis/data/all_stems_corrected.csv")

## Subset DBH data to create a data frame of remeasured stems ####
dbh2014 <- stems %>% filter(dbh>0, year == 2014)
length(unique(dbh2014$tag))
# 4158, indicating there are 294 stems with no dbh measured in 2014

dbh2003 <- stems %>% select(plot, tag, dbh, year) %>% 
      filter(dbh>0, year == 2003) # these 2 are the same as in 2005
# 11 observations

dbh2004 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2004)
# 3109 observations
anti_join(dbh2003, dbh2004, by = "tag")
dbh_a <- union(dbh2004, dbh2003)


dbh2005 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2005)
# 712 observations
anti_join(dbh2005, dbh_a, by = "tag") # this indicates that the dbh data for 2004 and 2005 are unique
dbh_a <- union(dbh_a, anti_join(dbh2005, dbh_a, by = "tag"))

dbh2006 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2006)
anti_join(dbh2006, dbh_a, by = "tag") # 2004, 2005, 2006 are independent sets
dbh_a <- union(dbh_a, anti_join(dbh2006, dbh_a, by = "tag"))

dbh2007 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2007)
anti_join(dbh2007, dbh_a, by = "tag") # there is one tag with remeasured dbh in 2007
dbh_a <- union(dbh_a, anti_join(dbh2007, dbh_a, by = "tag"))
length(unique(dbh_a$tag)) # check for repeated tags, if obs > value, then there is a duplicate

dbh2008 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2008)
anti_join(dbh2008, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2008, dbh_a, by = "tag"))

dbh2009 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2009)
anti_join(dbh2009, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2009, dbh_a, by = "tag"))

dbh2010 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2010)
anti_join(dbh2009, dbh_a, by = "tag") # no unique dbh values for 2010, so no new recruitment, or stems entering the study

dbh2011 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2011)
anti_join(dbh2011, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2011, dbh_a, by = "tag"))

dbh2012 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2012)
anti_join(dbh2012, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2012, dbh_a, by = "tag"))

dbh2014 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2014)
anti_join(dbh2014, dbh_a, by = "tag") # this inidcates that there are 113 stems that were measured for the first time in 2014

#### Explore DBH data: negative change in DBH ####
library(ggplot2)
# Scatterplot of first dbh measurement (y-axis) against year
qplot(year, dbh, data = dbh_a, geom="jitter")
# Scatterplot of most recent dbh measurement
qplot(year, dbh, data = dbh2014, geom = "jitter")

# Data frame of stems that were measured or remeasured in 2014
dbh <- right_join(dbh_a, dbh2014, by = c("plot","tag"))
summary(dbh)
sum(dbh$dbh.x, na.rm = TRUE)
sum(dbh$dbh.y, na.rm = TRUE)
sum(dbh$dbh.y, na.rm = TRUE) - sum(dbh$dbh.x, na.rm = TRUE)
# [1] 6345.03 centimeters of total growth in diameter for all stems remeasured in 2014

dbh <- rename(dbh, dbh1 = dbh.x, dbh2 = dbh.y, year_dbh1 = year.x, year_dbh2 = year.y) # This data frame has all the stems for which we have an initial measurement AND a remeasurement in 2014
dbh$delta_dbh <- dbh$dbh2 - dbh$dbh1
write.csv(dbh, "analysis/data/dbh_remeasures.csv")
summary(dbh)
# filter(dbh, delta_dbh < -20)

library(ggplot2)
qplot(dbh1, dbh2, data = dbh)
qplot(delta_dbh, dbh2, data = dbh, color = year_dbh1)

#Join dbh measure, remeasure and change data to `stems` data frame ####
dbh <- read.csv("analysis/data/dbh_2014_remeasures.csv")
stems <- left_join(stems, dbh, by = c("tag","plot"))
str(stems)
summary(stems)
qplot(delta_dbh, dbh2, data = stems %>% select(plot, cluster, tag, species, delta_dbh, dbh2, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead"), 
      color = status)

qplot(year, dbh, data = stems %>% filter(status == "Alive" | status == "Dead"), color = status, geom = "jitter")


stems %>% select(plot, cluster, tag, species, delta_dbh, year, year_dbh1, status) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(delta_dbh < -5, year == 2014, status == "Alive")
# I changed values for 9 of these 14 stems

stems %>% select(plot, cluster, tag, species, delta_dbh, year, year_dbh1, status) %>% group_by(plot, cluster, tag) %>% 
      filter(between(delta_dbh, -15, -10), year == 2014, status == "Alive")
# I changed values for 7 of 9 stems

stems %>% select(plot, cluster, tag, species, delta_dbh, year, year_dbh1, status) %>% group_by(plot, cluster, tag) %>% 
      filter(between(delta_dbh, -10, -7), year == 2014, status == "Alive")
# I changed values for 11 of 20 stems

stems %>% select(plot, cluster, tag, species, delta_dbh, year, year_dbh1, status) %>% group_by(plot, cluster, tag) %>% 
      filter(between(delta_dbh, -6.9, -5), year == 2014, status == "Alive")
# I changed values for 6 of 18 stems


#### Create a plot-level data frame of infection status ####
stems <- read.csv("analysis/data/all_stems_corrected.csv")
library(dplyr)

plots_umca <- stems %>%
      select(plot, date, species, slc, year, status) %>%
      group_by(plot, year, species) %>%
      filter(species == "UMCA", status == "Alive") %>%
      summarise(uninfected_bay_ct = length(which(slc==0)), infected_bay_ct = length(which(slc > 0)), tot_bay = length(species)) 

# plots_umca$infected <- ifelse(plots_umca$infected_bay_ct==0, 0, 1)
summary(plots_umca)
plots_umca <- droplevels(plots_umca)


length(which(plots_umca$ct_bay_NA > 0))
filter(plots_umca, infected == 0)
filter(plots_umca, infected_bay_ct == 0)

write.csv(plots_umca, "analysis/data/plots_umca_infection.csv")
# I sent this file to Francesco for informing the spread model


plot_qusp <- stems %>%
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

#### Checking large positive DBH changes ####
## So at this point I think I have corrected all the negative dbh change values that I legitimately could in the database. I have gone back and exported the query again and use this file below, and then recycling code from the beginning of this document.
stems <- read.csv("analysis/data/Stem_Summary_negDBH_corrected_201402.csv")
summary(stems)
str(stems)

stems$Date <- as.character(stems$Date)
class(stems$Date)
stems$date <- as.Date(stems$Date, format = "%m/%d/%Y")
library(lubridate)
stems$year <- year(stems$date)
summary(stems)
library(plyr)
library(dplyr)
filter(stems, is.na(year))
filter(stems, PlotID == "SECRET02", SpeciesID == "QUAG")

stems$year[is.na(stems$year)] <- 2012
filter(stems, is.na(date))
stems$date[is.na(stems$date)] <- "2012-05-22"
summary(stems)
stems[stems == -9999] <- NA
summary(stems)

stems <- rename(stems, plot = PlotID)
stems <- rename(stems, cluster = HostID, tag = TagNumber, 
                species = SpeciesID, status = StemStatus, 
                alive_class = AliveClass, dead_class= DeadClass,
                slc = SympLeafCount, canker = CankerPresent, dbh = DBH,
                lide_lf_symp = LideFoliarSymp, sod_dead = SOD_Dead, 
                location = Location)
str(stems)
stems <- rename(stems, notes = StemSymptomNotes)
summary(stems$status)
stems$status[stems$status=="Nearly dead"] <- "Alive"
stems <- droplevels(stems)
class(stems)
stems <- as.tbl(stems)

#### Subset DBH data to create a data frame of remeasured stems: Round 2 ####

detach("package:lubridate", unload=TRUE)# b/c it has `union` function masking the one from `dplyr` that apparently behaves a little differently, creating a list instead of a data frame

dbh2014 <- stems %>% filter(dbh>0, year == 2014)
length(unique(dbh2014$tag))
# 4158, indicating there are 294 stems with no dbh measured in 2014

dbh2003 <- stems %>% select(plot, tag, dbh, year, date) %>% 
      filter(dbh>0, year == 2003) # these 2 are the same as in 2005
# 11 observations

dbh2004 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2004)
# 3109 observations
anti_join(dbh2003, dbh2004, by = "tag")
dbh_a <- union(dbh2004, dbh2003)

dbh2005 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2005)
# 704 observations
anti_join(dbh2005, dbh_a, by = "tag") # this indicates that the dbh data for 2004 and 2005 are unique
dbh_a <- union(dbh_a, anti_join(dbh2005, dbh_a, by = "tag"))

dbh2006 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2006)
anti_join(dbh2006, dbh_a, by = "tag") # 2004, 2005, 2006 are independent sets
dbh_a <- union(dbh_a, anti_join(dbh2006, dbh_a, by = "tag"))

dbh2007 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2007)
anti_join(dbh2007, dbh_a, by = "tag") # there is one tag with remeasured dbh in 2007
dbh_a <- union(dbh_a, anti_join(dbh2007, dbh_a, by = "tag"))
length(unique(dbh_a$tag)) # check for repeated tags, if obs > value, then there is a duplicate

dbh2008 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2008)
anti_join(dbh2008, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2008, dbh_a, by = "tag"))

dbh2009 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2009)
anti_join(dbh2009, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2009, dbh_a, by = "tag"))

dbh2010 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2010)
anti_join(dbh2009, dbh_a, by = "tag") # no unique dbh values for 2010, so no new recruitment, or stems entering the study

dbh2011 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2011)
anti_join(dbh2011, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2011, dbh_a, by = "tag"))

dbh2012 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2012)
anti_join(dbh2012, dbh_a, by = "tag")
dbh_a <- union(dbh_a, anti_join(dbh2012, dbh_a, by = "tag"))

dbh2014 <- stems %>% select(plot, tag, dbh, year, date) %>%
      filter(dbh>0, year == 2014)
anti_join(dbh2014, dbh_a, by = "tag") # this inidcates that there are 113 stems that were measured for the first time in 2014
anti_join(dbh_a, dbh2014, by = "tag")# this indicates that there were 284 stems previously measured that were not part of the remeasurement in 2014. This looks in part to be due to plot abandonments/decomissions, but may also be due to other reasons.

## Some plots of DBH measurements ####
library(ggplot2)
# Scatterplot of first dbh measurement (y-axis) against year
qplot(tag, dbh, data = dbh_a, geom="jitter", facets = year ~.)

# Scatterplot of most recent dbh measurement
qplot(tag, dbh, data = dbh2014, geom = "jitter", facets = year ~.)
qplot(dbh, data = dbh2014)

# Data frame of stems that were measured or remeasured in 2014
dbh <- right_join(dbh_a, dbh2014, by = c("plot","tag"))
summary(dbh)

dbh <- rename(dbh, dbh1 = dbh.x, dbh2 = dbh.y, year_dbh1 = year.x, year_dbh2 = year.y) # This data frame has all the stems for which we have an initial measurement AND a remeasurement in 2014
dbh$delta_dbh <- dbh$dbh2 - dbh$dbh1
write.csv(dbh, "analysis/data/dbh_2014_remeasures.csv")
summary(dbh)
sum(dbh$delta_dbh, na.rm=T)
# filter(dbh, delta_dbh < -20)


#### Check the positive DBH change outliers of stems remeasured in 2014 ####
dbh <- read.csv("analysis/data/dbh_2014_remeasures.csv")
summary(dbh)
library(ggplot2)
qplot(dbh1, dbh2, data = dbh)
qplot(delta_dbh, dbh2, data = dbh, color = year_dbh1)
summary(dbh)
dbh <- as.tbl(dbh)

# Join dbh measure, remeasure, and change data to `stems` data frame
dbh <- select(dbh, -X) %>%
      rename(date_dbh1 = date.x, date_dbh2 = date.y)
summary(dbh)
summary(stems)
stems <- left_join(stems, dbh, by = c("tag","plot"))
str(stems)
summary(stems)

# Exploratory plot of dbh data
qplot(delta_dbh, dbh1, data = stems %>% select(plotid, cluster, tag, species, delta_dbh, dbh1, year, status) %>% 
            group_by(plotid, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead", species == "UMCA"), 
      color = status)
# A visual check shows that there are a handful of stems with growth > 20cm

#### Suggestions from Margaret about exploring DBH changes ####
# I emailed Richard and Margaret about thresholds for where I should really look more closely at the measurements. Richard noted that there isn't much out there on growth rates for these species. I did find something done for tree species of the Northeastern U.S. <http://www.fs.fed.us/ne/newtown_square/publications/research_papers/pdfs/scanned/OCR/ne_rp649.pdf>

# Margaret suggested a few things, including other ways of looking at the data:
# 1. "Have you tried looking at relative growth instead of the absolute increment growth?  Either DBH2/DBH1, or the instantaneous rate ln(DBH2/DBH1) / (time2-time1).  20 cm growth means a lot more to a 5 cm tree than a 50 cm tree… either way, it's still big.  I'd also consider this species by species, because I bet some grow much faster than others."

stems$dbh2_1_ratio <- stems$dbh2/stems$dbh1# absolute increment growth
stems$inst_grwth_rate <- log(stems$dbh2/stems$dbh1) / (stems$year_dbh2 - stems$year_dbh1)
#stems$time_diff <- (stems$year_dbh2 - stems$year_dbh1)
summary(stems)
library(lubridate)
stems$doy <- yday(stems$date)# create a day of year variable
write.csv(stems, "analysis/data/stems_growth.csv")
summary(stems$status)
summary(stems$tag)

qplot(dbh2_1_ratio, data = stems %>% select(plotid, tag, species, dbh2_1_ratio, dbh1, year, status, doy) %>% 
            group_by(plotid, tag) %>% 
            filter(status == "Alive" | status == "Dead", dbh1 >= 2.0),
      color = status)

qplot(stems$dbh2_1_ratio, main = "Relative Growth")

qplot(dbh1, dbh2_1_ratio, data = stems %>% select(plotid, cluster, tag, species, dbh2_1_ratio, dbh1, year, status, dbh) %>% 
            group_by(plotid, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead", dbh1 >= 2.0), 
      color = status, main = "Change in DBH vs. Original DBH")

qplot(dbh1, inst_grwth_rate, data = stems %>% select(plotid, cluster, tag, species, inst_grwth_rate, dbh1, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead", dbh1 >= 2.0), 
      color = status, main = "Instantaneous Growth Rate vs. Original DBH")

qplot(dbh1, inst_grwth_rate, data = stems %>% select(plotid, cluster, tag, species, inst_grwth_rate, dbh1, year_dbh1, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(species == "QUAG", year == 2014, status == "Dead" | status == "Alive", dbh1 >= 2.0), 
      color = status, main = "QUAG Instantaneous Growth Rate vs Original DBH")

qplot(dbh1, inst_grwth_rate, data = stems %>% select(plotid, cluster, tag, species, inst_grwth_rate, dbh1, year_dbh1, year, status) %>% 
            group_by(plotid, cluster, tag) %>%
            filter(species == "UMCA", year == 2014, status == "Dead" | status == "Alive", dbh1 >= 2.0), 
      color = status, main = "UMCA Instantaneous Growth Rate vs Original DBH")

boxplot(dbh ~ year, data = stems %>% filter(species == "UMCA", status == "Alive", year == 2005|year==2012|year == 2014))

stems %>% filter(species == "UMCA", status == "Alive", year == 2005|year==2012|year == 2014)

u.2005 <- stems %>%
      filter(species == "UMCA", status == "Alive", year == 2005|year==2014)

t.test(dbh~year, data = stems %>% filter(species == "UMCA", status == "Alive", year_dbh1 == 2005|year_dbh2 == 2014))
t.test(dbh~year, data = stems %>% filter(species == "QUAG", status == "Alive", year == 2005|year == 2014))
t.test(dbh~year, data = stems %>% filter(species == "QUAG", status == "Dead", year == 2005|year == 2014))
t.test(dbh~year, data = stems %>% filter(species == "QUKE", status == "Alive", year == 2005|year == 2014))
t.test(dbh~year, data = stems %>% filter(species=="QUAG"|species == "QUKE", status == "Alive", year == 2005|year == 2014))



stems %>% select(plot, cluster, tag, species, year, year_dbh1, status, delta_dbh, dbh1, dbh2, notes) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(delta_dbh > 20, year == 2014, status == "Alive")
# 8 stems

stems %>% select(plot, tag, cluster, species, year, year_dbh1, status, dbh, delta_dbh, dbh1, dbh2, notes) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(plot == "LAZY05", species == "QUKE", status == "Alive")

stems %>% select(plot, tag, cluster, species, year, year_dbh1, status, dbh, delta_dbh, dbh1, dbh2, notes) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(plot == "BOWES02", species == "UMCA", status == "Alive")

stems %>% arrange(desc(delta_dbh)) %>% 
      select(plot, cluster, tag, species, delta_dbh, dbh1, dbh2, year, year_dbh1, status) %>% group_by(plot, cluster, tag) %>% 
      filter(between(delta_dbh, 10, 20), year == 2014, status == "Alive")
# 39 stems
summary(stems)
rm(list = ls(pattern = "dbh20*"))

## Exploring Symptomatic Leafcount Data ####
stems <- read.csv("analysis/data/stems_growth.csv")
str(stems)
summary(stems)
stems$date <- as.Date(stems$date)
library(plyr);library(dplyr)
stems <- as.tbl(stems)
summary(stems)
stems <- select(stems, -X, -Date)
summary(stems$date)
summary(stems$tag)
filter(stems, tag == 4194)

bay_laurel <- stems %>%
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
#### Make the same corrections to the `stems` data frame
stems$date[stems$tag == 4194 & stems$date == "2014-04-18" & stems$slc == 60] <- "2012-05-14"
stems <- stems[!(stems$tag==4194 & stems$dbh2==2.0),]# remove rows based on conditions of tag equaling 4194 and dbh2 equaling 2.0
write.csv(stems, "analysis/data/stems_growth.csv")# overwrite the incorrect data

#### Read corrected data set back in
stems <- read.csv("analysis/data/stems_growth.csv")
summary(stems)
str(stems)
stems$date <- as.Date(stems$date)

summary(bay_laurel)# there are 1075 `NA` values for `slc`
summary(filter(bay_laurel, is.na(slc)))# it appears that most are either dead or missing stems at the time of observation
summary(filter(bay_laurel, status == "Alive" & is.na(slc)))# 105 are alive with NA for slc
alive.0.slc <- (filter(stems, species == "UMCA" & status == "Alive" & is.na(slc)))
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

