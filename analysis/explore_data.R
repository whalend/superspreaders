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
summary(stems$status) # Now there are no blanks, 36864 Alive records, 2482 dead records

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

#### Explore DBH data: change in DBH ####
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

#Join dbh measure, remeasure and change data to `stems` data frame
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
      select(plot, Date, species, slc, year, status) %>%
      group_by(plot, year) %>%
      filter(species == "UMCA", status == "Alive") %>%
      summarise(uninfected_bay_ct = length(which(slc==0)), infected_bay_ct = length(which(slc > 0)), tot_bay = length(species)) 

plots_umca$infected <- ifelse(plots_umca$infected_bay_ct==0, 0, 1)
summary(plots_umca)

length(which(plots_umca$ct_bay_NA > 0))
filter(plots_umca, infected == 0)
filter(plots_umca, infected_bay_ct == 0)

write.csv(plots_umca, "analysis/data/plots_umca_infection.csv")
# I sent this file to Francesco for informing the spread model

#### Checking large positive DBH changes ####
## So at this point I think I have corrected all the negative dbh change values that I legitimately could in the database. I have gone back and exported the query again and use this file below, recycling code from the beginning of this document.
stems <- read.csv("analysis/data/Stem_Summary_negDBH_corrected_201402.csv")
summary(stems)

stems$Date <- as.character(stems$Date)
class(stems$Date)
library(lubridate)
stems$year <- year(as.Date(stems$Date, format = "%m/%d/%Y"))
summary(stems)
library(plyr)
library(dplyr)
stems_sub <- filter(stems, is.na(year))
stems_sub2 <- filter(stems, PlotID == "SECRET02", SpeciesID == "QUAG", year == 2012)

stems$year[is.na(stems$year)] <- 2012
summary(stems)
stems[stems == -9999] <- NA
summary(stems)

library(dplyr)
stems <- rename(stems, plot = PlotID)
stems <- rename(stems, cluster = ClusterID, tag = TagNumber, 
                species = SpeciesID, status = StemStatus, 
                alive_class = AliveClass, dead_class= DeadClass,
                slc = SympLeafCount, canker = CankerPresent, dbh = DBH,
                lide_lf_symp = LideFoliarSymp, sod_dead = SOD_Dead, 
                location = Location)
str(stems)
stems <- rename(stems, notes = StemSymptomNotes)
summary(stems)


#### Subset DBH data to create a data frame of remeasured stems: Round 2 ####
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
# 704 observations
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

library(ggplot2)
# Scatterplot of first dbh measurement (y-axis) against year
qplot(year, dbh, data = dbh_a, geom="jitter")
# Scatterplot of most recent dbh measurement
qplot(year, dbh, data = dbh2014, geom = "jitter")

# Data frame of stems that were measured or remeasured in 2014
dbh <- right_join(dbh_a, dbh2014, by = c("plot","tag"))
summary(dbh)

dbh <- rename(dbh, dbh1 = dbh.x, dbh2 = dbh.y, year_dbh1 = year.x, year_dbh2 = year.y) # This data frame has all the stems for which we have an initial measurement AND a remeasurement in 2014
dbh$delta_dbh <- dbh$dbh2 - dbh$dbh1
write.csv(dbh, "analysis/data/dbh_remeasures.csv")
summary(dbh)
sum(dbh$delta_dbh, na.rm=T)
# filter(dbh, delta_dbh < -20)


#### Check the positive DBH change outliers ####

library(ggplot2)
qplot(dbh1, dbh2, data = dbh)
qplot(delta_dbh, dbh2, data = dbh, color = year_dbh1)

# Join dbh measure, remeasure, and change data to `stems` data frame
dbh <- read.csv("analysis/data/dbh_remeasures.csv")
dbh <- select(dbh, -X)
stems <- left_join(stems, dbh, by = c("tag","plot"))
str(stems)
summary(stems)

# Exploratory plot of dbh data
qplot(delta_dbh, dbh1, data = stems %>% select(plot, cluster, tag, species, delta_dbh, dbh1, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead"), 
      color = status)
# A visual check shows that there are a handful of stems with growth > 20cm
# I emailed Richard and Margaret about thresholds for where I should really look more closely at the measurements. Richard noted that there isn't much out there on growth rates for these species. Margaret suggested a few things, including other ways of looking at the data:
# 1. "Have you tried looking at relative growth instead of the absolute increment growth?  Either DBH2/DBH1, or the instantaneous rate ln(DBH2/DBH1) / (time2-time1).  20 cm growth means a lot more to a 5 cm tree than a 50 cm tree… either way, it's still big.  I'd also consider this species by species, because I bet some grow much faster than others."
stems$dbh2_dbh1 <- stems$dbh2/stems$dbh1
stems$inst_grwth_rate <- log(stems$dbh2/stems$dbh1) / (stems$year_dbh2 - stems$year_dbh1)
stems$time_diff <- (stems$year_dbh2 - stems$year_dbh1)
summary(stems)

qplot(stems$dbh2_dbh1, main = "Relative Growth")

qplot(dbh1, dbh2_dbh1, data = stems %>% select(plot, cluster, tag, species, dbh2_dbh1, dbh1, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead"), 
      color = status, main = "Change in DBH vs. original DBH")

qplot(dbh1, inst_grwth_rate, data = stems %>% select(plot, cluster, tag, species, inst_grwth_rate, dbh1, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(year == 2014, status == "Alive" | status == "Dead"), 
      color = status, main = "Instantaneous Growth Rate vs. Original DBH")

qplot(dbh1, inst_grwth_rate, data = stems %>% select(plot, cluster, tag, species, inst_grwth_rate, dbh1, year_dbh1, year, status) %>% 
            group_by(plot, cluster, tag) %>%
            filter(species == "QUAG", year == 2014, status == "Dead" | status == "Alive"), 
      color = status, main = "QUAG Instantaneous Growth Rate vs Original DBH")


stems %>% select(plot, cluster, tag, species, year, year_dbh1, status, delta_dbh, dbh1, dbh2, notes) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(delta_dbh > 20, year == 2014, status == "Alive")
# 8 stems, 

stems %>% select(plot, tag, cluster, species, year, year_dbh1, status, dbh, delta_dbh, dbh1, dbh2, notes) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(plot == "LAZY05", species == "QUKE", status == "Alive")

stems %>% select(plot, tag, cluster, species, year, year_dbh1, status, dbh, delta_dbh, dbh1, dbh2, notes) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(plot == "BOWES02", species == "UMCA", status == "Alive")

stems %>% select(plot, cluster, tag, species, delta_dbh, dbh1, dbh2, year, year_dbh1, status) %>% group_by(plot, cluster, tag) %>% 
      filter(between(delta_dbh, 10, 20), year == 2014, status == "Alive")
# 39 stems



## Exploring Symptomatic Leafcount Data ####
str(stems)
summary(stems$status)
bay_slc <- stems %>%
      select(plot, Date, cluster, tag, species, status, dbh, year, slc) %>% 
      filter(species == "UMCA" & status == "Alive")
str(bay_slc)
bay_slc$Date <- as.Date(bay_slc$Date, format = "%m/%d/%Y")
summary(bay_slc)
bay_slc <- droplevels(bay_slc)

library(ggplot2)
qplot(slc, tag, data = bay_slc, facet = year ~ .)


### ### ### ### ### ### ### ### ### ### ### ### ### ### #
### Older code that can and will likely be abandoned ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### # 
library(foreign)
plot_dat1<-read.dbf("analysis/data/UMCA_slc_dbh_04-12.dbf")
str(plot_dat1)
summary(plot_dat1)
head(plot_dat1)

plot_dat1$avg_slc_04_12<-with(plot_dat1, rowMeans(cbind(avgslc04,avgslc05,avgslc06,avgslc07,avgslc08,avgslc09,avgslc10,avgslc11,avgslc12)))

years<-c("2004","2005","2006","2007","2008","2009","2010","2011","2012")
with(plot_dat1, boxplot(avgslc04,avgslc05,avgslc06,avgslc07,avgslc08,
                        avgslc09,avgslc10,avgslc11,avgslc12,names=years,xlab="Year",ylab="Mean SLC"))


library(reshape)
plot_dat1<-rename(plot_dat1,c(PLOT_ID="plotid"))
library(dplyr)
plot_dat1<-left_join(plot_dat1,plot_202_rain,by="plotid")

avg_slc_year<-with(plot_dat1, colMeans(cbind(avgslc04,avgslc05,avgslc06,avgslc07,avgslc08,
                                             avgslc09,avgslc10,avgslc11,avgslc12)))


plot(years,avg_slc_year,xlab="Year",ylab="Mean SLC")

avg_rain_year<-with(plot_202_rain, colMeans(cbind(amou2004,amou2005,amou2006,amou2007,amou2008,
                                                  amou2009,amou2010,amou2011,amou2012)))
rainy_days<-with(plot_202_rain, colMeans(cbind(freq2004,freq2005,freq2006,freq2007,freq2008,
                                               freq2009,freq2010,freq2011,freq2012)))

df_rain_slc<-as.data.frame(rbind(avg_rain_year,avg_slc_year))

###### These may be useful functions for reference #### 
matplot(cbind(avg_rain_year,avg_slc_year),type="b")
matplot(cbind(rainy_days,avg_slc_year),type="b")



library(plotrix)
twoord.stackplot(lx=years,rx=years,ldata=avg_slc_year,rdata=avg_rain_year,lcol="red",rcol="blue",
                 ltype="b",rtype="b",rylab="Mean SLC",lylab="Mean Rainfall",xlab="Year")

library(ggplot2)
qplot(years,avg_slc_year,geom=c("point","line"))

ggplot() +
      geom_line(data=bay_dat1,aes(x=years,y=colMeans(cbind(amou2004,amou2005,amou2006,amou2007,amou2008,
                                                           amou2009,amou2010,amou2011,amou2012)),
                                  colour="blue")) +
      geom_line(data=bay_dat1,aes(x=years,y=colMeans(cbind(cbind(avgslc04,avgslc05,avgslc06,avgslc07,avgslc08,
                                                                 avgslc09,avgslc10,avgslc11,avgslc12))),
                                  colour="red"))

