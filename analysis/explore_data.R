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

stems$Date <- as.character(stems$Date)
class(stems$Date)
library(lubridate)
stems$year <- year(as.Date(stems$Date, format = "%m/%d/%Y"))
summary(stems)
stems$year[is.na(stems$year)] <- 2012
stems[stems == -9999] <- NA
summary(stems)

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

library(ggplot2)
qplot(year, sum(dbh), data = stems)

dbh <- filter(stems, dbh>0, year == 2014)
length(unique(dbh$year))

length(unique(stems$tag))
# 4452 unique tag numbers in the data

dbh2014 <- stems %>% filter(dbh>0, year == 2014)
length(unique(dbh2014$tag))
# 4158, indicating there are 294 stems with no dbh measured in 2014

dbh2003 <- stems %>% select(plot, tag, dbh, year) %>% 
      filter(dbh>0, year == 2003) # these 2 are the same as in 2005
rm(dbh2003)

dbh2004 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2004, plot == "ANN02")
# 3109 observations

dbh2005 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2005)
# 712 observations
anti_join(dbh2004, dbh2005, by = "tag") # this indicates that the dbh data for 2004 and 2005 are unique
dbh1 <- union(dbh2004, dbh2005)

dbh2006 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2006)
summary(anti_join(dbh2006, dbh1, by = "tag")) # 2004, 2005, 2006 are independent sets
dbh1 <- union(dbh1, dbh2006)

dbh2007 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2007)
anti_join(dbh2007, dbh1, by = "tag") # there is one tag with remeasured dbh in 2007
inner_join(dbh2007, dbh1, by = "tag")
dbh1 <- union(dbh1, anti_join(dbh2007, dbh1, by = "tag"))
length(unique(dbh1$tag))

dbh2008 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2008)
anti_join(dbh2008, dbh1, by = "tag")
dbh1 <- union(dbh1, anti_join(dbh2008, dbh1, by = "tag"))

dbh2009 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2009)
anti_join(dbh2009, dbh1, by = "tag")
dbh1 <- union(dbh1, anti_join(dbh2009, dbh1, by = "tag"))

dbh2010 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2010)
anti_join(dbh2009, dbh1, by = "tag") # no unique dbh values for 2010, so no new recruitment, or stems entering the study

dbh2011 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2011)
anti_join(dbh2011, dbh1, by = "tag")
dbh1 <- union(dbh1, anti_join(dbh2011, dbh1, by = "tag"))

dbh2012 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2012)
anti_join(dbh2012, dbh1, by = "tag")
dbh1 <- union(dbh1, anti_join(dbh2012, dbh1, by = "tag"))

dbh2014 <- stems %>% select(plot, tag, dbh, year) %>%
      filter(dbh>0, year == 2014)
anti_join(dbh2014, dbh1, by = "tag") # this inidcates that there are 113 stems that were measured for the first time in 2014

library(ggplot2)
# Scatterplot of first dbh measurement (y-axis) against year
qplot(year, dbh, data = dbh1, geom="jitter")
# Scatterplot of most recent dbh measurement
qplot(year, dbh, data = dbh2014, geom = "jitter")

dbh <- right_join(dbh1, dbh2014, by = c("plot","tag"))
summary(dbh)
sum(dbh$dbh.x, na.rm = TRUE)
sum(dbh$dbh.y, na.rm = TRUE)
sum(dbh$dbh.y, na.rm = TRUE) - sum(dbh$dbh.x, na.rm = TRUE)
# [1] 5781.91 centimeters of total growth in diameter for all stems remeasured in 2014

dbh <- rename(dbh, dbh1 = dbh.x, dbh2 = dbh.y, year_dbh1 = year.x, year_dbh2 = year.y) # This data frame has all the stems for which we have an initial measurement AND a remeasurement in 2014
dbh$delta_dbh <- dbh$dbh2 - dbh$dbh1
write.csv(dbh, "analysis/data/dbh_remeasures.csv")
summary(dbh)
filter(dbh, delta_dbh < -20)

library(ggplot2)
qplot(dbh1, dbh2, data = dbh)
qplot(delta_dbh, dbh2, data = dbh)

#Join dbh measure, remeasure and change data to `stems` data frame
stems <- left_join(stems, dbh, by = c("tag","plot"))
str(stems)
summary(stems)
filter(stems, tag == 2285)
stems %>% select(plot, cluster, tag, species, delta_dbh, year, year_dbh1, status) %>% 
      group_by(plot, cluster, tag) %>% 
      filter(delta_dbh < -15, year == 2014, status == "Alive")
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

