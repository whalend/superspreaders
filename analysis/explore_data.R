--------
# Exploring UMCA stem infection through time #
# Superspreaders chapter of dissertation      
--------

######Set working directory, read in all stem data, examine data structure ####
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
