# Tagged stem DBH data for Jim to implement in FVS
### This uses the data where I corrected some of the negative growth values. There are still some with negative growth, but I had no way to provide reasonable replacement values. I haven't finished checking the positive growth outliers yet.

#### Read in data with corrected DBH ####
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

#### Subset and structure data
library(dplyr)
library(tidyr)
fvs_dbh <- stems %>% select(plot, year, tag, species, status, dbh) %>%
      filter(dbh > 0)
fvs_dbh <- as.tbl(fvs_dbh)
class(fvs_dbh)
fvs_dbh
fvs_dbh  %>% spread(species, dbh)
fvs_dbh[5391:5392,]
fvs_dbh[c(5413:5414,5424:5425),]
fvs_dbh2 <- fvs_dbh[-c(5414,5424:5425,5392),]
fvs_dbh2 <- fvs_dbh2 %>% spread(species, dbh)
write.csv(fvs_dbh2, "deliverables/fvs_tagged_dbh.csv")
