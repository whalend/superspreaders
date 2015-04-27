## Explore data from vegetation transects and untagged host tree DBH ####

untagged_dbh <- read.csv("analysis/data/alt_species_dbh.csv")
str(alt_dbh)
summary(untagged_dbh)
untagged_dbh$date <- as.Date(untagged_dbh$Date, format = "%m/%d/%Y")

library(dplyr)
filter(untagged_dbh, is.na(date))
filter(untagged_dbh, PlotID == "SECRET02")
untagged_dbh$date[is.na(untagged_dbh$date)] <- "2012-05-22"

untagged_dbh <- as.tbl(untagged_dbh)
untagged_dbh <- untagged_dbh %>%
      rename(visit_id = PlotVisit_PlotVisitID, plot_id = PlotID) %>%
      select(-Date)


veg_data <- read.csv("analysis/data/veg_transects_qry.csv")
str(veg_data)
summary(veg_data)
