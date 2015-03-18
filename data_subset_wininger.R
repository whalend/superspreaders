## Data subset for Kerry Wininger, SSU

### "Would you be able to do me the favor of passing along any and all tree-level data you have from FOP, JROTH, and MROTH from last year?" *last year is referring to 2014

load("~/GitHub/superspreaders/.RData")
load("~/GitHub/superspreaders/stems_plots.Rdata")
str(stems)
str(plots)
str(dbh2014)

### Filter stem data to subset of requested plots and year
library(dplyr)

kw_data <- stems %>%
      filter(plot %in% c("FOP04","FOP324","FOP333","FOP340","JROTH01","JROTH02","JROTH03","JROTH04","JROTH05","MROTH01","MROTH02","MROTH03"), year == 2014)
str(kw_data)
summary(kw_data)

### "Plot temperature and precipitation points would be good too. Thanks for sending along some of this already." The temperature data is ready to go easily enough at the hourly level, but I still have some work to do on the precipitation data.




