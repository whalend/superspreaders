## Data subset for Kerry Wininger, SSU

### "Would you be able to do me the favor of passing along any and all tree-level data you have from FOP, JROTH, and MROTH from last year?" *last year is referring to 2014

load("~/GitHub/superspreaders/.RData")
load("~/GitHub/superspreaders/stems_plots.Rdata")
str(stems)
str(plots)
str(dbh2014)

### Filter stem data to subset of requested plots and year
library(dplyr)

kw_stem_data <- stems %>%
      filter(plot %in% c("FOP04","FOP324","FOP333","FOP340","JROTH01","JROTH02","JROTH03","JROTH04","JROTH05","MROTH01","MROTH02","MROTH03"), year == 2014)
str(kw_stem_data)
summary(kw_stem_data)
write.csv(kw_stem_data, "deliverables/kw_stem_data.csv")

### "Plot temperature and precipitation points would be good too. Thanks for sending along some of this already." The temperature data is ready to go easily enough at the hourly level, but I still have some work to do on the precipitation data.

### Read in the long format temperature data, it's a big file so might take a little longer than usual
long_ppt <- readRDS("~/Documents/microclimate_tonini/Temperature/eof_pred_long.rds")
str(long_ppt)
kw_temp_data <- long_ppt %>% select(-posix_time) %>% filter(plotid %in% c("FOP04","FOP324","FOP333","FOP340","JROTH01","JROTH02","JROTH03","JROTH04","JROTH05","MROTH01","MROTH02","MROTH03"))
summary(kw_temp_data)
write.csv(kw_temp_data, "deliverables/kw_temp_data.csv")
