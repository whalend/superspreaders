# Pathogen ~ Biomass relationship for parameterizing LANDIS-II extension
# From Francesco's email: "if we come up with a decent relationship to use biomass related to pathogen production we can improve the current model"
load("~/GitHub/superspreaders/data_201508.RData")

summary(bay_laurel)
summary(stems_dbh)
summary(stems)

library(ggplot2); library(dplyr)

qplot(dbh1, data = bay_laurel %>% filter(status == "Alive", dbh1 >= 2)) +
      geom_vline(xintercept = c(10.2, 13.34))

qplot(dbh2, data = bay_laurel %>% filter(status == "Alive", dbh2 >= 2)) +
      geom_vline(xintercept = c(11.1, 14.64), colour = "red", linetype = "longdash")
qplot(log(dbh2), data = bay_laurel %>% filter(status == "Alive", dbh2 >= 2))

qplot(slc, data = bay_laurel %>% filter(status == "Alive", dbh2 >= 2, year > 2003))
qplot(log1p(slc), data = bay_laurel %>% filter(status == "Alive", dbh2 >= 2, year > 2003))

p1 <- ggplot(filter(bay_laurel, status == "Alive", dbh2 > 2, year > 2003), 
             aes(x = log(dbh2), y = log1p(slc))) +
      geom_point() +
      facet_grid(year ~ .)

quantile(filter(bay_laurel, dbh2 >=2)$slc, na.rm = T)
quantile(filter(bay_laurel, dbh2 >= 2)$dbh2, 0.9, na.rm = T)

log1p(filter(bay_laurel, dbh2 >=2)$slc, na.rm = T)
