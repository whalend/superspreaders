#' # Path Modeling of Climatic Influence on Spillover
#' I'm using the `piecewiseSEM package, which operationalized the methods from Shipley 2009 & 2013 in particular for fitting and assessing path models that enable applying hierarchical modeling methods.
#' 
#+ load data and packages ####
load("pathmodel_data.RData")
library(plyr); library(dplyr); library(nlme); library(lme4)
library(piecewiseSEM)

#+ panel functions for pairs plots
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      # correlation coefficient
      r <- cor(x, y)
      txt <- format(c(r, 0.123456789), digits = digits)[1]
      txt <- paste("r= ", txt, sep = "")
      text(0.5, 0.6, txt)
      
      # p-value calculation
      p <- cor.test(x, y)$p.value
      txt2 <- format(c(p, 0.123456789), digits = digits)[1]
      txt2 <- paste("p= ", txt2, sep = "")
      if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
      text(0.5, 0.4, txt2)
}

panel.hist <- function(x, ...)
{
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(usr[1:2], 0, 1.5) )
      h <- hist(x, plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

#+ examine oak plots data ####
summary(oak_sod)
oak_sod$plotid <- as.factor(oak_sod$plotid)
str(oak_sod)

## join wet-hours data
summary(wet_hrs)# since these are counts the NAs are in fact zeroes
wet_hrs[is.na(wet_hrs)] <- 0# introduce a lot of zeroes for some columns
oak_sod <- left_join(
      oak_sod, select(wet_hrs, -X), by = c("plotid","sample_year"))
oak_sod$plotid <- as.factor(oak_sod$plotid)

#+ pairs plots of oak plots data ####
pairs(na.omit(select(oak_sod, starts_with("rain"), ends_with("wet"))), 
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Oak Plots Rainfall Variables")
# The regression variables are the least correlated within totals & days.
# The 2D/3D interpolations are strongly correlated with each other, and strongly correlated with with the Voronoi polygon extractions.
# The days & total precipitation variables are fairly strongly correlated.
# Total wet hours & hours below 10 & 14 are strongly correlated with rainy days.
# Wet hours below 14 & total hours are correlated with total rainfall r = 0.66 - 0.76
# Wet hours below 10 is correlated with total rainfall r = 0.61 - 0.68

pairs(select(oak_sod, ends_with("rs"), ends_with("ds"), ends_with("wet")),
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Oak Plots Temperature Variables")
# Overall, correlations among variables are relatively weak with a few strong pairs. For specific pairs of variables it may be okay to use one from the dry season and one from the rainy season.
pairs(select(oak_sod, ends_with("ds"), ends_with("wet")),
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Oak Plots Temperature Variables")
# 

pairs(select(oak_sod, ends_with("rs"), starts_with("rain"), ends_with("wet")),
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Oak Plots RS Temperature and Rainfall Variables")
# For the most part rainfall and temperature variables are fairly weakly correlated. 
# Rainy days interpolations are correlated with hours below 10, hours 14-20, and average maximum temperature

pairs(select(oak_sod, ends_with("ds"), starts_with("rain")),
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Oak Plots DS Temperature and Rainfall Variables")
# Dry season temperature variables are only very weakly correlated with rainfall variables.

pairs(na.omit(select(oak_sod, inf_oak_ct, tot_bay, tot_lfct, avg_lfct, H.2005, us.H.2005, hrsblw14_wet, hrsblw10_wet, hrs1020_wet, avgtmax_wet, avgtmin_wet, hrs_blw10rs, avg_tmax_rs, avg_tmin_rs, hrs_abv25ds, avg_tmax_ds, avg_tmin_ds, rain_tot_v, rain_tot_2d, twi15m)),
      lower.panel = panel.cor, diag.panel = panel.hist)

#' # Oak Plots Path Models
#+ repeated measures plot-level path model ####
pairs(na.omit(select(oak_sod, inf_oak_ct, tot_bay, tot_lfct, avg_lfct, H.2005, H.2014, us.H.2005, us.H.2011, hrs_14_20_rs, avg_tmax_rs, avg_tmin_rs, avg_tmax_ds, avg_tmin_ds, rain_tot_v, rain_tot_2d, twi15m)),
      lower.panel = panel.cor, diag.panel = panel.hist)

# RS temperature hours 14-20, voronoi rain total, OS diversity 2005
## create binomial variable for oak infection
# oak_sod$oak.inf <- cbind(oak_sod$inf_oak_ct, oak_sod$uninf_oak_ct)

## subset data frame and transform/rescale variables
oakplots.sub <- select(oak_sod, plotid:tot_lfct, candens15m, twi15m, elev_m, veght_m, hrs_14_20_rs:us_rich11, rain_tot_2d, rain_tot_v, rain_days_v, ends_with("wet"))
summary(oakplots.sub)
pairs(select(oakplots.sub, -plotid), lower.panel = panel.cor, diag.panel = panel.hist)

# binomial variable for infected/uninfected oaks in each plot
oakplots.sub$oak.inf <- cbind(oakplots.sub$inf_oak_ct, oakplots.sub$uninf_oak_ct)
summary(filter(oak_sod, is.na(tot_lfct)))# lfct NAs are 0 bay plots
oakplots.sub$tot_lfct[is.na(oakplots.sub$tot_lfct)] <- 0
oakplots.sub$avg_lfct[is.na(oakplots.sub$avg_lfct)] <- 0

oakplots.sub$rain_tot_v.cm <- oakplots.sub$rain_tot_v/10
oakplots.sub$rain_tot_2d.cm <- oakplots.sub$rain_tot_2d/10
oakplots.sub$rain_tot_3d.cm <- oakplots.sub$rain_tot_3d/10

oakplots.sub$dys_14_20_rs <- oakplots.sub$hrs_14_20_rs/24
oakplots.sub$dys_blw10_rs <- oakplots.sub$hrs_blw10rs/24
oakplots.sub$dys_abv25_rs <- oakplots.sub$hrs_abv25rs/24
# oakplots.sub$tot_lfct.10 <- oakplots.sub$tot_lfct/10
# oakplots.sub$tot_lfct.100 <- oakplots.sub$tot_lfct/100
oakplots.sub$tot_lfct.log <- log1p(oakplots.sub$tot_lfct)
oakplots.sub$tot_bay.log <- log1p(oakplots.sub$tot_bay)
oakplots.sub$avg_lfct.log <- log1p(oakplots.sub$avg_lfct)
# oakplots.sub$sample_year <- as.factor(oakplots.sub$sample_year)

oakplots.modlist <- list(
      lme(tot_bay.log ~ twi15m + H.2005, random = 1|sample_year/plotid, 
            data = oakplots.sub, na.action = na.omit),
      
      lme(H.2005 ~ twi15m, random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      lme(dys_14_20_rs ~ tot_bay.log + twi15m, random = ~1|sample_year/plotid,
           data = oakplots.sub, na.action = na.omit),
      
      lme(tot_lfct.log ~ tot_bay.log + dys_14_20_rs + rain_tot_v.cm,
          random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      glmer(oak.inf ~ tot_lfct.log + dys_14_20_rs + rain_tot_v.cm + H.2005 + (1|sample_year) + (1|plotid),
            data = oakplots.sub, family = binomial(link = "logit"), 
            na.action = na.omit)
)

pairs(na.omit(select(oakplots.sub, inf_oak_ct, tot_bay, tot_bay.log, tot_lfct.log, avg_lfct, H.2005, H.2014, us.H.2005, us.H.2011, dys_14_20_rs, avg_tmax_rs, avgtmin_rs, avg_tmin_ds, avg_tmax_ds, rain_tot_v, twi15m)),
      lower.panel = panel.cor, diag.panel = panel.hist)

sem.fit(oakplots.modlist, data = oakplots.sub)
# Unacceptable fit at a p = 0.05 threshold (p = 0)
# Missing paths indicated: lfct ~ diversity and temp ~ diversity

sem.model.fits(oakplots.modlist)
# Most of the variation is explained by the random effects

(coef.table <- sem.coefs(oakplots.modlist, oakplots.sub))
(coef.table.std <- sem.coefs(oakplots.modlist, oakplots.sub, 
                             standardize = 'scale'))


#+ Add indicated missing paths ####
oakplots.modlist2 <- list(
#       glmer(tot_bay ~ twi15m + H.2005 + (1|plotid) + (1|sample_year), 
#             data = oakplots.sub, family = "poisson", na.action = na.omit),
#       
      lme(tot_bay.log ~ twi15m + H.2005, random = ~1|sample_year/plotid,
          data = oakplots.sub, na.action = na.omit),
      
      lme(H.2005 ~ twi15m, random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      lme(dys_14_20_rs ~ tot_bay.log + twi15m, random = ~1|sample_year/plotid,
          data = oakplots.sub, na.action = na.omit),
      
      lme(tot_lfct.log ~ tot_bay.log + dys_14_20_rs + rain_tot_v.cm + H.2005 + twi15m,
          random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      glmer(oak.inf ~ tot_lfct.log + dys_14_20_rs + rain_tot_v.cm + H.2005 + (1|plotid) + (1|sample_year),
            data = oakplots.sub, family = binomial(link = "logit"), 
            na.action = na.omit)
)

sem.lavaan(oakplots.modlist2, oakplots.sub)
# Unacceptable fit using Gaussian assumptions, p = 0

sem.fit(oakplots.modlist2, oakplots.sub)
# Unacceptable fit at p = 0.05 threshold (p = 0) using glmer bay count
## AIC: 151.59
## Missing paths: leaf count ~ diversity (p=0), temperature ~ diversity (p=.007)
# Unacceptable fit at p = 0.05 threshold (p = 0) using lme log bay count
## AIC: 123.97
## Missing paths: leaf count ~ diversity (p=0), leaf count ~ twi (p=.03)
### adding diversity: AIC = 91.59, p = .016, still missing twi (p=.0006)
### add twi: AIC = 78.72, p = 0.39, no missing paths with p < 0.05

sem.model.fits(oakplots.modlist2)
#' The random effects account for the majority of the variability in all models, except for the fixed effects do a decent job for predicting the total leaf count (marginal R^2 = 0.64).

#+ coefficients from revised model ####
(coef.table <- sem.coefs(oakplots.modlist2, oakplots.sub))

(coef.table.std <- sem.coefs(oakplots.modlist2, oakplots.sub, 
                             standardize = 'scale'))

pairs(na.omit(select(oakplots.sub, inf_oak_ct, tot_bay, tot_bay.log, tot_lfct.log, avg_lfct.log, H.2005, H.2014, us.H.2005, us.H.2011, dys_14_20_rs, avg_tmax_rs, avg_tmin_rs, avg_tmin_ds, avg_tmax_ds, rain_tot_v, twi15m)),
      lower.panel = panel.cor, diag.panel = panel.hist)
#'
#'
#' Change temperature variable to average for the rainy season, the leaf count variable to average for the plot, and the
#+ path model using rainy season average #### 
oakplots.modlist <- list(
#       glmer(tot_bay ~ twi15m + H.2005 + (1|sample_year) + (1|plotid), 
#             data = oakplots.sub, family = "poisson", na.action = na.omit),
      
      lme(tot_bay.log ~ twi15m + H.2005, random = ~1|sample_year/plotid,
          data = oakplots.sub, na.action = na.omit),
      
      lme(H.2005 ~ twi15m, random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      lme(avg_tmax_rs ~ tot_bay.log + twi15m, random = ~1|sample_year/plotid,
          data = oakplots.sub, na.action = na.omit),
      
      lme(avg_lfct.log ~ tot_bay.log + avg_tmax_rs + rain_tot_2d.cm + twi15m + H.2005,
          random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      glmer(oak.inf ~ avg_lfct.log + avg_tmax_rs + rain_tot_2d.cm + H.2005 + (1|sample_year) + (1|plotid),
            data = oakplots.sub, family = binomial(link = "logit"), 
            na.action = na.omit)
)

sem.missing.paths(oakplots.modlist, oakplots.sub)
sem.fit(oakplots.modlist, oakplots.sub)
sem.model.fits(oakplots.modlist)
sem.coefs(oakplots.modlist, oakplots.sub)
sem.coefs(oakplots.modlist, oakplots.sub, standardize = 'scale')


# Path model using dry season temperature variable ####
## model set with sample year as random slope

#' I went through these with Beth on 2/18/2016 and we looked at the difficult to interpret residuals together. The patterning is likely due to the large number of 0 values for symptomatic leaf count, number of bay laurel, and number of infected oak stems. Specifically, the residuals are off from the fitted values because of the many zeroes in the data, seen especially for the values near zero. The models seem to be properly parameterized with the available data, so Beth did not express too much concern with the residuals in the end. Sometimes funny looking residuals are just the way the data are, in this case, a lot of zero values.
#' 
#' Beth made a couple of suggestions
#' 
#'  1. Standardize the variables and rerun the models. Double check that the `scale` function does what it's supposed to do: rescale each variable to a mean of zero and variance of one. This will put the estimates into standard deviations from zero, but should not change the p-values.
#'  2. Look at the different components of each multivariate regression model as bivariate regressions. This can help understand the contributions of each variable to the model and their relationship with the response.
#'  3. Look up the calculation of Pearson's, and specifically how they are interpreted for a binomial model.

# Pearson Residuals in Binomial Models ####
#' Pearson's residuals are reported for the binomial model, and Beth wanted me to check how these are calculated to make sure the interpretation was correct. A google search lead me here, [](http://data.princeton.edu/wws509/notes/c3s8.html), which says Pearson residuals are, "*the difference between observed and fitted values divided by an estimate of the standard deviation of the observed value*", and, "They take into account the fact that different observations have different variances, but they make no allowance for additional variation arising from estimation of the parameters, in the way studentized residuals in classical linear models do." So, they are not fully standardized residuals, but this is the common standardization for binomial models. 
#' 
#' For base GLM there are five types of residuals available, c("deviance", "pearson", "working","response", "partial")

m1 <- lme(tot_bay.log ~ twi15m + H.2005, random = ~1|sample_year,
          data = oakplots.sub, na.action = na.omit)
summary(m1)
plot(m1, main = "Total # Bay Laurel ~ TWI + Diveristy Residuals")

m2 <- lme(H.2005 ~ twi15m, random = ~1|sample_year, data = oakplots.sub,
          na.action = na.omit)
summary(m2)
plot(m2, main = "Diversity ~ TWI Residuals")

m3 <- lme(avg_tmax_ds ~ tot_bay.log + twi15m, random = ~1|sample_year,
          data = oakplots.sub, na.action = na.omit)
summary(m3)
plot(m3)
# m3a <- lmer(avg_tmax_ds ~ tot_bay.log + twi15m + (1|sample_year),
#             data = oakplots.sub, na.action = na.omit)
# summary(m3a); AIC(m3a)
# plot(m3a)

m4 <- lme(tot_lfct.log ~ tot_bay.log + avg_tmax_ds + rain_tot_2d.cm + twi15m + H.2005, random = ~1|sample_year, data = oakplots.sub, 
          na.action = na.omit)
plot(m4)
summary(m4)
E2 <- resid(m4, type = "normalized")# equivalent of "standardized" in plot
# F2 <- fitted(m4)
# plot(F2, E2)
plot(E2 ~ tot_bay.log, data = oakplots.sub, ylab = "resids", xlab = "total # bay laurel")# this seems the most likely source of the pattern
plot(E2 ~ avg_tmax_ds, data = oakplots.sub, ylab = "resids", xlab = "temp")
plot(E2 ~ rain_tot_2d.cm, data = oakplots.sub, ylab = "resids", xlab = "rain")
plot(E2 ~ twi15m, data = oakplots.sub, ylab = "resids", xlab = "twi")
plot(E2 ~ H.2005, data = oakplots.sub, ylab = "resids", xlab = "Shannon's Index")
boxplot(E2 ~ sample_year, data = oakplots.sub, ylab = "resids", xlab = "year")
## These plots indicate that the trend in the residuals is not a result of the predictors, at least none of them are contributing the substantial amount abserved.

m4a <- lme(tot_lfct.log ~ tot_bay.log + avg_tmax_ds + rain_tot_2d.cm + twi15m + H.2005, random = ~1|sample_year/plotid, data = oakplots.sub, 
          na.action = na.omit)
# I believe this is allowing the intercept to vary among groups within sample year
plot(m4a)
summary(m4a)

m4b <- lmer(tot_lfct.log ~ tot_bay.log + avg_tmax_ds + rain_tot_2d.cm + twi15m + H.2005 + (1|sample_year/plotid), data = oakplots.sub)
# this formulation is allowing the intercept to vary among each group
summary(m4b)
plot(m4b)
image(getME(m4b, name = "L"))
image((getME(m4b, name = "A")))

m5 <- glmer(oak.inf ~ tot_lfct.log + avg_tmax_ds + rain_tot_2d.cm + H.2005 + (1|sample_year),
            data = oakplots.sub, family = binomial(link = "logit"), 
            na.action = na.omit)
summary(m5)
plot(m5)
summary(oakplots.sub)
# lfct: 0-8.8, tmax: 19.4-31.4, rain: 199.5-1529.1, H.2005: 0-1.6
# based on these ranges I'm supposing that the rain variable is throwing things off re: large eigenvalue

oakplots.sub$rain_tot_2d.dm <- oakplots.sub$rain_tot_2d.cm/10
m5a <- glmer(oak.inf ~ tot_lfct.log + avg_tmax_ds + rain_tot_2d.dm + H.2005 + (1|sample_year),
             data = oakplots.sub, family = binomial(link = "logit"), 
             na.action = na.omit)
# still fails to converge at set tolerance level, but avoids large eigenvalue
summary(m5a)# essentially the same fixed effects values taking into account rainfall transformation
plot(m5a)# same as for model m5

# add an observation-level (plot) random effect to account for additional variation...
m5b <- glmer(oak.inf ~ tot_lfct.log + avg_tmax_ds + rain_tot_2d.cm + H.2005 + (1|sample_year) + (1|plotid),
             data = oakplots.sub, family = binomial(link = "logit"), 
             na.action = na.omit)
summary(m5b)# notable changes to estimates and p-values; lots of variance attributed to plot-level random effect
plot(m5b)# this is...artsy? still with the nonlinear negative trend

image(getME(m5a, name = "L"))
image(getME(m5b, name = "L"))


# Fit Path Model
oakplots.modlist <- list(m1,m2,m3,m4,m5)

sem.fit(oakplots.modlist, oakplots.sub)
sem.model.fits(oakplots.modlist)
sem.coefs(oakplots.modlist, oakplots.sub)
sem.coefs(oakplots.modlist, oakplots.sub, standardize = 'scale')


# Oak Stem-Level Model ####
names(stems)
stems$plotid <- tolower(stems$plotid)
stems$species <- tolower(stems$species)


oak_stems <- stems %>% 
      select(plotid:status, dbh, dbh1, canker, sod_dead, location, year) %>%
      filter(status == "Alive" | status == "Dead", species != "lide", species != "unknown sp.", species != "umca")
oak_stems <- droplevels(oak_stems)
unique(oak_stems$species)

head(oak_stems,20)
summary(oak_stems)
tail(filter(oak_stems, is.na(dbh1)),20)#
filter(oak_stems, location == "Out")
oak_stems$infected <- ifelse(oak_stems$canker == 0, 0, 1)
names(oakplots.sub)

# join plot-level metrics to stem-level data
oak_stems <- left_join(oak_stems, 
                       select(oakplots.sub, plotid, sample_year, infected_bay_ct:tot_lfct, tot_bay.log, tot_lfct.log, twi15m, veght_m, avg_tmax_ds, hrs_abv20ds, avgtmin_wet, hrs1020_wet, H.2005, J.2005, os_rich05, rain_tot_v.cm),
                       by = c("plotid", "year"="sample_year"))
summary(oak_stems)






# Host Stem-Level Data ####
summary(stems)
names(stems)
stems.sub <- select(stems, plotid:status, dbh, dbh1, slc, canker, lide_lf_symp:location, date, year)
head(stems.sub,20)
summary(stems.sub)
stems.sub <- stems.sub %>% filter(status == "Alive" | status == "Dead", species != "LIDE", species != "unknown sp.")
stems.sub$plotid <- tolower(stems.sub$plotid)
stems.sub$species <- tolower(stems.sub$species)
stems.sub <- droplevels(stems.sub)
unique(stems.sub$species)
# Assign 'dbh1' values to missing 'dbh' values...?

stems.sub <- filter(stems.sub, plotid != "duer01", plotid != "ponti01", plotid != "sweet01")
stems.sub$plotid[stems.sub$plotid=="sugar02"] <- "sugar30"

names(plots_env)
stems.sub <- left_join(stems.sub, select(plots_env, plotid, candens15m, twi15m, elev_m, veght_m), by = "plotid")
summary(stems.sub)
filter(stems.sub, is.na(elev_m))
# stems.sub <- filter(stems.sub, plotid != "duer01", plotid != "ponti01", plotid != "sweet01")
# stems.sub$plotid[stems.sub$plotid=="sugar02"] <- "sugar30"

names(os.us.diversity)
stems.sub <- left_join(stems.sub, select(os.us.diversity, plotid, H.2005, J.2005, H.2014, J.2014, os_rich05, os_rich14), by = "plotid")
summary(stems.sub)
unique(filter(stems.sub, is.na(J.2014))$plotid)

names(summer_temps)
stems.sub <- left_join(stems.sub, select(summer_temps, plotid, sample_year, avg_tmax_ds, hrs_abv20ds, hrs_abv25ds), by = c("plotid", "year"="sample_year"))
summary(stems.sub)

names(wet_hrs)
stems.sub <- left_join(stems.sub, select(wet_hrs, plotid, sample_year, avgtmax_wet, avgtmin_wet, hrs1020_wet), by = c("plotid", "year"="sample_year"))
summary(stems.sub)

