#' # Path Modeling of Climatic Influence on Spillover
#' I'm using the `piecewiseSEM package, which operationalized the methods from Shipley 2009 & 2013 in particular for fitting and assessing path models that enable applying hierarchical modeling methods.
#' 
#+ load data and packages ####
load("pathmodel_data.RData")
library(plyr); library(dplyr); library(piecewiseSEM)
library(nlme); library(lme4); library(lmerTest)

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

#+ pairs plots of oak plots data ####
pairs(na.omit(select(oak_sod, starts_with("rain"))), 
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Oak Plots Rainfall Variables")
# The regression variables are the least correlated within totals & days.
# The 2D/3D interpolations are strongly correlated with each other, and strongly correlated with with the Voronoi polygon extractions.
# The days & total precipitation variables are fairly strongly correlated.

pairs(select(oak_sod, ends_with("rs"), ends_with("ds")),
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Oak Plots Temperature Variables")
# Overall, correlations among variables are relatively weak with a few strong pairs. For specific pairs of variables it may be okay to use one from the dry season and one from the rainy season.

pairs(select(oak_sod, ends_with("rs"), starts_with("rain")),
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Oak Plots RS Temperature and Rainfall Variables")
# For the most part rainfall and temperature variables are fairly weakly correlated. 
# Rainy days interpolations are correlated with hours below 10, hours 14-20, and average maximum temperature

pairs(select(oak_sod, ends_with("ds"), starts_with("rain")),
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Oak Plots DS Temperature and Rainfall Variables")
# Dry season temperature variables are only very weakly correlated with rainfall variables.

#' # Oak Plots Path Models
#+ repeated measures plot-level path model ####
pairs(na.omit(select(oak_sod, inf_oak_ct, tot_bay, tot_lfct, avg_lfct, H.2005, H.2014, us.H.2005, us.H.2011, hrs_14_20_rs, avg_tmax_rs, avg_tmin_rs, avg_tmax_ds, avg_tmin_ds, rain_tot_v, rain_tot_2d, twi15m)),
      lower.panel = panel.cor, diag.panel = panel.hist)

# RS temperature hours 14-20, voronoi rain total, OS diversity 2005
## create binomial variable for oak infection
oak_sod$oak.inf <- cbind(oak_sod$inf_oak_ct, oak_sod$uninf_oak_ct)
## subset data frame and transform/rescale variables
oakplots.sub <- select(oak_sod, plotid:tot_lfct, candens15m, twi15m, elev_m, veght_m, hrs_14_20_rs:rain_tot_v, rain_tot_psm:rain_days_v)
summary(oakplots.sub)
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

oakplots.modlist <- list(
      glmer(tot_bay ~ twi15m + H.2005 + (1|sample_year) + (1|plotid), 
            data = oakplots.sub, family = "poisson", na.action = na.omit),
      
      lme(H.2005 ~ twi15m, random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      lme(dys_14_20_rs ~ tot_bay + twi15m, random = ~1|sample_year/plotid,
           data = oakplots.sub, na.action = na.omit),
      
      lme(tot_lfct.log ~ tot_bay + dys_14_20_rs + rain_tot_v.cm,
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
oakplots.modlist <- list(
      #       glmer(tot_bay ~ twi15m + H.2005 + (1|sample_year) + (1|plotid), 
      #             data = oakplots.sub, family = "poisson", na.action = na.omit),
      
      lme(tot_bay.log ~ twi15m + H.2005, random = ~1|sample_year/plotid,
          data = oakplots.sub, na.action = na.omit),
      
      lme(H.2005 ~ twi15m, random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      lme(avg_tmax_ds ~ tot_bay.log + twi15m, random = ~1|sample_year/plotid,
          data = oakplots.sub, na.action = na.omit),
      
      lme(avg_lfct.log ~ tot_bay.log + avg_tmax_ds + rain_tot_2d.cm + twi15m + H.2005,
          random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      glmer(oak.inf ~ avg_lfct.log + avg_tmax_ds + rain_tot_2d.cm + H.2005 + (1|sample_year) + (1|plotid),
            data = oakplots.sub, family = binomial(link = "logit"), 
            na.action = na.omit)
)

sem.fit(oakplots.modlist, oakplots.sub)
sem.model.fits(oakplots.modlist)
sem.coefs(oakplots.modlist, oakplots.sub)
sem.coefs(oakplots.modlist, oakplots.sub, standardize = 'scale')
