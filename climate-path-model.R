#' # Path Modeling of Climatic Influence on Spillover
#' I'm using the `piecewiseSEM package, which operationalized the methods from Shipley 2009 & 2013 in particular for fitting and assessing path models that enable applying hierarchical modeling methods.
#' 
#+ load data and packages ####
load("pathmodel_data.RData")
library(plyr); library(dplyr); library(lme4); library(lmerTest)
library(nlme); library(lavaan); library(piecewiseSEM); library(pbkrtest)

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
      # p <- cor.test(x, y)$p.value
      # txt2 <- format(c(p, 0.123456789), digits = digits)[1]
      # txt2 <- paste("p= ", txt2, sep = "")
      # if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
      # text(0.5, 0.4, txt2)
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
# oak_sod$plotid <- as.factor(oak_sod$plotid)
str(oak_sod)
anyDuplicated(oak_sod)

#+ join wet-hours and lagged rainfall data ####
summary(wet_hrs)# since these are counts the NAs are in fact zeroes
wet_hrs[is.na(wet_hrs)] <- 0# introduce a lot of zeroes for some columns
anyDuplicated(wet_hrs)
oak_sod <- left_join(
      oak_sod, wet_hrs, by = c("plotid","sample_year"))
anyDuplicated(oak_sod)
# oak_sod$plotid <- as.factor(oak_sod$plotid)

summary(rs2d_lag1)
oak_sod <- left_join(
      oak_sod, select(rs2d_lag1, -rs2d_total), by = c("plotid","sample_year"))

# oak_sod$plotid <- as.factor(oak_sod$plotid)

summary(rs_temps_lag)
oak_sod <- left_join(
      oak_sod, 
      select(rs_temps_lag, plotid, sample_year, contains("tminus"), contains("t1")), by = c("plotid","sample_year"))

summary(summer_temps_lag)
oak_sod <- left_join(
      oak_sod,
      select(summer_temps_lag, plotid, sample_year, contains("tminus"), contains("t1")), by = c("plotid","sample_year"))

summary(wet_hrs_lag)
oak_sod <- left_join(
      oak_sod,
      select(wet_hrs_lag, plotid, sample_year, contains("tminus"), contains("t1")), by = c("plotid","sample_year"))

# oak_sod$plotid <- as.factor(oak_sod$plotid)

## drop 3d, PRISM, and regression rainfall variables based on discussion with Ryan; use a simple 2d spline interpolation or the Voronoi polygon values
# oak_sod <- select(oak_sod, -rain_tot_3d, -rain_tot_r, -rain_tot_psm)
# oak_sod <- filter(oak_sod, plotid != "ponti01", plotid != "sweet01", plotid != "bush01")

#' While thinking out loud with Monica yesterday I was puzzling over how bay laurel abundance relates to disease prevalence. In this case, where it's the number of stems in given area (the plot) abundance is essentially stem density, i.e. number of stems per plot area (225 m^2). Total basal area of UMCA and number of UMCA stems in the plot should be fairly strongly correlated, but I think this allows for a different interpretation of the relationships. Basal area is a clearer indicator of species dominance in a plot, since it is directly measuring the area occupied by a species, which therefore can't be occupied by another species. This means that plots with high bay laurel basal area are likely to have fewer other tree species, because bay laurel is occupying the available space. I am trying to make the distinction between this concept, and how bay stem density/abundance should be intepreted. It is possible to have a high density of small stems and for this high density of small stems to occupy a small portion of the plot area. We already know that total symptomatic leaf count and the number of bay laurel stems are strongly correlated, especially as stage of invasion progresses. I expect the same for basal area and symptomatic leaf count, though maybe a little less strong. This is because a couple of very large stems could produce a large basal area, but relatively small leaf count since we max out at around 200 for 60 second leaf counts. Perhaps this isn't "fair" representation of the epidemiological processes? 
#' 
#' A single large tree supporting a lot of infection can have a substantial epidemiological role in the plot, and likely beyond the plot. I think total leaf count underrepresents the role of these stems, whereas the average leaf count may underrepresent the role of many small stems in producing a high total inoculum load. I guess the leads to the question of whether I want to use the average inoculum load on bay laurel in the plot, capturing the role of a few large trees with potentially high inoculum load, or the total inoculum load on bay laurel, which ignores the role of a few large trees while capturing some of the effect of stem density. 
#' 
#+ calculate & add "in-plot" living umca basal area ####
umca_ba <- stems %>% filter(year > 2003, species == "umca", stem_status == "Alive", location == "In") %>% select(plotid, year, dbh_entry) %>% group_by(plotid, year) %>% summarise(umca_basal_area = sum(pi*dbh_entry))
oak_sod <- left_join(oak_sod, ungroup(umca_ba), by = c("plotid","sample_year" = "year"))

unique(filter(oak_sod, is.na(umca_basal_area))$plotid)
oak_sod$umca_basal_area[is.na(oak_sod$umca_basal_area)] <- 0
oak_sod$umca_basal_area[oak_sod$plotid == "yahng02" & oak_sod$sample_year == 2010] <- NA
# oak_sod$plotid <- as.factor(oak_sod$plotid)
# oak_sod <- mutate(oak_sod, avg_lfct = tot_lfct/tot_bay)
unique(filter(oak_sod, is.na(avg_lfct))$plotid)
oak_sod$avg_lfct[is.na(oak_sod$avg_lfct)] <- 0
oak_sod$avg_lfct[oak_sod$plotid == "yahng02" & oak_sod$sample_year == 2010] <- NA

#+ calculate and add importance value ####
# impt_vals <- rbind(stems %>% filter(location == "In", dbh >= 2.0, stem_status == "Alive") %>% select(plotid, species, dbh, stem_status, year), 
#                    untagged_stems %>% filter(stem_status == "Alive") %>% select(-date)) %>% 
#       group_by(plotid, species, year) %>% 
#       summarise(sp_abund = length(species), 
#                 sp_basal_area = sum(sum(pi*dbh^2/40000, na.rm = T))) %>% 
#       group_by(plotid, year) %>% 
#       mutate(total_abund = sum(sp_abund), 
#              total_basal_area = sum(sp_basal_area), 
#              imp_val =  
#                    0.5 * ( (sp_abund / total_abund) + 
#                                  (sp_basal_area / total_basal_area) ) )
# unique(impt_vals$year)
# impt_vals$year[impt_vals$year < 2005] <- 2005
# impt_vals <- select(ungroup(impt_vals), plotid, species, year, imp_val) %>% 
#       filter(species == "umca") %>% 
#       rename(umca_impval = imp_val) %>% 
#       select(-species)
# ## bad idea to try and use importance value because it isn't consistently calculated; we don't have all species for each year
# 
# oak_sod <- left_join(oak_sod, impt_vals, by = c("plotid", "sample_year"="year"))

# pairs(na.omit(select(oak_sod, inf_oak_ct, tot_bay, tot_lfct, avg_lfct, umca_basal_area, umca_impval, candens15m)),
#       lower.panel = panel.cor, diag.panel = panel.hist)

#+ pairs plots of oak plots data ####
# pairs(na.omit(select(oak_sod, inf_oak_ct, tot_bay, tot_lfct, avg_lfct, hrsblw14_wet, hrsblw10_wet, hrs1020_wet, hrsabv20_wet, avgtmax_wet, avgtmin_wet, rs2d_total, ss2d_total, rs_days, ss_days)), 
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Current Year Rainfall Variables")
# 
# pairs(na.omit(select(oak_sod, inf_oak_ct, tot_bay, tot_lfct, avg_lfct, starts_with("avg_rs2d"), contains("rs2d_t"))), 
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Lagged Rainfall Variables")
# 
# pairs(na.omit(select(oak_sod, inf_oak_ct, tot_bay, tot_lfct, avg_lfct, hrsblw14_wet, avgtmax_wet, hrs1020_wet, rs2d_total, ss2d_total, rs_days, ss_days, rs2d_total_t_minus1, avg_rs2d_t_minus1)), 
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Rainfall Variables")
# The regression variables are the least correlated within totals & days.
# The 2D/3D interpolations are strongly correlated with each other, and strongly correlated with with the Voronoi polygon extractions.
# The days & total precipitation variables are fairly strongly correlated.
# Total wet hours & hours below 10 & 14 are strongly correlated with rainy days.
# Wet hours below 14 & total hours are correlated with total rainfall r = 0.66 - 0.76
# Wet hours below 10 is correlated with total rainfall r = 0.61 - 0.68
# To me, this indicates the possibility of using the wet-hours count combination variable of hours below 14 C in place of rainfall. This is only weakly correlated with the average temperature variables calculated for wet days. The number of hours below 14 C on wet days definitely captures the number of rainy days, as it should (r = 0.94).

# pairs(na.omit(select(oak_sod, inf_oak_ct, tot_lfct, avg_lfct, umca_basal_area, ends_with("rs"), ends_with("ds"), hrsblw14_wet, avgtmin_wet, avgtmax_wet)),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Current Temperature Variables")
# Overall, correlations among variables are relatively weak with a few strong pairs. For specific pairs of variables it may be okay to use one from the dry season and one from the rainy season.

# pairs(na.omit(select(oak_sod, inf_oak_ct, inf_oak_2plus_ct, tot_lfct, avg_lfct, ends_with("ds"), ends_with("wet"))),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Wet-Days & Dry Season Temperature Variables")
# 
# pairs(na.omit(select(oak_sod, inf_oak_ct, tot_lfct, ends_with("rs"), starts_with("rain"), ends_with("wet"))),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots RS Temperature and Rainfall Variables")
# For the most part rainfall and temperature variables are fairly weakly correlated. 
# Rainy days interpolations are correlated with hours below 10, hours 14-20, and average maximum temperature

# pairs(na.omit(select(oak_sod, inf_oak_ct, tot_lfct, ends_with("ds"), contains("total"))),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots DS Temperature and Rainfall Variables")
# Dry season temperature variables are only very weakly correlated with rainfall variables.

## Correlations of lagged variables 
# pairs(na.omit(select(oak_sod, inf_oak_ct, inf_oak_2plus_ct, tot_lfct, contains("t1"))),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Average Lagged Variables")
# pairs(na.omit(select(oak_sod, inf_oak_ct, inf_oak_2plus_ct, tot_lfct, contains("tminus1"))),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Lagged Variables t-1")
# pairs(na.omit(select(oak_sod, inf_oak_ct, inf_oak_2plus_ct, tot_lfct, contains("tminus2"))),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Lagged Variables t-2")
# pairs(na.omit(select(oak_sod, inf_oak_ct, inf_oak_2plus_ct, tot_lfct, contains("rs2d"))),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Lagged Rainfall Variables")


# pairs(na.omit(select(oak_sod, inf_oak_ct, tot_bay, tot_lfct, umca_basal_area, H.2005, H.2014, hrsblw14_wet, hrsblw10_wet, hrs1020_wet, avgtmax_wet, avgtmin_wet, hrs_blw10rs, avg_tmax_rs, avg_tmin_rs, hrs_abv25ds, avg_tmax_ds, avg_tmin_ds, rs2d_total, ss2d_total, twi15m)),
#       lower.panel = panel.cor, diag.panel = panel.hist)

#' # Oak Plots Path Models
#' These are models of that terminate in disease prevalence of oak species at the plot-level. 
#' 
#' I am going to model each year independently, and use a "repeated measures" approach that accounts for correlation among plots sampled during the same season by applying the sample year as a random effect. This is effectively a single model implemented to account for the same plots being sampled repeatedly, i.e. once each year.
#' 
#+ repeated measures plot-level path model ####
# pairs(na.omit(select(oak_sod, inf_oak_ct, tot_bay, tot_lfct, H.2005, H.2014, us.H.2005, us.H.2011, hrsblw14_wet, hrs1020_wet, avg_tmax_ds, avg_tmin_ds, hrs_abv25ds, rs2d_total, rs_v_total, twi15m)),
#       lower.panel = panel.cor, diag.panel = panel.hist)

# RS temperature hours 14-20, voronoi rain total, OS diversity 2005
## create binomial variable for oak infection
# oak_sod$oak.inf <- cbind(oak_sod$inf_oak_ct, oak_sod$uninf_oak_ct)

#+ subset data frame and transform/rescale/standardize variables ####
oakplots.sub <- select(oak_sod, plotid:inf_oak_ct, inf_oak_2plus_ct, avg_lfct:week, umca_basal_area, candens15m:veght_m, contains("wet"), contains("rs"), contains("ds"), contains("tminus1"), H.2005, H.2014, J.2005, J.2014, os_rich05, os_rich14, contains("rs2d"), contains("ss2d"))
oakplots.sub <- select(oakplots.sub, -contains("t_minus2"), -contains("t_minus3"), -contains("t1t2"), -contains("t1t2t3"), -contains("tminus2"), -contains("tminus3"))
summary(oakplots.sub)

# filter(oakplots.sub, is.na(elev_m))
# filter(oakplots.sub, is.na(avg_tmax_rs))
# filter(oakplots.sub, is.na(us_rich05))# one plot, never has infected oak
# filter(oakplots.sub, is.na(os_rich05))# this shouldn't have NAs b/c there should always be at least one species...

# pairs(na.omit(select(oakplots.sub, -plotid:-uninf_oak_ct)), lower.panel = panel.cor, diag.panel = panel.hist)

# binomial variable for infected/uninfected oaks in each plot
oakplots.sub$oak.inf <- cbind(oakplots.sub$inf_oak_ct, oakplots.sub$uninf_oak_ct)
# summary(filter(oak_sod, is.na(tot_lfct)))# lfct NAs are 0 bay plots
# oakplots.sub$tot_lfct[is.na(oakplots.sub$tot_lfct)] <- 0
# oakplots.sub$avg_lfct[is.na(oakplots.sub$avg_lfct)] <- 0

oakplots.sub$rs2d.cm <- oakplots.sub$rs2d_total/10
oakplots.sub$ss2d.cm <- oakplots.sub$ss2d_total/10
oakplots.sub$avg_rs_tminus1.cm <- oakplots.sub$avg_rs2d_t_minus1/10
# oakplots.sub$avg_rs_tminus2.cm <- oakplots.sub$avg_rs2d_t_minus2/10
oakplots.sub$tot_rs_tminus1.cm <- oakplots.sub$rs2d_total_t_minus1/10
# oakplots.sub$tot_rs_tminus2.cm <- oakplots.sub$rs2d_total_t_minus2/10
# oakplots.sub$avg_rs_t1t2.cm <- oakplots.sub$avg_rs2d_t1t2/10
# oakplots.sub$rain_tot_2d.cm <- oakplots.sub$rain_tot_2d/10
# oakplots.sub$rain_tot_3d.cm <- oakplots.sub$rain_tot_3d/10

# oakplots.sub$dysblw14_wet <- oakplots.sub$hrsblw14_wet/24
oakplots.sub$dys_abv25ds <- oakplots.sub$hrs_abv25ds/24
# oakplots.sub$dys_blw10_rs <- oakplots.sub$hrs_blw10rs/24
# oakplots.sub$dys_abv25_rs <- oakplots.sub$hrs_abv25rs/24
# oakplots.sub$tot_lfct.10 <- oakplots.sub$tot_lfct/10
# oakplots.sub$tot_lfct.100 <- oakplots.sub$tot_lfct/100
oakplots.sub$tot_lfct.log <- log1p(oakplots.sub$tot_lfct)
oakplots.sub$tot_bay.log <- log1p(oakplots.sub$tot_bay)
oakplots.sub$avg_lfct.log <- log1p(oakplots.sub$avg_lfct)
oakplots.sub$umca_ba.log <- log1p(oakplots.sub$umca_basal_area)
oakplots.sub$twi15m.log <- log(oakplots.sub$twi15m)
oakplots.sub$avg_hrs1422_wet_t_t1.log <- log1p(oakplots.sub$avg_hrs1422_wet_t_t1)
oakplots.sub$hrs1422_wet_tminus1.log <- log1p(oakplots.sub$hrs1422_wet_tminus1)
# oakplots.sub$avg_lfct.log <- log1p(oakplots.sub$avg_lfct)
# oakplots.sub$sample_year <- as.factor(oakplots.sub$sample_year)
# oakplots.sub$oak.inf.2plus <- cbind(oakplots.sub$inf_oak_2plus_ct, (oak_sod$tot_oak_ct - oakplots.sub$inf_oak_2plus_ct))

summary(oakplots.sub)
anyDuplicated(oakplots.sub)
oakplots.sub[1565:1570,]


## rough path-model 1, how does this work? ####
oakplots.modlist1 <- list(
      lme(tot_bay.log ~ twi15m + H.2005, random = ~1|sample_year/plotid, 
            data = oakplots.sub, na.action = na.omit),
      
      lme(H.2005 ~ twi15m, random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      lme(dys_abv25ds ~ tot_bay.log + twi15m, random = ~1|sample_year/plotid,
           data = oakplots.sub, na.action = na.omit),
      
      lme(tot_lfct.log ~ tot_bay.log + dys_abv25ds + ss2d.cm + twi15m,
          random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      glmer(oak.inf ~ tot_lfct.log + dys_abv25ds + ss2d.cm + H.2005 + (1|sample_year) + (1|plotid),
            data = oakplots.sub, family = binomial(link = "logit"), 
            na.action = na.omit)
)

# pairs(na.omit(select(oakplots.sub, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, H.2005, H.2014, dys_abv25ds, avg_tmax_rs, avg_tmax_ds, ss2d.cm, twi15m, twi15m.log)),
#       lower.panel = panel.cor, diag.panel = panel.hist)

sem.fit(oakplots.modlist1, data = oakplots.sub)
# Unacceptable fit at a p = 0.05 threshold (p = 0)
# Missing paths indicated: lfct ~ diversity and temp ~ rainfall
sem.fit(oakplots.modlist1, data = oakplots.sub, 
        corr.errors = c("dys_abv25ds ~~ ss2d.cm",
                        "tot_lfct.log ~~ H.2005"))
sem.model.fits(oakplots.modlist1)
# Most of the variation is explained by the random effects

(coef.table <- sem.coefs(oakplots.modlist1, oakplots.sub,
                         corr.errors = c("dys_abv25ds ~~ ss2d.cm",
                                         "tot_lfct.log ~~ H.2005")))
(coef.table.std <- sem.coefs(oakplots.modlist1, oakplots.sub, 
                             standardize = 'scale'))
#'


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
#' 
#'
#' 
#' ## Single Season Models
#' In this section I'm going to paramterized single season models, i.e. model each sample season separately instead of in a repeated measures/random-effect framework. I'm going to start with 2005 because I deemed that to be the sesaon when plot-establishment was completed. Something along the lines of 600+ tagged stems were added to the database that year, many more than during subsequent visits.
#'


#+ make dataframe for each season ####
names(oakplots.sub)
oakplots.2005 <- oakplots.sub %>% filter(sample_year == 2005) %>% 
      select(plotid, oak.inf, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, umca_ba.log, candens15m, H.2005:os_rich14, contains("wet"), contains("ds"), contains("rs"), twi15m, twi15m.log)
summary(oakplots.2005)
anyDuplicated(oakplots.2005)
oakplots.2005[200:212,]
# pairs(select(oakplots.2005, -plotid, -oak.inf, -oak.inf.2plus),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Data 2005 Sample Season")

oakplots.2006 <- (oakplots.sub %>% filter(sample_year == 2006) %>% 
                        select(plotid, oak.inf, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, umca_ba.log, candens15m, H.2005:os_rich14, contains("wet"), contains("ds"), contains("rs"), twi15m, twi15m.log))
anyDuplicated(oakplots.2006)
# pairs(select(oakplots.2006, -plotid, -oak.inf, -oak.inf.2plus),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Data 2006 Sample Season")

oakplots.2007 <- (oakplots.sub %>% filter(sample_year == 2007) %>% 
                        select(plotid, oak.inf, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, umca_ba.log, candens15m, H.2005:os_rich14, contains("wet"), contains("ds"), contains("rs"), twi15m, twi15m.log))
anyDuplicated(oakplots.2007)
# pairs(select(oakplots.2007, -plotid, -oak.inf, -oak.inf.2plus),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Data 2007 Sample Season")

oakplots.2008 <- (oakplots.sub %>% filter(sample_year == 2008) %>% 
                        select(plotid, oak.inf, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, umca_ba.log, candens15m, H.2005:os_rich14, contains("wet"), contains("ds"), contains("rs"), twi15m, twi15m.log))
anyDuplicated(oakplots.2008)
# pairs(select(oakplots.2008, -plotid, -oak.inf, -oak.inf.2plus),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Data 2008 Sample Season")

oakplots.2009 <- (oakplots.sub %>% filter(sample_year == 2009) %>% 
                        select(plotid, oak.inf, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, umca_ba.log, candens15m, H.2005:os_rich14, contains("wet"), contains("ds"), contains("rs"), twi15m, twi15m.log))
anyDuplicated(oakplots.2009)
# pairs(select(oakplots.2009, -plotid, -oak.inf, -oak.inf.2plus),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Data 2009 Sample Season")

oakplots.2010 <- (oakplots.sub %>% filter(sample_year == 2010) %>% 
                        select(plotid, oak.inf, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, umca_ba.log, candens15m, H.2005:os_rich14, contains("wet"), contains("ds"), contains("rs"), twi15m, twi15m.log))
anyDuplicated(oakplots.2010)
# pairs(na.omit(select(oakplots.2010, -plotid, -oak.inf, -oak.inf.2plus)),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Data 2010 Sample Season")

oakplots.2011 <- (oakplots.sub %>% filter(sample_year == 2011) %>% 
                        select(plotid, oak.inf, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, umca_ba.log, candens15m, H.2005:os_rich14, contains("wet"), contains("ds"), contains("rs"), twi15m, twi15m.log))
anyDuplicated(oakplots.2011)
# pairs(select(oakplots.2011, -plotid, -oak.inf, -oak.inf.2plus),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Data 2011 Sample Season")

oakplots.2012 <- (oakplots.sub %>% filter(sample_year == 2012) %>% 
                        select(plotid, oak.inf, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, umca_ba.log, candens15m, H.2005:os_rich14, contains("wet"), contains("ds"), contains("rs"), twi15m, twi15m.log))
anyDuplicated(oakplots.2012)
# pairs(select(oakplots.2012, -plotid, -oak.inf, -oak.inf.2plus),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Data 2012 Sample Season")

oakplots.2014 <- (oakplots.sub %>% filter(sample_year == 2014) %>% 
                        select(plotid, oak.inf, inf_oak_ct, tot_bay.log, tot_lfct.log, avg_lfct.log, umca_ba.log, candens15m, H.2005:os_rich14, contains("wet"), contains("ds"), contains("rs"), twi15m, twi15m.log))
anyDuplicated(oakplots.2014)
# pairs(select(oakplots.2014, -plotid, -oak.inf, -oak.inf.2plus),
#       lower.panel = panel.cor, diag.panel = panel.hist,
#       main = "Oak Plots Data 2014 Sample Season")

save.image("climate_path_model_data.RData")
#'


#+ path model 2005 ####
oakplots2005.modlist1 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2005),
      m2 <- lm(tot_bay.log ~ twi15m.log, data = oakplots.2005),
      m3 <- lm(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log, data = oakplots.2005),
      m6 <- lm(candens15m ~ tot_bay.log, data = oakplots.2005),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2005),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2005),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2005)
)
par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
E3 <- resid(m3)
plot(E3 ~ tot_bay.log, data = oakplots.2005)
plot(E3 ~ hrs_abv25ds, data = oakplots.2005)
plot(E3 ~ twi15m.log, data = oakplots.2005)

plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")

# variogram analysis of leaf count residuals ####
oak_plots_shp <- plots_shp[tolower(plots_shp$PLOT_ID) %in% unique(oak_sod$plotid),]
oak_plots_shp$PLOT_ID <- tolower(oak_plots_shp$PLOT_ID)
oakplots.2005 <- left_join(oakplots.2005,
                           select(oak_plots_shp@data, PLOT_ID, X, Y),
                           by = c("plotid"="PLOT_ID"))
coords2005 <- as.matrix(oakplots.2005[,c("X","Y")])
leafcount <- oakplots.2005$tot_lfct.log
# plot(coords2005, pch = 1, cex = wethours/500, col = "darkgreen",
#      xlab = "Easting (m)", ylab = "Northing (m)")
summary(m3)
library(spBayes); library(classInt); library(RColorBrewer); library(geoR)

lfct.resid <- resid(m3)
max.dist <- 0.33 * max(iDist(coords2005))# iDist is from spBayes package
bins <- 100

vario.lfct <- variog(coords = coords2005, 
                     data = leafcount,
                     uvec = (seq(0, max.dist, length = bins)))
str(vario.lfct)
fit.lfct <- variofit(
      vario.lfct, 
      ini.cov.pars = c(300, 200/-log(0.5)), 
      cov.model = "exponential", 
      minimisation.function = "optim",
      weights = "equal"
)

vario.lfct.resid <- variog(coords = coords2005, 
                           data = lfct.resid,
                           uvec = (seq(0, max.dist, length = bins)))
str(vario.lfct.resid)
fit.lfct.resid.vario <- variofit(
      vario.lfct.resid, 
      #ini.cov.pars = c(300, 200/-log(0.5)), 
      cov.model = "exponential", 
      minimisation.function = "optim",
      weights = "equal"
                                 )
par(mfrow = c(1,2))
plot(vario.lfct)
lines(fit.lfct)
abline(h = fit.lfct$nugget, col = "blue")# nugget
abline(h = fit.lfct$cov.pars[1] + fit.lfct$nugget, col = "green")# sill
abline(v = -log(.05) * fit.lfct$cov.pars[2], col = "red3")# range

plot(vario.lfct.resid)
lines(fit.lfct.resid.vario)
abline(h = fit.lfct.resid.vario$nugget, col = "blue")# nugget
abline(h = fit.lfct.resid.vario$cov.pars[1] + fit.lfct.resid.vario$nugget, col = "green")# sill
abline(v = -log(.05) * fit.lfct.resid.vario$cov.pars[2], col = "red3")# range

library(gstat); library(spdep); library(maptools)
sp.data.2005 <- as.data.frame(oakplots.2005)
coordinates(sp.data.2005) <- c("X","Y")
vario.lfct2 <- variogram(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log,
                         data = sp.data.2005, cutoff = max.dist,
                         width = 50, alpha = (0:3) * 45)
fit.vario.lfct2 <- fit.variogram(vario.lfct2, 
                                 vgm(psill = 3, model = "Exp", 
                                     range = 11000, nugget = 1.5))
print(plot(vario.lfct2, fit.vario.lfct2))

sel <- plot(variogram(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log,
                      data = sp.data.2005, cutoff = 3000,
                      width = 50, cloud = T), digitize = T)
plot(sel, sp.data.2005)


sp.data.2005@data$lfct.resid <- lfct.resid
bubble(sp.data.2005, zcol = "lfct.resid")

# nearest neighbor distances ####
sp.data.2005.kn2 <- knn2nb(knearneigh(coordinates(sp.data.2005), k = 2),
                           row.names=sp.data.2005@data$plotid)
plot(sp.data.2005, main = "K-nearest-neighbors = 2")
plot(sp.data.2005.kn2, coords2005, add = T)
sp.data.2005.kn2.w <- nb2listw(sp.data.2005.kn2)

sp.data.2005.kn3 <- knn2nb(knearneigh(coordinates(sp.data.2005), k = 3),
                           row.names=sp.data.2005@data$plotid)
plot(sp.data.2005, main = "K-nearest-neighbors = 3")
plot(sp.data.2005.kn3, coords2005, add = T)
sp.data.2005.kn3.w <- nb2listw(sp.data.2005.kn3)

sp.data.2005.kn4 <- knn2nb(knearneigh(coordinates(sp.data.2005), k = 4),
                           row.names=sp.data.2005@data$plotid)
plot(sp.data.2005, main = "K-nearest-neighbors = 4")
plot(sp.data.2005.kn4, coords2005, add = T)
sp.data.2005.kn4.w <- nb2listw(sp.data.2005.kn4, style = "B")

sp.data.2005.kn5 <- knn2nb(knearneigh(coordinates(sp.data.2005), k = 5),
                           row.names=sp.data.2005@data$plotid)
plot(sp.data.2005, main = "K-nearest-neighbors = 5")
plot(sp.data.2005.kn5, coords2005, add = T)
sp.data.2005.kn5.w <- nb2listw(sp.data.2005.kn5, style = "B")


moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn2.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn2.w)

moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn3.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn3.w)

moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn4.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn4.w)
moran.plot(sp.data.2005$lfct.resid, listw = sp.data.2005.kn5.w)

moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn5.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn5.w)

moran.test(sp.data.2005$tot_lfct.log, listw = sp.data.2005.kn5.w)
moran.plot(sp.data.2005$tot_lfct.log, listw = sp.data.2005.kn5.w)
geary.test(sp.data.2005$tot_lfct.log, listw = sp.data.2005.kn5.w)

sp.data.2005.kn1 <- knn2nb(knearneigh(coordinates(sp.data.2005), k = 1),
                           row.names=sp.data.2005@data$plotid)
sp.data.2005.kn1.w <- nb2listw(sp.data.2005.kn1)
dist <- unlist(nbdists(sp.data.2005.kn1, coords2005))
summary(dist)
max.k1 <- max(dist)
moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn1.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kn1.w)

sp.data.2005.kd1 <- dnearneigh(coords2005, d1 = 0, d2 = 0.75*max.k1, 
                               row.names = sp.data.2005@data$plotid)
sp.data.2005.kd1.w <- nb2listw(sp.data.2005.kd1, zero.policy=T)
plot(sp.data.2005)
plot(sp.data.2005.kd1, coords2005, add = T)
# moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd1.w)

sp.data.2005.kd2 <- dnearneigh(coords2005, d1 = 0, d2 = max.k1, 
                               row.names = sp.data.2005@data$plotid)
sp.data.2005.kd2.w <- nb2listw(sp.data.2005.kd2, style = "B")
plot(sp.data.2005)
plot(sp.data.2005.kd2, coords2005, add = T)
moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd2.w)
moran.plot(sp.data.2005$lfct.resid, listw = sp.data.2005.kd2.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd2.w)

sp.data.2005.kd3 <- dnearneigh(coords2005, d1 = 0, d2 = 1.5*max.k1, 
                               row.names = sp.data.2005@data$plotid)
sp.data.2005.kd3.w <- nb2listw(sp.data.2005.kd3, style = "B")
plot(sp.data.2005)
plot(sp.data.2005.kd3, coords2005, add = T)
moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd3.w)
moran.plot(sp.data.2005$lfct.resid, listw = sp.data.2005.kd3.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd3.w)

lm.morantest(m3, sp.data.2005.kd3.w)
lm.LMtests(m3, sp.data.2005.kd3.w, test = "all")
## indicates error & lag dependence, but stronger error dependence

sp.data.2005.kd4 <- dnearneigh(coords2005, d1 = 0, d2 = 2*max.k1, 
                               row.names = sp.data.2005@data$plotid)
sp.data.2005.kd4.w <- nb2listw(sp.data.2005.kd4, style = "B")
moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd4.w)
moran.plot(sp.data.2005$lfct.resid, listw = sp.data.2005.kd4.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd4.w)

sp.data.2005.kd5 <- dnearneigh(coords2005, d1 = 0, d2 = 2.5*max.k1, 
                               row.names = sp.data.2005@data$plotid)
sp.data.2005.kd5.w <- nb2listw(sp.data.2005.kd5, style = "B")
moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd5.w)
moran.plot(sp.data.2005$lfct.resid, listw = sp.data.2005.kd5.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd5.w)
globalG.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd5.w)

sp.data.2005.kd6 <- dnearneigh(coords2005, d1 = 0, d2 = 3*max.k1, 
                               row.names = sp.data.2005@data$plotid)
sp.data.2005.kd6.w <- nb2listw(sp.data.2005.kd6, style = "B")
plot(sp.data.2005)
plot(sp.data.2005.kd6, coords2005, add = T)
moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd6.w)
moran.plot(sp.data.2005$lfct.resid, listw = sp.data.2005.kd6.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd6.w)

sp.data.2005.kd7 <- dnearneigh(coords2005, d1 = 0, d2 = 3.25*max.k1, 
                               row.names = sp.data.2005@data$plotid)
sp.data.2005.kd7.w <- nb2listw(sp.data.2005.kd7, style = "B")
plot(sp.data.2005)
plot(sp.data.2005.kd7, coords2005, add = T)
moran.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd7.w)
moran.plot(sp.data.2005$lfct.resid, listw = sp.data.2005.kd7.w)
geary.test(sp.data.2005$lfct.resid, listw = sp.data.2005.kd7.w)
localmoran(sp.data.2005$lfct.resid, listw = sp.data.2005.kd7.w)
localmoran.exact(m3, nb = sp.data.2005.kd7)


# umca basal area 2005 model ####
oakplots2005.modlist2 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2005),
      m2 <- lm(umca_ba.log ~ twi15m.log, data = oakplots.2005),
      m3 <- lm(tot_lfct.log ~ umca_ba.log + hrs_abv25ds + twi15m.log, data = oakplots.2005),
      m6 <- lm(candens15m ~ umca_ba.log, data = oakplots.2005),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2005),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2005),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2005)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Dominance ~")
E2 <- resid(m2)
plot(E2 ~ twi15m.log, data = oakplots.2005)
## distinct non-normality induced by the zero values for bay laurel basal area
plot(m3, main = "Leaf Symptoms ~")
E3 <- resid(m3)
plot(E3 ~ tot_bay.log, data = oakplots.2005)
plot(E3 ~ hrs_abv25ds, data = oakplots.2005)
plot(E3 ~ twi15m.log, data = oakplots.2005)
## patterns in residuals produced by zero values for bay laurel leaf counts

plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature")

sem.fit(oakplots2005.modlist1, data = oakplots.2005)
sem.fit(oakplots2005.modlist2, data = oakplots.2005)
## the model set using bay laurel basal area produces a better fit by about 3.9 points of AICc
# sem.fit(oakplots2005.modlist, data = oakplots.2005,
#         corr.errors = "dys_abv25ds ~~ avg_rs_tminus1.cm")
sem.model.fits(oakplots2005.modlist1)
sem.model.fits(oakplots2005.modlist2)
(coef.table.2005 <- sem.coefs(oakplots2005.modlist1, oakplots.2005))
(coef.table.2005.scale <- sem.coefs(oakplots2005.modlist1, oakplots.2005, standardize = 'scale'))
sem.coefs(oakplots2005.modlist2, oakplots.2005)
sem.coefs(oakplots2005.modlist2, oakplots.2005, standardize = 'scale')
# sem.coefs(oakplots2005.modlist, oakplots.2005, standardize = 'scale', 
#           corr.errors = "dys_abv25ds ~~ avg_rs_tminus1.cm")


#+ path model 2006 ####
oakplots2006.modlist1 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2006),
      m2 <- lm(tot_bay.log ~ twi15m.log, data = oakplots.2006),
      m3 <- lm(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log, data = oakplots.2006),
      m6 <- lm(candens15m ~ tot_bay.log, data = oakplots.2006),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2006),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2006),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2006)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")

oakplots2006.modlist2 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2006),
      m2 <- lm(umca_ba.log ~ twi15m.log, data = oakplots.2006),
      m3 <- lm(tot_lfct.log ~ umca_ba.log + hrs_abv25ds + twi15m.log, data = oakplots.2006),
      m6 <- lm(candens15m ~ umca_ba.log, data = oakplots.2006),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2006),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2006),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2006)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


sem.model.fits(oakplots2006.modlist1)
sem.model.fits(oakplots2006.modlist2)

sem.fit(oakplots2006.modlist1, data = oakplots.2006)
sem.fit(oakplots2006.modlist2, data = oakplots.2006)

(coef.table.2006 <- sem.coefs(oakplots2006.modlist1, data = oakplots.2006))
(coef.table.2006.scale <- sem.coefs(oakplots2006.modlist1, data = oakplots.2006, standardize = 'scale'))
sem.coefs(oakplots2006.modlist2, data = oakplots.2006)
sem.coefs(oakplots2006.modlist2, data = oakplots.2006, standardize = 'scale')

# sem.fit(oakplots2006.modlist, data = oakplots.2006,
#         corr.errors = "dys_abv25ds ~~ avg_rs_tminus1.cm")
# sem.coefs(oakplots2006.modlist, data = oakplots.2006,
#           corr.errors = "dys_abv25ds ~~ avg_rs_tminus1.cm")


#+ path model 2007 ####
oakplots2007.modlist1 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2007),
      m2 <- lm(tot_bay.log ~ twi15m.log, data = oakplots.2007),
      m3 <- lm(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log, data = oakplots.2007),
      m6 <- lm(candens15m ~ tot_bay.log, data = oakplots.2007),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2007),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2007),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2007)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
E3 <- resid(m3)
oakplots.2007[58,]# outlier with quite a few UMCA, but zero leaf count
plot(E3 ~ tot_bay.log, data = oakplots.2007)
# there may be something non-linear going on with bay laurel abundance
plot(E3 ~ hrs_abv25ds, data = oakplots.2007)
plot(E3 ~ twi15m.log, data = oakplots.2007)

plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")

oakplots2007.modlist2 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2007),
      m2 <- lm(umca_ba.log ~ twi15m.log, data = oakplots.2007),
      m3 <- lm(tot_lfct.log ~ umca_ba.log + hrs_abv25ds + twi15m.log, data = oakplots.2007),
      m6 <- lm(candens15m ~ umca_ba.log, data = oakplots.2007),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2007),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2007),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2007)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


sem.fit(oakplots2007.modlist1, oakplots.2007)
sem.fit(oakplots2007.modlist2, oakplots.2007)

(coef.table.2007 <- sem.coefs(oakplots2007.modlist1, oakplots.2007))
(coef.table.2007.scale <- sem.coefs(oakplots2007.modlist1, oakplots.2007, standardize = 'scale'))
sem.coefs(oakplots2007.modlist2, oakplots.2007)
sem.coefs(oakplots2007.modlist2, oakplots.2007, standardize = 'scale')


#+ path model 2008 ####
oakplots2008.modlist1 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2008),
      m2 <- lm(tot_bay.log ~ twi15m.log, data = oakplots.2008),
      m3 <- lm(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log, data = oakplots.2008),
      m6 <- lm(candens15m ~ tot_bay.log, data = oakplots.2008),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2008),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2008),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2008)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
E3 <- resid(m3)
plot(E3 ~ tot_bay.log, data = oakplots.2008)
## strong negative relationship between UMCA abundance & residuals
## perhaps due to zero values
plot(E3 ~ hrs_abv25ds, data = oakplots.2008)
plot(E3 ~ twi15m.log, data = oakplots.2008)

plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


oakplots2008.modlist2 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2008),
      m2 <- lm(umca_ba.log ~ twi15m.log, data = oakplots.2008),
      m3 <- lm(tot_lfct.log ~ umca_ba.log + hrs_abv25ds + twi15m.log, data = oakplots.2008),
      m6 <- lm(candens15m ~ umca_ba.log, data = oakplots.2008),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2008),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2008),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2008)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


sem.fit(oakplots2008.modlist1, oakplots.2008)
sem.fit(oakplots2008.modlist2, oakplots.2008)

(coef.table.2008 <- sem.coefs(oakplots2008.modlist1, oakplots.2008))
(coef.table.2008.scale <- sem.coefs(oakplots2008.modlist1, oakplots.2008, standardize = 'scale'))
sem.coefs(oakplots2008.modlist2, oakplots.2008)
sem.coefs(oakplots2008.modlist2, oakplots.2008, standardize = 'scale')


#+ path model 2009 ####
oakplots2009.modlist1 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2009),
      m2 <- lm(tot_bay.log ~ twi15m.log, data = oakplots.2009),
      m3 <- lm(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log, data = oakplots.2009),
      m6 <- lm(candens15m ~ tot_bay.log, data = oakplots.2009),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2009),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2009),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2009)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
E3 <- resid(m3)
plot(E3 ~ tot_bay.log, data = oakplots.2009)
## strong negative relationship between UMCA abundance & residuals
## perhaps due to the zero values
plot(update(m3, .~., data = filter(oakplots.2009, tot_lfct.log > 0)))
plot(E3 ~ hrs_abv25ds, data = oakplots.2009)
plot(E3 ~ twi15m.log, data = oakplots.2009)

plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


oakplots2009.modlist2 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2009),
      m2 <- lm(umca_ba.log ~ twi15m.log, data = oakplots.2009),
      m3 <- lm(tot_lfct.log ~ umca_ba.log + hrs_abv25ds + twi15m.log, data = oakplots.2009),
      m6 <- lm(candens15m ~ umca_ba.log, data = oakplots.2009),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2009),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2009),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2009)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


sem.fit(oakplots2009.modlist1, oakplots.2009)
sem.fit(oakplots2009.modlist2, oakplots.2009)

sem.coefs(oakplots2009.modlist1, oakplots.2009)
sem.coefs(oakplots2009.modlist1, oakplots.2009, standardize = 'scale')
sem.coefs(oakplots2009.modlist2, oakplots.2009)
sem.coefs(oakplots2009.modlist2, oakplots.2009, standardize = 'scale')


#+ path model 2010 ####
oakplots2010.modlist1 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2010),
      m2 <- lm(tot_bay.log ~ twi15m.log, data = oakplots.2010),
      m3 <- lm(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log, data = oakplots.2010),
      m6 <- lm(candens15m ~ tot_bay.log, data = oakplots.2010),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2010),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2010),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2010)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
E3 <- resid(m3)
plot(E3 ~ tot_bay.log, data = na.omit(oakplots.2010))
## strong negative relationship between UMCA abundance & residuals
## perhaps due to the zero values
plot(update(m3, .~., data = filter(na.omit(oakplots.2009), tot_lfct.log > 0)))
plot(E3 ~ dys_abv25ds, data = na.omit(oakplots.2010))
plot(E3 ~ twi15m.log, data = na.omit(oakplots.2010))

plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")

oakplots2010.modlist2 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2010),
      m2 <- lm(umca_ba.log ~ twi15m.log, data = oakplots.2010),
      m3 <- lm(tot_lfct.log ~ umca_ba.log + hrs_abv25ds + twi15m.log, data = oakplots.2010),
      m6 <- lm(candens15m ~ umca_ba.log, data = oakplots.2010),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2010),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2010),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2010)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")

sem.fit(oakplots2010.modlist1, oakplots.2010)
sem.fit(oakplots2010.modlist2, oakplots.2010)

sem.coefs(oakplots2010.modlist1, oakplots.2010)
sem.coefs(oakplots2010.modlist1, oakplots.2010, standardize = 'scale')
sem.coefs(oakplots2010.modlist2, oakplots.2010)
sem.coefs(oakplots2010.modlist2, oakplots.2010, standardize = 'scale')


#+ path model 2011 ####
oakplots2011.modlist1 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2011),
      m2 <- lm(tot_bay.log ~ twi15m.log, data = oakplots.2011),
      m3 <- lm(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log, data = oakplots.2011),
      m6 <- lm(candens15m ~ tot_bay.log, data = oakplots.2011),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2011),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2011),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2011)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
E3 <- resid(m3)
plot(E3 ~ tot_bay.log, data = na.omit(oakplots.2011))
## strong negative relationship between UMCA abundance & residuals
## perhaps due to the zero values
plot(update(m3, .~., data = filter(na.omit(oakplots.2009), tot_lfct.log > 0)))
plot(E3 ~ dys_abv25ds, data = na.omit(oakplots.2011))
plot(E3 ~ twi15m.log, data = na.omit(oakplots.2011))

plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


oakplots2011.modlist2 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2011),
      m2 <- lm(umca_ba.log ~ twi15m.log, data = oakplots.2011),
      m3 <- lm(tot_lfct.log ~ umca_ba.log + hrs_abv25ds + twi15m.log, data = oakplots.2011),
      m6 <- lm(candens15m ~ umca_ba.log, data = oakplots.2011),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2011),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2011),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2011)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


sem.fit(oakplots2011.modlist1, oakplots.2011)
sem.fit(oakplots2011.modlist2, oakplots.2011)

sem.coefs(oakplots2011.modlist1, oakplots.2011)
sem.coefs(oakplots2011.modlist1, oakplots.2011, standardize = 'scale')
sem.coefs(oakplots2011.modlist2, oakplots.2011)
sem.coefs(oakplots2011.modlist2, oakplots.2011, standardize = 'scale')


#+ path model 2012 ####
oakplots2012.modlist1 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2012),
      m2 <- lm(tot_bay.log ~ twi15m.log, data = oakplots.2012),
      m3 <- lm(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log, data = oakplots.2012),
      m6 <- lm(candens15m ~ tot_bay.log, data = oakplots.2012),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2012),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2012),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2012)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
E3 <- resid(m3)
plot(E3 ~ tot_bay.log, data = na.omit(oakplots.2012))
## strong negative relationship between UMCA abundance & residuals
## perhaps due to the zero values
plot(update(m3, .~., data = filter(na.omit(oakplots.2009), tot_lfct.log > 0)))
plot(E3 ~ dys_abv25ds, data = na.omit(oakplots.2012))
plot(E3 ~ twi15m.log, data = na.omit(oakplots.2012))

plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


oakplots2012.modlist2 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2012),
      m2 <- lm(umca_ba.log ~ twi15m.log, data = oakplots.2012),
      m3 <- lm(tot_lfct.log ~ umca_ba.log + hrs_abv25ds + twi15m.log, data = oakplots.2012),
      m6 <- lm(candens15m ~ umca_ba.log, data = oakplots.2012),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2012),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2012),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2012)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


sem.fit(oakplots2012.modlist1, oakplots.2012)
sem.fit(oakplots2012.modlist2, oakplots.2012)

sem.coefs(oakplots2012.modlist1, oakplots.2012)
sem.coefs(oakplots2012.modlist1, oakplots.2012, standardize = 'scale')
sem.coefs(oakplots2012.modlist2, oakplots.2012)
sem.coefs(oakplots2012.modlist2, oakplots.2012, standardize = 'scale')


#+ path model 2014 ####
oakplots2014.modlist1 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2014),
      m2 <- lm(tot_bay.log ~ twi15m.log, data = oakplots.2014),
      m3 <- lm(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + twi15m.log, data = oakplots.2014),
      m6 <- lm(candens15m ~ tot_bay.log, data = oakplots.2014),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2014),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2014),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2014)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
E3 <- resid(m3)
plot(E3 ~ tot_bay.log, data = na.omit(oakplots.2014))
## strong negative relationship between UMCA abundance & residuals
## perhaps due to the zero values
plot(update(m3, .~., data = filter(na.omit(oakplots.2009), tot_lfct.log > 0)))
plot(E3 ~ dys_abv25ds, data = na.omit(oakplots.2014))
plot(E3 ~ twi15m.log, data = na.omit(oakplots.2014))

plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


oakplots2014.modlist2 <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.2014),
      m2 <- lm(umca_ba.log ~ twi15m.log, data = oakplots.2014),
      m3 <- lm(tot_lfct.log ~ umca_ba.log + hrs_abv25ds + twi15m.log, data = oakplots.2014),
      m6 <- lm(candens15m ~ umca_ba.log, data = oakplots.2014),
      m7 <- lm(avg_hrs1422_wet_t_t1 ~ candens15m, data = oakplots.2014),
      m4 <- lm(hrs_abv25ds ~ candens15m + twi15m.log, data = oakplots.2014),
      m5 <- glm(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1 + H.2005 + candens15m, family = binomial(link = "logit"), data = oakplots.2014)
)

par(mfrow=c(2,2))
plot(m1, main = "Diversity ~")
plot(m2, main = "UMCA Abundance ~")
plot(m3, main = "Leaf Symptoms ~")
plot(m4, main = "Dry Season Temperature ~")
plot(m5, main = "Disease Prevalence ~")
plot(m6, main = "Canopy Density ~")
plot(m7, main = "Wet-Days Temperature ~")


sem.fit(oakplots2014.modlist1, oakplots.2014)
sem.fit(oakplots2014.modlist2, oakplots.2014)

sem.coefs(oakplots2014.modlist1, oakplots.2014)
sem.coefs(oakplots2014.modlist1, oakplots.2014, standardize = 'scale')
sem.coefs(oakplots2014.modlist2, oakplots.2014)
sem.coefs(oakplots2014.modlist2, oakplots.2014, standardize = 'scale')
#'


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
      
      lme(tot_lfct.log ~ tot_bay.log + avg_tmax_rs + ss2d.cm + twi15m + H.2005,
          random = ~1|sample_year/plotid, data = oakplots.sub,
          na.action = na.omit),
      
      glmer(oak.inf ~ tot_lfct.log + avg_tmax_rs + ss2d.cm + H.2005 + (1|sample_year) + (1|plotid),
            data = oakplots.sub, family = binomial(link = "logit"), 
            na.action = na.omit)
)

sem.missing.paths(oakplots.modlist, oakplots.sub)
sem.fit(oakplots.modlist, oakplots.sub, 
        corr.errors = "avg_tmax_rs ~~ ss2d.cm")
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


#+ subset and rescale variables ####
names(oakplots.sub)
oakplots.sub2 <- oakplots.sub %>% 
      select(plotid:umca_basal_area, oak.inf, tot_bay.log, umca_ba.log, tot_lfct.log, hrs_abv25ds, avg_hrsblw14_wet_tminus1, hrs1422_wet_tminus1, hrs1422_wet_tminus1.log, avg_hrs1422_wet_t_t1, avg_hrs1422_wet_t_t1.log, candens15m, H.2005, H.2014, twi15m, twi15m.log, -inf_oak_2plus_ct) %>% 
      filter(is.na(H.2005)==F, is.na(tot_bay)==F, is.na(avg_hrs1422_wet_t_t1)==F)
summary(oakplots.sub2)

## rescale variables
# oakplots.sub2$inf_bay.scl <- scale((oakplots.sub2$infected_bay_ct))
oakplots.sub2$tot_bay.scl <- scale(oakplots.sub2$tot_bay)
oakplots.sub2$tot_bay.log.scl <- scale(oakplots.sub2$tot_bay.log)
oakplots.sub2$umca_basal_area.scl <- scale(oakplots.sub2$umca_basal_area)
oakplots.sub2$umca_ba.log.scl <- scale(oakplots.sub2$umca_ba.log)
oakplots.sub2$tot_lfct.scl <- scale(oakplots.sub2$tot_lfct)
oakplots.sub2$tot_lfct.log.scl <- scale(oakplots.sub2$tot_lfct.log)
# oakplots.sub2$avg_lfct.scl <- scale(oakplots.sub2$avg_lfct)
oakplots.sub2$hrs_abv25ds.scl <- scale(oakplots.sub2$hrs_abv25ds)
oakplots.sub2$avg_hrsblw14_wet_tminus1.scl <- scale(oakplots.sub2$avg_hrsblw14_wet_tminus1)
oakplots.sub2$hrs1422_wet_tminus1.log.scl <- scale(oakplots.sub2$hrs1422_wet_tminus1.log)
# oakplots.sub2$avg_hrs1422_wet_t_t1.scl <- scale(oakplots.sub2$avg_hrs1422_wet_t_t1)
oakplots.sub2$avg_hrs1422_wet_t_t1.log.scl <- scale(oakplots.sub2$avg_hrs1422_wet_t_t1.log)

## static variables
oakplots.sub2$candens15m.scl <- scale(oakplots.sub2$candens15m)
oakplots.sub2$shannons2005.scl <- scale(oakplots.sub2$H.2005)
oakplots.sub2$shannons2014.scl <- scale(oakplots.sub2$H.2014)
oakplots.sub2$twi15m.scl <- scale(oakplots.sub2$twi15m)
oakplots.sub2$twi15m.log.scl <- scale(oakplots.sub2$twi15m.log)

summary(oakplots.sub2)


names(oakplots.sub2)
pairs(select(oakplots.sub2, ends_with("scl")),
      lower.panel = panel.cor, diag.panel = panel.hist,
      main = "Scaled and/or Transformed Variables")
#'

#+ path model pieces ####
m1 <- lme(tot_bay.log ~ twi15m, 
             random = ~1|sample_year,
             data = oakplots.sub2, na.action = na.omit)
summary(m1)
plot(m1, main = "Total # Bay Laurel ~ TWI + Diveristy Residuals")

m1a <- lme(tot_bay.log ~ twi15m.log, 
           random = ~1|sample_year,
           data = oakplots.sub2, na.action = na.omit)
summary(m1a)
plot(m1a, main = "Total # Bay Laurel ~ TWI Residuals")

m1scl <- lme(tot_bay.log.scl ~ twi15m.scl, 
          random = ~1|sample_year,
          data = oakplots.sub2, na.action = na.omit)
summary(m1scl)
plot(m1scl, main = "Total # Bay Laurel ~ TWI + Diveristy Residuals")

m1scl.a <- lme(tot_bay.log.scl ~ twi15m.log.scl, 
             random = ~1|sample_year,
             data = oakplots.sub2, na.action = na.omit)
summary(m1scl.a)
plot(m1scl.a, main = "Total # Bay Laurel ~ TWI + Diveristy Residuals")
m1 <- m1scl.a



m2 <- lme(H.2005 ~ twi15m.log.scl, 
          random = ~1|sample_year, data = oakplots.sub2,
          na.action = na.omit)
summary(m2)
plot(m2, main = "Diversity ~ TWI Residuals")

m2scl <- lme(shannons2005.scl ~ twi15m.log.scl, 
          random = ~1|sample_year, data = oakplots.sub2,
          na.action = na.omit)
summary(m2scl)
plot(m2scl, main = "Diversity ~ TWI Residuals")
m2 <- m2scl

m2a <- lme(H.2005 ~ twi15m.log, 
          random = ~1|sample_year, data = oakplots.sub2,
          na.action = na.omit)
summary(m2a)
plot(m2a, main = "Diversity ~ TWI Residuals")


m3 <- lme(avg_tmax_ds.scl ~ tot_bay.log.scl + twi15m.log.scl + shannons2005.scl+ rain_days.scl, 
             random = ~1|sample_year,
             data = oakplots.sub2, na.action = na.omit)
summary(m3)
plot(m3)

m3a <- lme(avg_tmax_ds ~ tot_bay.log + twi15m.log + H.2005 + rain_days_v, 
          random = ~1|sample_year,
          data = oakplots.sub2, na.action = na.omit)
summary(m3a)
plot(m3a)

# m3a <- lmer(avg_tmax_ds ~ tot_bay.log + twi15m + (1|sample_year),
#             data = oakplots.sub, na.action = na.omit)
# summary(m3a); AIC(m3a)
# plot(m3a)

m4 <- lme(tot_lfct.log.scl ~ tot_bay.log.scl + avg_tmax_ds.scl + rain_days.scl + twi15m.log.scl + shannons2005.scl, 
          random = ~1|sample_year, data = oakplots.sub2, 
          na.action = na.omit)
plot(m4)
summary(m4)
E2 <- resid(m4, type = "normalized")# equivalent of "standardized" in plot
F2 <- fitted(m4)
plot(F2, E2)
plot(E2 ~ tot_bay.log.scl, data = oakplots.sub2, ylab = "resids", xlab = "total # bay laurel")# this seems the most likely source of the pattern
plot(E2 ~ avg_tmax_ds.scl, data = oakplots.sub2, ylab = "resids", xlab = "temp")
plot(E2 ~ rain_days.scl, data = oakplots.sub2, ylab = "resids", xlab = "rain")
plot(E2 ~ twi15m.log.scl, data = oakplots.sub2, ylab = "resids", xlab = "twi")
plot(E2 ~ shannons2005.scl, data = oakplots.sub2, ylab = "resids", xlab = "Shannon's Index")
boxplot(E2 ~ sample_year, data = oakplots.sub2, ylab = "resids", xlab = "year")
## These plots indicate that the trend in the residuals is not a result of the predictors, at least none of them are contributing the substantial amount abserved.

m4a <- lme(tot_lfct.log ~ tot_bay.log + avg_tmax_ds + rain_days_v + twi15m.log + H.2005, 
           random = ~1|sample_year/plotid, data = oakplots.sub2, 
          na.action = na.omit)
# I believe this is allowing the intercept to vary among groups within sample year

m4a <- lme(tot_lfct.log ~ tot_bay.log + avg_tmax_ds + rain_days_v + twi15m.log + H.2005, 
           random = ~1|sample_year, data = oakplots.sub2, 
           na.action = na.omit)
plot(m4a)
summary(m4a)

m4b <- lmer(tot_lfct.log.scl ~ tot_bay.log.scl + avg_tmax_ds.scl + rain_days.scl + twi15m.log.scl + shannons2005.scl + (1|sample_year) + (1|plotid), 
            data = oakplots.sub2)
# this formulation is allowing the intercept to vary among each group
summary(m4b)
plot(m4b)
image(getME(m4b, name = "L"))
image((getME(m4b, name = "A")))

m5 <- glmer(oak.inf ~ tot_lfct.log.scl + tot_bay.log.scl + avg_tmax_ds.scl + shannons2005.scl + (1|sample_year),
            data = oakplots.sub2, family = binomial(link = "logit"), 
            na.action = na.omit)
summary(m5)
plot(m5)

m5a <- glmer(oak.inf ~ tot_lfct.log + tot_bay.log + avg_tmax_ds + H.2005 + (1|sample_year),
            data = oakplots.sub2, family = binomial(link = "logit"), 
            na.action = na.omit)
summary(m5a)
plot(m5a)

summary(oakplots.sub2)
# lfct: 0-8.8, tmax: 19.4-31.4, rain: 199.5-1529.1, H.2005: 0-1.6
# based on these ranges I'm supposing that the rain variable is throwing things off re: large eigenvalue
## After standardizing (scale) all the variables I no longer get an error when running this model.


# oakplots.sub$rain_tot_2d.dm <- oakplots.sub$rain_tot_2d.cm/10
# m5a <- glmer(oak.inf ~ tot_lfct.log + avg_tmax_ds + rain_tot_2d.dm + H.2005 + (1|sample_year),
#              data = oakplots.sub, family = binomial(link = "logit"), 
#              na.action = na.omit)
# # still fails to converge at set tolerance level, but avoids large eigenvalue
# summary(m5a)# essentially the same fixed effects values taking into account rainfall transformation
# plot(m5a)# same as for model m5
# 
# # add an observation-level (plot) random effect to account for additional variation...
# m5b <- glmer(oak.inf ~ tot_lfct.log + avg_tmax_ds + rain_tot_2d.cm + H.2005 + (1|sample_year) + (1|plotid),
#              data = oakplots.sub, family = binomial(link = "logit"), 
#              na.action = na.omit)
# summary(m5b)# notable changes to estimates and p-values; lots of variance attributed to plot-level random effect
# plot(m5b)# this is...artsy? still with the nonlinear negative trend
# 
# image(getME(m5a, name = "L"))
# image(getME(m5b, name = "L"))

# Fit Path Model
oakplots.modlist <- list(m1,m2,m3,m4,m5)
sem.lavaan(oakplots.modlist, oakplots.sub2)
sem.fit(oakplots.modlist, oakplots.sub2)
## p-value is significant; too many significant missing paths
## add correlated error structure
# sem.fit(oakplots.modlist, oakplots.sub2, corr.errors = 
#               c("avg_tmax_ds.scl ~~ rain_days.scl"))
# sem.model.fits(oakplots.modlist)
sem.coefs(oakplots.modlist, oakplots.sub2)


oakplots.modlist2 <- list(m1a,m2a,m3a,m4a,m5a)
sem.fit(oakplots.modlist2, oakplots.sub2)
# sem.fit(oakplots.modlist2, oakplots.sub2, corr.errors = 
#               c("avg_tmax_ds ~~ rain_days_v"))
sem.model.fits(oakplots.modlist2)

sem.coefs(oakplots.modlist2, oakplots.sub2)
sem.coefs(oakplots.modlist2, oakplots.sub2, standardize = "scale")

coef.table <- (sem.coefs(oakplots.modlist2, oakplots.sub2))

#+ revised disease prevalence repeated measures path model ####
disease.prevalence.modlist <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.sub2),
      
      m2 <- lme(tot_bay.log ~ twi15m.log, random = ~1|sample_year, data = oakplots.sub2),
      
      m3 <- lme(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + H.2005 + twi15m.log, random = ~1|sample_year, data = oakplots.sub2),
      
      m4 <- lme(hrs_abv25ds ~ avg_hrs1422_wet_t_t1.log + tot_bay.log + twi15m.log, random = ~1|sample_year, data = oakplots.sub2),
      
      m6 <- lme(avg_hrs1422_wet_t_t1.log ~ tot_bay.log + H.2005, random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1.log + hrs_abv25ds + tot_bay.log + H.2005 + (1|sample_year), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)

disease.prevalence.modlist.scl <- list(
      m1 <- lm(shannons2005.scl ~ twi15m.log.scl, data = oakplots.sub2),
      
      m2 <- lme(tot_bay.log.scl~ twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m3 <- lme(tot_lfct.log.scl ~ tot_bay.log.scl + hrs_abv25ds.scl + shannons2005.scl + twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m4 <- lme(hrs_abv25ds.scl ~ avg_hrs1422_wet_t_t1.scl.log + tot_bay.log.scl + twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m6 <- lme(avg_hrs1422_wet_t_t1.scl.log ~ tot_bay.log.scl + shannons2005.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log.scl + avg_hrs1422_wet_t_t1.scl.log + hrs_abv25ds.scl + tot_bay.log.scl + shannons2005.scl + (1|sample_year), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)

sem.fit(disease.prevalence.modlist, oakplots.sub2)
sem.fit(disease.prevalence.modlist.scl, oakplots.sub2
        #corr.errors = "hrs_abv25ds.scl ~~ avg_hrs1422_wet_t_t1.scl"
)

sem.coefs(disease.prevalence.modlist, oakplots.sub2)
sem.coefs(disease.prevalence.modlist, oakplots.sub2, standardize = "scale")
sem.coefs(disease.prevalence.modlist.scl, oakplots.sub2)


disease.prevalence.ws.modlist <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.sub2),
      
      m2 <- lme(tot_bay.log ~ twi15m.log, random = ~1|sample_year, data = oakplots.sub2),
      
      m3 <- lme(tot_lfct.log ~ tot_bay.log + H.2005 + twi15m.log, random = ~1|sample_year, data = oakplots.sub2),
      
      # m4 <- lme(hrs_abv25ds ~ avg_hrs1422_wet_t_t1.log + tot_bay.log + twi15m.log, random = ~1|sample_year, data = oakplots.sub2),
      
      m6 <- lme(avg_hrs1422_wet_t_t1.log ~ tot_bay.log + H.2005, random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log + avg_hrs1422_wet_t_t1.log + tot_bay.log + H.2005 + (1|sample_year), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)
sem.fit(disease.prevalence.ws.modlist, oakplots.sub2)
sem.coefs(disease.prevalence.ws.modlist, oakplots.sub2, standardize = "scale")

disease.prevalence.ws.modlist.scl <- list(
      m1 <- lm(shannons2005.scl ~ twi15m.log.scl, data = oakplots.sub2),
      
      m2 <- lme(tot_bay.log.scl~ twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m3 <- lme(tot_lfct.log.scl ~ tot_bay.log.scl + shannons2005.scl + twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      # m4 <- lme(hrs_abv25ds.scl ~ avg_hrs1422_wet_t_t1.scl.log + tot_bay.log.scl + twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m4 <- lme(avg_hrs1422_wet_t_t1.scl.log ~ tot_bay.log.scl + shannons2005.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log.scl + avg_hrs1422_wet_t_t1.scl.log +  tot_bay.log.scl + shannons2005.scl + (1|sample_year), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)


disease.prevalence.ds.modlist <- list(
      m1 <- lm(H.2005 ~ twi15m.log, data = oakplots.sub2),
      
      m2 <- lmer(tot_bay.log ~ twi15m.log + (1|sample_year) + (1|plotid), data = oakplots.sub2),
      
      m3 <- lmer(tot_lfct.log ~ tot_bay.log + hrs_abv25ds + H.2005 + twi15m.log + (1|sample_year) + (1|plotid), data = oakplots.sub2),
      
      m4 <- lmer(hrs_abv25ds ~ tot_bay.log + twi15m.log + (1|sample_year) + (1|plotid), data = oakplots.sub2),
      
      # m6 <- lme(avg_hrs1422_wet_t_t1..log ~ tot_bay.log., random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log + hrs_abv25ds + H.2005 + (1|sample_year) + (1|plotid), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)
sem.fit(disease.prevalence.ds.modlist, oakplots.sub2)
sem.coefs(disease.prevalence.ds.modlist, oakplots.sub2, standardize = "none")

## dry season temperature model year/plotid random effects ####
disease.prevalence.ds.modlist.scl <- list(
      m1 <- lm(shannons2005.scl ~ twi15m.log.scl, data = oakplots.sub2),
      
      m2 <- lmer(tot_bay.log.scl~ twi15m.log.scl + (1|sample_year) + (1|plotid), data = oakplots.sub2),
      
      m3 <- lmer(tot_lfct.log.scl ~ tot_bay.log.scl + hrs_abv25ds.scl + shannons2005.scl + twi15m.log.scl + (1|sample_year) + (1|plotid), data = oakplots.sub2),
      
      m4 <- lmer(hrs_abv25ds.scl ~ tot_bay.log.scl + twi15m.log.scl + (1|sample_year) + (1|plotid), data = oakplots.sub2),
      
      # m6 <- lme(avg_hrs1422_wet_t_t1.scl.log ~ tot_bay.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log.scl + hrs_abv25ds.scl + shannons2005.scl + (1|sample_year) + (1|plotid), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)
sem.fit(disease.prevalence.ds.modlist.scl, oakplots.sub2)
sem.coefs(disease.prevalence.ds.modlist.scl, oakplots.sub2)

dp.ds.lme.modlist.scl <- list(
      m1 <- lm(shannons2005.scl ~ twi15m.log.scl, data = oakplots.sub2),
      
      m2 <- lme(tot_bay.log.scl~ twi15m.log.scl, random = ~1|sample_year/plotid, data = oakplots.sub2),
      
      m3 <- lme(tot_lfct.log.scl ~ tot_bay.log.scl + hrs_abv25ds.scl + shannons2005.scl + twi15m.log.scl, random = ~1|sample_year/plotid, data = oakplots.sub2),
      
      m4 <- lme(hrs_abv25ds.scl ~ tot_bay.log.scl + twi15m.log.scl, random = ~1|sample_year/plotid, data = oakplots.sub2),
      
      # m6 <- lme(avg_hrs1422_wet_t_t1.scl.log ~ tot_bay.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log.scl + hrs_abv25ds.scl + shannons2005.scl + (1|sample_year) + (1|plotid), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)
sem.fit(dp.ds.lme.modlist.scl, oakplots.sub2)
sem.coefs(dp.ds.lme.modlist.scl, oakplots.sub2)

## dry season temperature without bay laurel density ####
disease.prevalence.ds.rmbay.modlist.scl <- list(
      m1 <- lm(shannons2005.scl ~ twi15m.log.scl, data = oakplots.sub2),
      
      # m2 <- lme(umca_ba.log.scl~ twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m3 <- lme(tot_lfct.log.scl ~ hrs_abv25ds.scl + shannons2005.scl + twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m4 <- lme(hrs_abv25ds.scl ~ twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      # m6 <- lme(avg_hrs1422_wet_t_t1.scl.log ~ umca_ba.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log.scl + hrs_abv25ds.scl + shannons2005.scl + (1|sample_year), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)

## model fitting ####
sem.fit(disease.prevalence.ws.modlist.scl, oakplots.sub2)
sem.fit(disease.prevalence.ds.modlist.scl, oakplots.sub2)
sem.model.fits(disease.prevalence.ds.modlist.scl)
sem.fit(disease.prevalence.ds.rmbay.modlist.scl, oakplots.sub2)
sem.coefs(disease.prevalence.ds.modlist.scl, oakplots.sub2)
sem.coefs(disease.prevalence.ds.rmbay.modlist.scl, oakplots.sub2)
sem.coefs(disease.prevalence.ws.modlist.scl, oakplots.sub2)


## dry season temperature with year as fixed effect ####
oakplots.sub2$sample_year.scl <- scale(oakplots.sub2$sample_year)

disease.prevalence.ds.yearfixed.modlist.scl <- list(
      m1 <- lm(shannons2005.scl ~ twi15m.log.scl, data = oakplots.sub2),
      
      m2 <- lme(tot_bay.log.scl~ twi15m.log.scl + sample_year.scl, random = ~1|plotid, data = oakplots.sub2),
      
      m3 <- lme(tot_lfct.log.scl ~ tot_bay.log.scl + hrs_abv25ds.scl + shannons2005.scl + twi15m.log.scl  + sample_year.scl, random = ~1|plotid, data = oakplots.sub2),
      
      m4 <- lmer(hrs_abv25ds.scl ~ tot_bay.log.scl + twi15m.log.scl + sample_year.scl + (1|plotid), data = oakplots.sub2),
      
      # m6 <- lme(avg_hrs1422_wet_t_t1.scl.log ~ tot_bay.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log.scl + hrs_abv25ds.scl + shannons2005.scl + sample_year.scl + (1|plotid), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)
sem.fit(disease.prevalence.ds.yearfixed.modlist.scl, oakplots.sub2)


## dry season canopy density model ####
disease.prevalence.ds.candens.modlist.scl <- list(
      m1 <- lm(shannons2005.scl ~ twi15m.log.scl, data = oakplots.sub2),
      
      # m2 <- lme(candens15m.scl~ twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m3 <- lme(tot_lfct.log.scl ~ hrs_abv25ds.scl + candens15m.scl + shannons2005.scl + twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m4 <- lme(hrs_abv25ds.scl ~ candens15m.scl + twi15m.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      # m6 <- lme(avg_hrs1422_wet_t_t1.scl.log ~ umca_ba.log.scl, random = ~1|sample_year, data = oakplots.sub2),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log.scl + hrs_abv25ds.scl + shannons2005.scl + candens15m.scl + (1|sample_year), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)
sem.fit(disease.prevalence.ds.candens.modlist.scl, oakplots.sub2)
sem.model.fits(disease.prevalence.ds.candens.modlist.scl)
(coef.table <- sem.coefs(disease.prevalence.ds.candens.modlist.scl, oakplots.sub2))
sem.plot(disease.prevalence.ds.candens.modlist.scl, oakplots.sub2)
sem.plot(coef.table)

## simplified path model: TWI-->dry-season temp-->slc-->oak infection ####
dp.twi.ds.slc.modlist <- list(
      lmer(hrs_abv25ds.scl ~ twi15m.log.scl + (1|sample_year) + (1|plotid), data = oakplots.sub2),
      lmer(tot_lfct.log.scl ~ hrs_abv25ds.scl + twi15m.log.scl + (1|sample_year) + (1|plotid), data = oakplots.sub2),
      glmer(oak.inf ~ tot_lfct.log.scl + (1|sample_year) + (1|plotid), family = binomial(link = "logit"), data = oakplots.sub2, na.action = na.omit)
)
sem.fit(dp.twi.ds.slc.modlist, oakplots.sub2)
sem.model.fits(dp.twi.ds.slc.modlist)


# Oak Stem-Level Model ####
names(stems)
# stems$plotid <- tolower(stems$plotid)
# stems$species <- tolower(stems$species)

oak_stems <- stems %>% 
      select(plotid:stem_status, canker, sod_killed, dbh_entry, location, year) %>%
      filter(stem_status == "Alive" | stem_status == "Dead", species != "lide", species != "unknown sp.", species != "umca", location == "In")
oak_stems <- droplevels(oak_stems)
unique(oak_stems$species)

head(oak_stems,20)
summary(oak_stems)
tail(filter(oak_stems, is.na(initial_dbh)),20)#
filter(oak_stems, is.na(initial_dbh))
oak_stems$infected <- ifelse(oak_stems$canker == 0, 0, 1)
names(oakplots.sub2)
summary(oakplots.sub2)

# join plot-level metrics to stem-level data ####
oak_stems2005_plus <- left_join(oak_stems, 
                       oakplots.sub2,
                       by = c("plotid", "year"="sample_year")) %>% 
      filter(is.na(tot_bay)==F, year >= 2005)
oak_stems2005_plus <- rename(oak_stems2005_plus, sample_year = year)
summary(oak_stems2005_plus)
# unique(filter(oak_stems2005_plus, is.na(tot_bay))$plotid)
length(unique(oak_stems2005_plus$tag))

# build SEM model pieces ####
oak.infection.modlist.scl <- list(
      m1 <- lme(shannons2005.scl ~ twi15m.log.scl, random = ~1|sample_year, data = oak_stems2005_plus),
      
      m2 <- lme(tot_bay.log.scl ~ twi15m.log.scl, random = ~1|sample_year, data = oak_stems2005_plus),
      
      m3 <- lmer(tot_lfct.log.scl ~ tot_bay.log.scl + hrs_abv25ds.scl + shannons2005.scl + twi15m.log.scl + (1|sample_year), data = oak_stems2005_plus),
      
      m4 <- lme(hrs_abv25ds.scl ~ avg_hrs1422_wet_t_t1.scl.log + tot_bay.log.scl + twi15m.log.scl, random = ~1|sample_year, data = oak_stems2005_plus),
      
      m6 <- lme(avg_hrs1422_wet_t_t1.scl.log ~ tot_bay.log.scl, random = ~1|sample_year, data = oak_stems2005_plus),
      
      m5 <- glmer(oak.inf ~ tot_lfct.log.scl + avg_hrs1422_wet_t_t1.scl.log + hrs_abv25ds.scl + tot_bay.log.scl + shannons2005.scl + (1|sample_year), family = binomial(link = "logit"), data = oak_stems2005_plus, na.action = na.omit)
)



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

