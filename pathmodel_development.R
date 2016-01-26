#' # Path Analysis Presented at 2015 IALE World Congress
#' Final details for the IALE presentation are in the `iale2015_sem.Rmd` file. 

#' ## Path Analysis of Plot-Level Oak Infection 2005
#' Here I am using the `oak_sod` data set which will be the oak infection (binomial) as the terminal response variable.

#+ oak_sod path model 2005
# load data
load("pathmodel_data_20160112.RData")
library(plyr); library(dplyr)
pairs(select(oak_sod, elev_m, starts_with("hrs"), starts_with("avg_t"), contains("rain_tot")), lower.panel = panel.cor, diag.panel = panel.hist, main = "Oak-SOD Rainfall & Temperature Data")

library(ggm)
names(oak_sod)
oak_path1 <- DAG(
      twi15m ~ elev_m,
      os_rich05 ~ twi15m,
      hrs_14_20_rs ~ elev_m + d2c + avg_psi + twi15m,
      tot_lfct ~ hrs_14_20_rs + os_rich05 + d2c,
      inf_oak_ct ~ tot_lfct + os_rich05 + hrs_14_20_rs + d2c
)
isAcyclic(oak_path1)
bu_oak_path1 <- basiSet(oak_path1)
bu_oak_path1

oak_path2 <- DAG(
      twi15m ~ elev15m,
      os_rich ~ twi15m,
      hrs_blw10 ~ elev15m + d2c + avg_psi + twi15m,
      tot_lfct ~ hrs_blw10 + os_rich + d2c,
      inf_oak_ct ~ tot_lfct + os_rich + hrs_blw10 + d2c
)
isAcyclic(oak_path2)

oak_2005 <- filter(oak_sod, sample_year == 2005)
summary(oak_2005)
# oak_2005 <- select(oak_2005, plotid, uninf_oak_ct, inf_oak_ct, tot_lfct, hrs_14_20, hrs_blw10, hrs_abv25, avg_tmax, avg_tmin, rain05_2d, rain05_v, elev15m, twi15m, avg_psi, os_rich, d2c)
filter(oak_2005, is.na(oak_2005$tot_lfct))

oak_2005$tot_lfct[is.na(oak_2005$tot_lfct)] <- 0
filter(oak_2005, is.na(os_rich05))
# oak_2005 <- na.omit(oak_2005)
oak_2005$plotid <- as.factor(oak_2005$plotid)

oak2005.modlist <- list(
      lme(twi15m ~ elev_m, random = ~1|plotid, data = oak_2005, 
          na.action = na.omit),
      lme(os_rich05 ~ twi15m, random = ~1|plotid, data = oak_2005,
          na.action = na.omit),
      lme(I(hrs_14_20_rs/10) ~ elev_m + I(d2c/10000) + avg_psi + twi15m, random = ~1|plotid,
          data = oak_2005, na.action = na.omit),
      glmer(I(tot_lfct/10) ~ I(hrs_14_20_rs/10) + os_rich05 + I(d2c/10000) + (1|plotid), family = "poisson", data = oak_2005, na.action = na.omit),
      glmer(cbind(inf_oak_ct,uninf_oak_ct) ~ I(tot_lfct/10) + os_rich05 + I(hrs_14_20_rs/10) + I(d2c/10000) + (1|plotid), family = binomial(link = "logit"), data = oak_2005, na.action = na.omit) 
)

sem.basis.set(oak2005.modlist)

system.time(oak2005.fit <- sem.fit(oak2005.modlist, oak_2005))
oak2005.fit

sem.coefs(oak2005.modlist, oak_2005)

sem.model.fits(oak2005.modlist)


# This model makes Gaussian assumptions, which should be incorrect for these data
shipley.test(oak_path1, cov(select(oak_2005, elev15m, twi15m, os_rich, hrs_14_20, d2c, avg_psi, tot_lfct, inf_oak_ct)), n = 170)
# The p-value of 0.08 indicates that there is not a severe violation of the conditional independence assumptions with the application of a Gaussian model

shipley.test(oak_path2, cov(select(oak_2005, elev15m, twi15m, os_rich, hrs_blw10, d2c, avg_psi, tot_lfct, inf_oak_ct)), n = 170)
# The p-value of .04 indicates that the conditional independence assumptions under the application of a Gaussian model are questionable.


#' Violation or questioning of conditional independence supports addressing the complexity of the data using different model distributions to deal with the data.

#+ assess conditional independence claims 2005
bu_oak_path1 # 15 independence claims means 15 models
library(nlme); library(lme4); library(MuMIn)

f1 <- lme(avg_psi ~ d2c, data = oak_2005, random = ~1|plotid)
summary(f1) # p = 0.5091
f2 <- lme(avg_psi ~ elev15m, data = oak_2005, random = ~1|plotid)
summary(f2) # p = 0.2013
f3 <- lme(avg_psi ~ twi15m + elev15m, data = oak_2005, random = ~1|plotid)
summary(f3) # twi15m; p = 0.0959
f4 <- lme(os_rich ~ avg_psi + twi15m, data = oak_2005, random = ~1|plotid)
summary(f4) # avg_psi; p = 0.1712
f5 <- lm(log1p(tot_lfct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2005)
# f5 <- lme(log1p(tot_lfct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2005, random = ~1|plotid)
# f5 <- glmer(tot_lfct ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + (1|plotid), data = oak_2005, family = "poisson")
summary(f5) # avg_psi; p = 0.5993
# f6 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2005, family = "binomial")
# f6 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2005, family = "quasibinomial")
f6 <- glmer(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct) + (1|plotid), data = oak_2005, family = "binomial")
summary(f6) # avg_psi; p = 0.08970
f7 <- lm(elev15m ~ d2c, data = oak_2005)
# f7 <- lme(elev15m ~ d2c, data = oak_2005, random = ~1|plotid)
summary(f7) # p = 0.03184
f8 <- lm(log(twi15m) ~ d2c + elev15m, data = oak_2005)
# f8 <- lme(log(twi15m) ~ d2c + elev15m, data = oak_2005, random = ~1|plotid)
summary(f8) # d2c; p = 0.3693
f9 <- lm(os_rich ~ d2c + twi15m, data = oak_2005)
# f9 <- lme(os_rich ~ d2c + twi15m, data = oak_2005, random = ~1|plotid)
summary(f9) # d2c; p = 0.9370
f10 <- lm(os_rich ~ elev15m + twi15m, data = oak_2005)
# f10 <- lme(os_rich ~ elev15m + twi15m, data = oak_2005, random = ~1|plotid)
summary(f10) # elev15m; p = 0.19174

# f11 <- lm(tot_lfct ~ scale(elev15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2005)
# f11 <- lm(log1p(tot_lfct) ~ scale(elev15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2005)
f11 <- glm(tot_lfct ~ scale(elev15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2005, family = "poisson")
chat <- deviance(f11)/df.residual(f11)
library(MuMIn)
x.quasipoisson<-function(...){
      res<-quasipoisson(...)
      res$aic<-poisson(...)$aic
      res
}
QAIC(update(f11, family = x.quasipoisson), chat = chat)
f11 <- glm(tot_lfct ~ scale(elev15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2005, family = "quasipoisson")
# f11 <- lme(log1p(tot_lfct) ~ scale(elev15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), random = ~1|plotid, data = oak_2005)
# f11 <- glmer(tot_lfct ~ scale(elev15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + (1|plotid), data = oak_2005, family = "poisson")
summary(f11) # elev15m; p = 0.360630

# f12 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(elev15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2005, family = "binomial")
# chat <- deviance(f12)/df.residual(f12)
# QAIC(update(f12, family = x.quasibinomial), chat = chat)
# f12 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(elev15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2005, family = "quasibinomial")
f12 <- glmer(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(elev15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct) + (1|plotid), data = oak_2005, family = "binomial")
summary(f12) # elevation; p = 0.44861

f13 <- lm(log1p(tot_lfct) ~ scale(twi15m) + scale(d2c) + scale(elev15m) + scale(hrs_14_20) + scale(os_rich), data = oak_2005)
# f13 <- lme(log1p(tot_lfct) ~ scale(tmi) + scale(d2c) + scale(elevation) + scale(hrs_14_20) + scale(os_rich), data = oak_2005, random = ~1|plotid)
# f13 <- glmer(tot_lfct ~ scale(tmi) + scale(d2c) + scale(elevation) + scale(hrs_14_20) + scale(os_rich) + (1|plotid), data = oak_2005, family = "poisson")
summary(f13) # twi15m; p = 0.801

# f14 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tmi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2005, family = "binomial")
# chat <- deviance(f14)/df.residual(f14)
# QAIC(update(f14, family = x.quasibinomial), chat = chat)
# f14 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tmi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2005, family = "quasibinomial")

f14 <- glmer(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(twi15m) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct) + (1|plotid), data = oak_2005, family = "binomial")
summary(f14) # twi15m; p = 0.61177

f15 <- lm(os_rich ~ hrs_14_20 + avg_psi + d2c + elev15m + twi15m, data = oak_2005)
# f15 <- lme(os_rich ~ hrs_14_20 + avg_psi + d2c + elevation + tmi, data = oak_2005, random = ~1|plotid)
summary(f15) # hrs_14_20; p = 0.08468
#'
#'
#+ Calculate C-statistic for 2005 data model
pvalues <- c(summary(f1)$tTable[2,5], summary(f2)$tTable[2,5], summary(f3)$tTable[2,5], summary(f4)$tTable[2,5], summary(f5)$coefficients[2,4], summary(f6)$coefficients[2,4], summary(f7)$tTable[2,5], summary(f8)$tTable[2,5], summary(f9)$tTable[2,5], summary(f10)$tTable[2,5], summary(f11)$coefficients[2,4], summary(f12)$coefficients[2,4], summary(f13)$coefficients[2,4], summary(f14)$coefficients[2,4], summary(f15)$tTable[2,5])

cstat <- -2 * sum(log(pvalues))
cstat
#'
#' Given the C-statistic of 26.9 and that 2k degrees of freedom would be 30 I think that this meets the criteria for supporting the conditional independence claims. So, this supports my path model, now I need to run the models to produce estimates for the path coefficients. This is the model set that I need to get the estimated path coefficients:
      
#'      1. os_rich ~ twi15m,
#'      2. hrs_14_20 ~ elev15m + d2c + avg_psi + twi15m,
#'      3. tot_lfct ~ hrs_14_20 + os_rich + d2c,
#'      4. infected oak ~ tot_lfct + os_rich + hrs_14_20 + d2c
#'      
#+ estimate path coefficients 2005
# f0 <- fitted(mlme1, level = 0)# Values from population model
# f1 <- fitted(mlme1, level = 1)# Within beach fitted values

# fit1 <- lm(tmi ~ scale(elevation), data = oak_2005)
fit1 <- lm(log(tmi) ~ scale(elevation), data = oak_2005) # best model
# fit1 <- lme(log(tmi) ~ scale(elevation), data = oak_2005, random = ~1|plotid, method = "ML")
summary(fit1)
par(mfrow = c(2,2))
plot(fit1)

# fit2 <- lm(hrs_14_20 ~ scale(elevation) + scale(d2c) + scale(avg_psi) + scale(tmi), data = oak_2005)
# fit2 <- lme(hrs_14_20 ~ scale(elevation) + scale(d2c) + scale(avg_psi) + scale(tmi), data = oak_2005, random = ~1|plotid) # best model
fit2 <- glmer(hrs_14_20 ~ scale(elevation) + scale(d2c) + scale(avg_psi) + scale(tmi) + (1|plotid), data = oak_2005, family = "poisson")
# fit2 <- glmer.nb(hrs_14_20 ~ scale(elevation) + scale(d2c) + scale(avg_psi) + scale(tmi) + (1|plotid), data = oak_2005)
summary(fit2)
plot(fit2)

# fit3 <- lm(tot_lfct ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c), data = oak_2005)
fit3 <- lm(log1p(tot_lfct) ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c), data = oak_2005)
fit3 <- glm(tot_lfct ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c), data = oak_2005, family = "poisson")
(chat <- deviance(fit3) / df.residual(fit3))
x.quasipoisson<-function(...){
      res<-quasipoisson(...)
      res$aic<-poisson(...)$aic
      res
}
library(MuMIn)
QAIC(update(fit3, family = x.quasipoisson), chat = chat) # 178.3
fit3 <- glm(tot_lfct ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c), data = oak_2005, family = "quasipoisson") # best model
# fit3 <- glmer(tot_lfct ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c) + (1|plotid), data = oak_2005, family = "poisson")
summary(fit3)

# fit4 <- lm((inf_oak_ct/(inf_oak_ct+uninf_oak_ct)) ~ scale(tot_lfct) + scale(os_rich) + scale(hrs_14_20) + scale(d2c), data = oak_2005)
# par(mfrow = c(2,2))
# plot(fit4)

# fit4 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tot_lfct) + scale(os_rich) + scale(hrs_14_20) + scale(d2c), data = oak_2005, family = "binomial")
# (chat <- deviance(fit4) / df.residual(fit4))
# x.quasibinomial<-function(...){
#   res<-quasibinomial(...)
#   res$aic<-binomial(...)$aic
#   res
# }
# library(MuMIn)
# QAIC(update(fit4, family = x.quasibinomial), chat = chat)
# fit4 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tot_lfct) + scale(os_rich) + scale(hrs_14_20) + scale(d2c), data = oak_2005, family = "quasibinomial")
# summary(fit4)

fit4 <- glmer(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tot_lfct) + scale(os_rich) + scale(d2c) + (1|plotid), data = oak_2005, family = "binomial") # best model
summary(fit4)

fit5 <- lm(os_rich ~ tmi, data = oak_2005); AIC(fit5)
# fit5 <- lme(os_rich ~ tmi, data = oak_2005, random = ~1|plotid)
summary(fit5)

# Trying to plot the fitted model by each group...
# f0 <- fitted(fit4, level = 0)
# f1 <- fitted(fit4, level = 1)
# 
# f0[I]
# 
# I <- order(oak_2005$tot_lfct) ; lfct <- sort(oak_2005$tot_lfct)
# plot(lfct, f0[I], lwd = 4, type = "l",
#      ylab = "oak infection", xlab = "total leaf count")
# for(i in 1:length(f1)){
#       x1 <- oak_2005$hrs_14_20[oak_2005$plotid == i]
#       y1 <- f1[oak_2005$hrs_14_20 == i]
#       K <- order(x1)
#       lines(sort(x1), y1[K], col = "red")
# }
# text(oak_2005$hrs_14_20, oak_2005$tot_lfct, oak_2005$plotid, cex = 0.9)


#' ## Path Analysis Oak Infection - 2014
#' Now I am going to fit the same path model to data from 2011.

#+ oak_sod path model 2014
oak_path1 <- DAG(
      tmi ~ elevation,
      os_rich ~ tmi,
      hrs_14_20 ~ elevation + d2c + avg_psi + tmi,
      tot_lfct ~ hrs_14_20 + os_rich + d2c,
      inf_oak_ct ~ tot_lfct + os_rich + hrs_14_20 + d2c
)

oak_2014 <- filter(oak_sod, year == 2014)
oak_2014 %>% select(plotid, tagged_rich, untagged_rich, os_rich) %>%
      filter(is.na(os_rich))
summary(oak_2014)
oak_2014 <- select(oak_2014, plotid, uninf_oak_ct, inf_oak_ct, tot_lfct, hrs_14_20, hrs_blw10, hrs_abv25, elevation, tmi, avg_psi, os_rich, d2c)
oak_2014$tot_lfct[is.na(oak_2014$tot_lfct)] <- 0
oak_2014 <- na.omit(oak_2014)
oak_2014$plotid <- as.factor(oak_2014$plotid)

shipley.test(oak_path1, cov(select(oak_2014, elevation, tmi, os_rich, hrs_14_20, d2c, avg_psi, tot_lfct, inf_oak_ct)), n = 168)
# Strongly significant p-value

# shipley.test(DAG(
#       tmi ~ elevation,
#       os_rich ~ tmi,
#       hrs_blw10 ~ elevation + d2c + avg_psi + tmi,
#       tot_lfct ~ hrs_blw10 + os_rich + d2c,
#       inf_oak_ct ~ tot_lfct + os_rich + hrs_blw10 + d2c), 
#       cov(select(oak_2014, elevation, tmi, os_rich, hrs_blw10, d2c, avg_psi, tot_lfct, inf_oak_ct)), 
#       n = 168)
#'
#'
#+ assess conditional independence claims 2014
f1 <- lme(avg_psi ~ d2c, data = oak_2014, random = ~1|plotid)
summary(f1) # p = 0.591
f2 <- lme(avg_psi ~ elevation, data = oak_2014, random = ~1|plotid)
summary(f2) # p = 0.2171
f3 <- lme(avg_psi ~ tmi + elevation, data = oak_2014, random = ~1|plotid)
summary(f3) # tmi; p = 0.1889
f4 <- lm(os_rich ~ avg_psi + tmi, data = oak_2014)
# f4 <- lme(os_rich ~ avg_psi + tmi, data = oak_2014, random = ~1|plotid)
summary(f4) # avg_psi; p = 0.156
f5 <- lm(log1p(tot_lfct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2014)
# f5 <- lme(log1p(tot_lfct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2014, random = ~1|plotid)
# f5 <- glmer(tot_lfct ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + (1|plotid), data = oak_2014, family = "poisson")
summary(f5) # avg_psi; p = 0.265

f6 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2014, family = "binomial")
chat <- deviance(f6)/df.residual(f6)
QAIC(update(f6, family = x.quasibinomial), chat = chat)
f6 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2014, family = "quasibinomial")
# f6 <- glmer(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(avg_psi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct) + (1|plotid), data = oak_2014, family = "binomial")
summary(f6) # avg_psi; p = 0.000722

f7 <- lm(log(elevation) ~ d2c, data = oak_2014)
# f7 <- lme(log(elevation) ~ d2c, data = oak_2014, random = ~1|plotid)
summary(f7) # p = 0.0101

f8 <- lm(log(tmi) ~ d2c + elevation, data = oak_2014)
# f8 <- lme(log(tmi) ~ d2c + elevation, data = oak_2014, random = ~1|plotid)
summary(f8) # d2c; p = 0.65795

f9 <- lm(log1p(os_rich) ~ d2c + tmi, data = oak_2014)
# f9 <- lme(log1p(os_rich) ~ d2c + tmi, data = oak_2014, random = ~1|plotid)
summary(f9) # d2c; p = 0.73

f10 <- lm(log1p(os_rich) ~ elevation + tmi, data = oak_2014)
# f10 <- lme(log1p(os_rich) ~ elevation + tmi, data = oak_2014, random = ~1|plotid)
summary(f10) # elevation; p = 0.126

# f11 <- lm(tot_lfct ~ scale(elevation) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2014)
f11 <- lm(log1p(tot_lfct) ~ scale(elevation) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2014)
f11 <- glm(tot_lfct ~ scale(elevation) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2014, family = "poisson")
chat <- deviance(f11)/df.residual(f11)
QAIC(update(f11, family = x.quasipoisson), chat = chat)
f11 <- glm((tot_lfct) ~ scale(elevation) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), data = oak_2014, family = "quasipoisson")
# f11 <- lme(log1p(tot_lfct) ~ scale(elevation) + scale(d2c) + scale(hrs_14_20) + scale(os_rich), random = ~1|plotid, data = oak_2014)
# f11 <- glmer(tot_lfct ~ scale(elevation) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + (1|plotid), data = oak_2014, family = "poisson")
summary(f11) # elevation; p = 0.0170

f12 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(elevation) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2014, family = "binomial")
chat <- deviance(f12)/df.residual(f12)
QAIC(update(f12, family = x.quasibinomial), chat = chat)
f12 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(elevation) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2014, family = "quasibinomial")
# f12 <- glmer(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(elevation) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct) + (1|plotid), data = oak_2014, family = "binomial")
summary(f12) # elevation; p = 0.00627

f13 <- lm(log1p(tot_lfct) ~ scale(tmi) + scale(d2c) + scale(elevation) + scale(hrs_14_20) + scale(os_rich), data = oak_2014)
f13 <- glm(tot_lfct ~ scale(tmi) + scale(d2c) + scale(elevation) + scale(hrs_14_20) + scale(os_rich), data = oak_2014, family = "poisson")
chat <- deviance(f13)/df.residual(f13)
QAIC(update(f13, family = x.quasipoisson), chat = chat)
f13 <- glm(tot_lfct ~ scale(tmi) + scale(d2c) + scale(elevation) + scale(hrs_14_20) + scale(os_rich), data = oak_2014, family = "quasipoisson")

# f13 <- lme(log1p(tot_lfct) ~ scale(tmi) + scale(d2c) + scale(elevation) + scale(hrs_14_20) + scale(os_rich), data = oak_2014, random = ~1|plotid)
# f13 <- glmer(tot_lfct ~ scale(tmi) + scale(d2c) + scale(elevation) + scale(hrs_14_20) + scale(os_rich) + (1|plotid), data = oak_2014, family = "poisson")
summary(f13) # tmi; p = 0.2562

f14 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tmi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2014, family = "binomial")
chat <- deviance(f14)/df.residual(f14)
QAIC(update(f14, family = x.quasibinomial), chat = chat)
f14 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tmi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct), data = oak_2014, family = "quasibinomial")
# f14 <- glmer(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tmi) + scale(d2c) + scale(hrs_14_20) + scale(os_rich) + scale(tot_lfct) + (1|plotid), data = oak_2014, family = "binomial")
summary(f14) # tmi; p = 0.7633

f15 <- lm(log1p(os_rich) ~ hrs_14_20 + avg_psi + d2c + elevation + tmi, data = oak_2014)
# f15 <- lme(log1p(os_rich) ~ hrs_14_20 + avg_psi + d2c + elevation + tmi, data = oak_2014, random = ~1|plotid)
summary(f15) # hrs_14_20; p = 0.1039

# Calculate C-statistic
pvalues <- c(summary(f1)$tTable[2,5], summary(f2)$tTable[2,5], summary(f3)$tTable[2,5], summary(f4)$coefficients[2,4], summary(f5)$coefficients[2,4], summary(f6)$coefficients[2,4], summary(f7)$coefficients[2,4], summary(f8)$coefficients[2,4], summary(f9)$coefficients[2,4], summary(f10)$coefficients[2,4], summary(f11)$coefficients[2,4], summary(f12)$coefficients[2,4], summary(f13)$coefficients[2,4], summary(f14)$coefficients[2,4], summary(f15)$coefficients[2,4])

cstat <- -2 * sum(log(pvalues))
cstat # 69.18
chisq.test(c(cstat,30))
qchisq(pvalues, df = 30)
chisq.test(cbind(pvalues, qchisq(pvalues, df = 30)))
#'
#' This C-statistic is much larger than the degrees of freedom, so I'm not confident in the conditional independencies of the structural model being met for the 2014 data. In fact, based on the C-test value from the earlier call to `shipleys.test` on these data this value is significantly different.
#'
#' The models testing for conditional independence indicate direct effects of PSI and elevation on oak infection, as well as a direct effect of elevation on symptomatic leaf count. This is indicated by the p-values showing strong statistical significance in the relationship between these variables. Below I have generated a new path model that accounts for the unassumed direct effects, which means they become assumed direct effects.
#'
#+ oak path 2014 v2
oak_path2 <- DAG(
tmi ~ elevation,
os_rich ~ tmi,
hrs_14_20 ~ elevation + d2c + avg_psi + tmi,
tot_lfct ~ hrs_14_20 + os_rich + d2c + elevation,
inf_oak_ct ~ tot_lfct + avg_psi + os_rich + hrs_14_20 + d2c + elevation
)
isAcyclic(oak_path2)
bu_oak_path2 <- basiSet(oak_path2)
bu_oak_path2

shipley.test(oak_path2, cov(select(oak_2014, elevation, tmi, os_rich, hrs_14_20, d2c, avg_psi, tot_lfct, inf_oak_ct)), n = 168)

plotGraph(oak_path2)
#'
#' This means commenting out models `f6` (oak ~ PSI), `f12` (oak ~ elevation), and `f11` (leaf count ~ elevation) in the `pvalues` summary.

#+ recalculate c-statistic 2014
pvalues <- c(
summary(f1)$tTable[2,5], summary(f2)$tTable[2,5], summary(f3)$tTable[2,5],
summary(f4)$coefficients[2,4], summary(f5)$coefficients[2,4], 
# summary(f6)$coefficients[2,4], 
summary(f7)$coefficients[2,4], summary(f8)$coefficients[2,4],
summary(f9)$coefficients[2,4], summary(f10)$coefficients[2,4],
# summary(f11)$coefficients[2,4], 
# summary(f12)$coefficients[2,4], 
summary(f13)$coefficients[2,4], summary(f14)$coefficients[2,4], summary(f15)$coefficients[2,4])

cstat <- -2 * sum(log(pvalues))
cstat # 36.56767
chisq.test(c(cstat,30))
#'
#' I think this would be evidence for wanting/needing to use the repeated measures approach across all the data to capture interannual variation in the relationships.
#' 
#' I am continuing with fitting the models to estimate the path coefficients for 2014 data based on this redrawn path model.
#' 
#+ estimate path coefficients 2014

# f0 <- fitted(mlme1, level = 0)# Values from population model
# f1 <- fitted(mlme1, level = 1)# Within beach fitted values

# fit1_2014 <- lm(tmi ~ scale(elevation), data = oak_2014)
fit1_2014 <- lm(log(tmi) ~ scale(elevation), data = oak_2014) # best model
# fit1_2014 <- lme(log(tmi) ~ scale(elevation), data = oak_2014, random = ~1|plotid)
summary(fit1_2014)
par(mfrow = c(2,2))
plot(fit1_2014)

# fit2_2014 <- lm(hrs_14_20 ~ scale(elevation) + scale(d2c) + scale(avg_psi) + scale(tmi), data = oak_2014)
# fit2_2014 <- lme(hrs_14_20 ~ scale(elevation) + scale(d2c) + scale(avg_psi) + scale(tmi), data = oak_2014, random = ~1|plotid) # best model
fit2_2014 <- glm(hrs_14_20 ~ scale(elevation) + scale(d2c) + scale(avg_psi) + scale(tmi), data = oak_2014, family = "poisson")
chat <- deviance(fit2_2014)/df.residual(fit2_2014)
QAIC(update(fit2_2014, family = x.quasipoisson), chat = chat)
fit2_2014 <- glm(hrs_14_20 ~ scale(elevation) + scale(d2c) + scale(avg_psi) + scale(tmi), data = oak_2014, family = "quasipoisson")
# fit2_2014 <- glmer(hrs_14_20 ~ scale(elevation) + scale(d2c) + scale(avg_psi) + scale(tmi) + (1|plotid), data = oak_2014, family = "poisson")
summary(fit2_2014)
plot(fit2_2014)

# fit3_2014 <- lm(tot_lfct ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c) + scale(elevation), data = oak_2014)
# fit3_2014 <- lm(log1p(tot_lfct) ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c) + scale(elevation), data = oak_2014)
fit3_2014 <- glm(tot_lfct ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c) + scale(elevation), data = oak_2014, family = "poisson")
(chat <- deviance(fit3_2014) / df.residual(fit3_2014))
x.quasipoisson<-function(...){
res<-quasipoisson(...)
res$aic<-poisson(...)$aic
res
}
library(MuMIn)
QAIC(update(fit3_2014, family = x.quasipoisson), chat = chat) # 178.3
fit3_2014 <- glm(tot_lfct ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c) + scale(elevation), data = oak_2014, family = "quasipoisson") # best model
# fit3_2014 <- glmer(tot_lfct ~ scale(hrs_14_20) + scale(os_rich) + scale(d2c) + scale(elevation) + (1|plotid), data = oak_2014, family = "poisson")
summary(fit3_2014)

# fit4_2014 <- lm((inf_oak_ct/(inf_oak_ct+uninf_oak_ct)) ~ scale(tot_lfct) + scale(os_rich) + scale(hrs_14_20) + scale(d2c) + scale(elevation) + scale(avg_psi), data = oak_2014)
par(mfrow = c(2,2))
plot(fit4_2014)

fit4_2014 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tot_lfct) + scale(os_rich) + scale(hrs_14_20) + scale(d2c) + scale(elevation) + scale(avg_psi), data = oak_2014, family = "binomial")
(chat <- deviance(fit4_2014) / df.residual(fit4_2014))
x.quasibinomial<-function(...){
res<-quasibinomial(...)
res$aic<-binomial(...)$aic
res
}
library(MuMIn)
QAIC(update(fit4_2014, family = x.quasibinomial), chat = chat)
fit4_2014 <- glm(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tot_lfct) + scale(os_rich) + scale(hrs_14_20) + scale(d2c) + scale(elevation) + scale(avg_psi), data = oak_2014, family = "quasibinomial")

# fit4_2014 <- glmer(cbind(inf_oak_ct, uninf_oak_ct) ~ scale(tot_lfct) + scale(os_rich) + scale(hrs_14_20) + scale(d2c) + (1|plotid), data = oak_2014, family = "binomial") # best model
summary(fit4_2014)

# fit5_2014 <- lm(os_rich ~ tmi, data = oak_2014); AIC(fit5_2014)
fit5_2014 <- lm(log(os_rich) ~ tmi, data = oak_2014); AIC(fit5_2014)
# fit5_2014 <- lme(os_rich ~ tmi, data = oak_2014, random = ~1|plotid)
# fit5_2014 <- lme(log(os_rich) ~ tmi, data = oak_2014, random = ~1|plotid)
summary(fit5_2014)


#' ## Importance Value
#' I also want to calculate some variables related to DBH, and in particular the "importance value" for each species. Importance value is calculated in the literature as:

#' `IV = 0.5 * (relative density / relative dominance)` where:
#'     
#'     - `Relative Density = # stems of species / total # of stems` and
#'     
#'     - `Relative Dominance = basal area of species / total basal area`
#'     
#+ calculate importance value ####
unique(untagged_dbh$species)
untagged_dbh
unique(stems$species)

tagged_dbh <- left_join(stems %>% 
filter(status == "Alive", location == "In") %>%
select(plotid, year, species, dbh) %>% 
group_by(plotid, species, year) %>% 
summarise(live_count = length(species), live_avg_dbh = mean(dbh), live_tot_dbh = sum(dbh)),
stems %>% 
filter(status == "Dead", location == "In") %>%
select(plotid, year, species, dbh) %>% 
group_by(plotid, species, year) %>% 
summarise(dead_count = length(species), dead_avg_dbh = mean(dbh), dead_tot_dbh = sum(dbh)),
by = c("plotid", "species", "year")
)
tagged_dbh <- as.tbl(tagged_dbh)
summary(tagged_dbh)

tagged_dbh
untagged_dbh

species_dbh <- rbind(tagged_dbh, untagged_dbh)
species_dbh
unique(species_dbh$species)
species_dbh %>% group_by(plotid, species, year)

stems_dbh <- species_dbh %>% group_by(plotid, year) %>% mutate(all_live_count = sum(live_count), all_live_dbh = sum(live_tot_dbh, na.rm = T), all_dead_count =  sum(dead_count, na.rm = T), all_dead_dbh = sum(dead_tot_dbh, na.rm = T))
stems_dbh
summary(stems_dbh)
filter(stems_dbh, all_live_dbh > 800)

stems_dbh <- stems_dbh %>% mutate(rel_dens = live_count / all_live_count, rel_dom = live_tot_dbh / all_live_dbh, imp_value = 0.5*(rel_dens + rel_dom))
summary(stems_dbh)
write.csv(stems_dbh, "analysis/data/species_dbh_impval_plot.csv")
#'
#+ pare down stems dataframe
stems %>% select(year, tag, species, status, notes) %>% filter(species == "umca") 
summary(stems) # Limited to "Alive/Dead"" status, "Missing" statuses excluded
unique(stems$species)
unique(stems$plotid)
summary(bay_laurel) # Still includes stems with "Missing" status
unique(bay_laurel$plot)

names(stems)
stemssub <- select(stems, plotid, year, cluster, tag, species, status, slc, canker, dbh, lide_lf_symp, sod_dead, location, dbh1, dbh2, delta_dbh, dbh2_1_ratio, inst_grwth_rate, doy, date)
summary(stemssub)
str(stemssub)

summary(oak_sod) # plots with oak stems
summary(umca_plots) # plots with bay laurel stems
summary(rain_data)
summary(temps_data)
summary(plots_env)
summary(veg_sub)
summary(richness)
summary(stems_dbh)

save.image("~/GitHub/superspreaders/data_201508.RData")

load("~/GitHub/superspreaders/data_201508.RData")
stemssub$cluster <- as.factor(stemssub$cluster)
stemssub$tag <- as.factor(stemssub$tag)
stemssub$canker <- as.factor(stemssub$canker)
stemssub$sod_dead <- as.factor(stemssub$sod_dead)
stemssub$date <- as.Date(stemssub$date)
str(stemssub)

