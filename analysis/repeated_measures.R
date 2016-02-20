---
#' Repeated measures analysis
---
#' Read in data for all stems that has been processed in "explore_data.R" ####
stems <- read.csv("analysis/data/stems_growth.csv")
str(stems)
summary(stems)

#' Load packages ####
library(plyr); library(dplyr); library(tidyr); library(ggplot2)

#' Examine stems data ####
names(stems)
stems <- stems %>%
      select(-X)# Remove the row names variable
stems$date <- as.Date(stems$date)
stems$year <- as.factor(stems$year)

#' Create bay laurel only data frame ####
umca <- as.tbl(stems) %>%
      select(plot, cluster, tag, species, status, slc, dbh, delta_dbh, dbh2_1_ratio, year, date) %>%
      filter(species == "UMCA", status == "Alive" | status == "Dead")
summary(umca)
umca <- droplevels(umca)

#' Develop plot level summaries data frame and join to stems data frame ####
umca <- left_join(umca, plot_umca <- umca %>%
      select(plot, cluster, tag, status, slc, year) %>% 
      group_by(plot, year) %>%
      filter(slc != "NA" & status == "Alive") %>%
      summarise(avg_slc = mean(slc), tot_slc = sum(slc), 
                live_dens = length(tag)), 
      by = c("plot", "year"))
summary(umca)# One NA in the summarized leaf count variables
filter(umca, is.na(avg_slc))
filter(umca, plot == "YAHNG02")# Looks like leaf count probably not done/recorded

#' Plot Level Repeated Measures Analysis Using Multi-level Mixed-Effects Regression ####

#' A couple of lines from Sarah's code to work from
# lmer.null <- lmer(log(cum.slc + 1) ~ 1 + (1 | year))
# lmer1 <- lmer(log(cum.slc + 1) ~ log(umca.ct + 1) + (1 | year))
# anova(lmer.null, lmer1)
# lmer2 <- lmer(log(cum.slc + 1) ~ log(umca.ct + 1) + (1 + log(umca.ct+1)|year), na.action=na.omit)
# anova(lmer1, lmer2)

#' Linear Mixed Effects Method ####

## Assumptions of linear mixed effects model:
##     1. b ~ N(0, D), random effects
##     2. E ~ N(0, var), errors "white noise"
##     3. These terms are independent in their distributions

## The use of random effects prevent systematic variation between groups from ending up in the residual errors ("white noise").

#' Starting with plot level data b/c simpler w/fewer nested groups
library(lme4)

#' Random Intercept Models ####
#' The random intercept model implies that there is one average curve (the population estimate) that may be shifted up or down by something normally distributed with variance d^2

#' Yi = Xi * B + Zi * bi + Ei
#'   bi ~ N(0, d^2), 
#'   Ei ~ N(0, covar)

#' The fitted lines for each group are parallel to the average curve

#' Simplest models, with only one random effect (random intercept)
m0 <- lmer(tot_slc ~ 1 + (1 | year), data = plot_umca)
m0a <- lmer(tot_slc ~ 1 + (1 | plot), data = plot_umca)
#' Add a second random effect where plot is nested within year
m0b <- lmer(tot_slc ~ 1 + (1 | year) + (1 | plot), data = plot_umca)
 
#' One predictor + random effect (intercept)
m1 <- lmer(tot_slc ~ live_dens + (1 | year), data = plot_umca)
m1a <- lmer(tot_slc ~ live_dens + (1 | plot), data = plot_umca)
#' One predictor + 2 random effects
m1b <- lmer(tot_slc ~ live_dens + (1 | year) + (1 | plot), data = plot_umca)

#' Random slope, random intercept model ####
#' The random effects term in the random slope, random intercept model accounts for the interaction between the predictor variable(s) and group in addition to the systematic variation between groups.
#'    - similar to linear regresstion with nominal & continuous variables, but requires a lot of parameters to be estimated if there are many groups (1 for each)
m2 <- lmer(tot_slc ~ live_dens + (1 + live_dens | year), data = plot_umca)
m2a <- lmer(tot_slc ~ live_dens + (1 | year) + (1 + live_dens | plot), data = plot_umca)
#' add the second random effect
m2b <- lmer(tot_slc ~ live_dens + (1 + live_dens | year) + (1 | plot), data = plot_umca)
m2c <- lmer(tot_slc ~ live_dens + (1 | year) + (1 + live_dens | plot), data = plot_umca)

#' Compare all models using Chi-square test of ANOVA ####
#' Refits model(s) with maximum likelihood instead of REML
anova(m0, m1, m2, m0a, m0b, m1a, m1b, m2a, m2b, m2c)
#' The random slope, random intercept model with the density:year interaction and the additional random effect is the best (lowest AIC), but the high values for the metrics may indicate some misspecification.

#' Plot the fitted models
par(mfrow = c(2,2))
plot(m0)
plot(m0a)
plot(m0b)
plot(m1)
plot(m1a)
plot(m1b)
plot(m2)
plot(m2a)
plot(m2b)
plot(m2c)
#' All of them show a lot of spread - expected because I know the distributions for leaf count and density aren't even close to normal
par(mfrow = c(2,2))
hist(plot_umca$tot_slc)
hist(plot_umca$live_dens)
#' Log transformation improves the distributions, but with some zeroes for symptomatic leaf count I need to add a constant. The `log1p` function computes log(1+x) accurately for |x| much less than one (automatically adds the constant 1) 
log_tot_slc <- log1p(plot_umca$tot_slc); hist(log_tot_slc)
log_live_dens <- log(plot_umca$live_dens); hist(log_live_dens)


#' Refit models with transformed variables ####
m0 <- lmer(log_tot_slc ~ 1 + (1 | year), data = plot_umca)
m0a <- lmer(log_tot_slc ~ 1 + (1 | plot), data = plot_umca)
m0b <- lmer(log_tot_slc ~ 1 + (1 | year) + (1 | plot), data = plot_umca)
m1 <- lmer(log_tot_slc ~ log_live_dens + (1 | year), data = plot_umca)
m1a <- lmer(log_tot_slc ~ log_live_dens + (1 | plot), data = plot_umca)
m1b <- lmer(log_tot_slc ~ log_live_dens + (1 | year) + (1 | plot), data = plot_umca)
m2 <- lmer(log_tot_slc ~ log_live_dens + (1 + log_live_dens | year), data = plot_umca)
m2a <- lmer(log_tot_slc ~ log_live_dens + (1 + log_live_dens | plot), data = plot_umca)
m2b <- lmer(log_tot_slc ~ log_live_dens + (1 + log_live_dens | year) + (1 | plot), data = plot_umca)
m2c <- lmer(log_tot_slc ~ log_live_dens + (1 | year) + (1 + log_live_dens | plot), data = plot_umca)

#models <- ls(pattern = "^m") # make vector of model object names; commented out b/c noticed that the row names are the model names of the anova
m_anova <- anova(m0, m1, m2, m0a, m0b, m1a, m1b, m2a, m2b, m2c)
m_anova$models <- row.names(m_anova)
#cbind(m_anova, models)
m_anova <- arrange(m_anova, AIC) # arrange by increasing AIC, so best model first
m_anova
anova(m2b, m1b, m2c) # Top 3 models

#' The random slope, random intercept model with the density:year interaction random effect and the plot random effect is the best, but only slightly better than the random intercept model with the year & plot intercept random effects, and significantly more better (based on change in AIC > 5) than the random slope, random intercept model with the density:plot interaction and year random effect.
anova(m1b, m2b) # Top 2 models, marginally "significant" difference
plot(m2b, main = "Model m2b")
plot(m1b, main = "Model m1b")
plot(m2c, main = "Model m2c")

#' Look at details of the 3 best models ####
m2b
m1b
m2c

#' View image of the Cholesky factor matrix "L" ####
#+ image Cholesky factor matrix “L" --------------------------------
image(getME(m2b, name = "L"))
image(getME(m1b, name = "L"))

image((getME(m1b, name = "Z")))
methods(class = "merMod")

#' Print Model Coefficients
coef(m2b)
coef(m1b)
coef(m2c)
#' Plot model residuals vs fitted values
plot((m2b))
plot((m1b))
plot((m2c))
#' Compute confidence intervals
confint(m2b)
confint(m1b)
confint(m2c)
#' Extract Residual Error
sigma(m0)
sigma(m1)
sigma(m2)

#' Plot model profiles using `lattice` package ####
library(lattice)
xyplot(profile(m2b), conf = c(50, 90, 95, 99)/100)
densityplot(profile(m2b))
splom(profile(m2b), conf = c(50, 90, 95, 99)/100)

#' Compare m2b to m1b based on predicted values
par(mfrow=c(1,3))
plot(predict(m1b), predict(m2b))
abline(lm(predict(m2b) ~ predict(m1b)), col = "red")

plot(predict(m2c), predict(m1b))
abline(lm(predict(m1b) ~ predict(m2c)), col = "red")

plot(predict(m2c), predict(m2b))
abline(lm(predict(m2b) ~ predict(m2c)), col = "red")

summary(lm(predict(m2b) ~ predict(m1b)))
summary(lm(predict(m1b) ~ predict(m2c)))
summary(lm(predict(m2b) ~ predict(m2c)))
