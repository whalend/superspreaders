
## Repeated measures analysis ####

#### Read in data for all stems that has been processed in "explore_data.R"
stems <- read.csv("analysis/data/stems_growth.csv")
str(stems)
summary(stems)


library(plyr); library(dplyr); library(tidyr); library(ggplot2)
names(stems)
stems <- stems %>%
      select(-X)# Remove the row names variable
stems$date <- as.Date(stems$date)
stems$year <- as.factor(stems$year)

### Create bay laurel only data frame ####
umca <- as.tbl(stems) %>%
      select(plot, cluster, tag, species, status, slc, dbh, delta_dbh, dbh2_1_ratio, year, date) %>%
      filter(species == "UMCA", status == "Alive" | status == "Dead")
summary(umca)
umca <- droplevels(umca)

# Develop plot level summaries data frame and join to stems data frame
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

#### Plot Level Repeated Measures Analysis Using Multi-level Mixed-Effects Regression ####

## Linear Mixed Effects Method ####
#### A couple of lines from Sarah's code to work from
# anova(lmer.null, lmer1)
# Models:
# lmer.null: log(cum.slc + 1) ~ 1 + (1 | year)
# lmer1: log(cum.slc + 1) ~ log(umca.ct + 1) + (1 | year)
#          Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# lmer.null  3 4378.9 4394.0 -2186.4   4372.9                             
# lmer1      4 3861.9 3882.1 -1927.0   3853.9 518.97      1  < 2.2e-16 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# lmer2 <- lmer(log(cum.slc + 1) ~ log(umca.ct + 1) + (1 + log(umca.ct+1)|year), na.action=na.omit)
# anova(lmer1, lmer2)
# Data: 
#  Models:
#  lmer1: log(cum.slc + 1) ~ log(umca.ct + 1) + (1 | year)
#  lmer2: log(cum.slc + 1) ~ log(umca.ct + 1) + (1 + log(umca.ct + 1) | 
#                                               lmer2:     year)
#       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# lmer1  4 3861.9 3882.1 -1927.0   3853.9                         
# lmer2  6 3865.2 3895.4 -1926.6   3853.2 0.7407      2     0.6905

# Starting with plot level data b/c simpler w/fewer nested groups
# Assumptions of linear mixed effects model:
#     1. Gaussian (normal) distribution
#     2. Errors ~ N(u, var)
library(lme4)
#### Simplest models, with only one random effect
m0 <- lmer(tot_slc ~ 1 + (1 | year), data = plot_umca)
m0a <- lmer(tot_slc ~ 1 + (1 | year) + (1 | plot), data = plot_umca)
m1 <- lmer(tot_slc ~ live_dens + (1 | year), data = plot_umca)
m2 <- lmer(tot_slc ~ live_dens + (live_dens | year), data = plot_umca)
anova(m0, m1, m2)

m0 <- lmer(log1p(tot_slc) ~ 1 + (1 | year), data = plot_umca)
m1 <- lmer(log1p(tot_slc) ~ log1p(live_dens) + (1 | year), data = plot_umca)
m2 <- lmer(log1p(tot_slc) ~ log1p(live_dens) + (log1p(live_dens) | year), data = plot_umca)
anova(m0, m1, m2) # Qualitatively similar results to SH's output
anova(m1, m2)

#### Print Model Coefficients
coef(m0) # intercept only (null) model
coef(m1) # varying intercept
coef(m2) # vayring intercept, varying slope
#### Plot model residuals
plot(residuals(m0))
plot(residuals(m1))
plot(residuals(m2))
#### Compute confidence intervals
confint(m0)
confint(m1)
confint(m2)
#### Extract Residual Error
sigma(m0)
sigma(m1)
sigma(m2)

#### Plot model profiles using `lattice` package
library(lattice)
xyplot(profile(m1), conf = c(50, 90, 95, 99)/100)
densityplot(profile(m1))
splom(profile(m1), conf = c(50, 90, 95, 99)/100)

#### Compare m1 to m2 based on predicted values
par(mfrow=c(1,3))
plot(predict(m1), predict(m2))
abline(lm(predict(m1) ~ predict(m2)), col = "red")

plot(predict(m0), predict(m1))
abline(lm(predict(m0) ~ predict(m1)), col - "red")

plot(predict(m0), predict(m2))
abline(lm(predict(m0) ~ predict(m2)), col - "red")

summary(lm(predict(m1) ~ predict(m2)))
summary(lm(predict(m0) ~ predict(m1)))
summary(lm(predict(m0) ~ predict(m2)))





