# Repeated measures analysis ####
setwd("~/GitHub/superspreaders/analysis")
stems <- read.csv("data/stems.csv")
str(stems)
summary(stems)

### Change the assigned missing data value to NA across entire data frame
stems[stems == -9999] <- NA

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

### Create bay laurel only data frame ####
umca <- as.tbl(stems) %>%
      select(plot, cluster, tag, species, status, slc, dbh, year) %>%
      filter(species == "UMCA")

# Develop plot level summaries data frame and join to stems data frame
umca <- left_join(umca, plot_umca <- umca %>%
      select(plot, cluster, tag, status, slc, year) %>% 
      group_by(plot, year) %>%
      filter(slc >= 0, status == "Alive") %>%
      summarise(avg_slc = mean(slc), tot_slc = sum(slc), 
                live_dens = length(tag)), 
      by = c("plot", "year"))
summary(umca)


###Linear Mixed Effects Method ####
# Starting with plot level data b/c simpler w/fewer nested groups
library(lme4)
m0 <- lmer(tot_slc ~ 1 + (1 | year), data = plot_umca)
m1 <- lmer(tot_slc ~ live_dens + (1 | year), data = plot_umca)
m2 <- lmer(tot_slc ~ live_dens + (live_dens | year), data = plot_umca)
anova(m0, m1, m2)

m0 <- lmer(log1p(tot_slc) ~ 1 + (1 | year), data = plot_umca)
m1 <- lmer(log1p(tot_slc) ~ log1p(live_dens) + (1 | year), data = plot_umca)
m2 <- lmer(log1p(tot_slc) ~ log1p(live_dens) + (live_dens | year), data = plot_umca)
anova(m0, m1, m2)

par(mfrow=c(1,2))
plot(predict(m1), predict(m2))
abline(lm(predict(m1) ~ predict(m2)), col = "red")

summary(lm(predict(m1) ~ predict(m2)))

anova(lmer.null, lmer1)
#Models:
#  lmer.null: log(cum.slc + 1) ~ 1 + (1 | year)
#lmer1: log(cum.slc + 1) ~ log(umca.ct + 1) + (1 | year)
#          Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#lmer.null  3 4378.9 4394.0 -2186.4   4372.9                             
#lmer1      4 3861.9 3882.1 -1927.0   3853.9 518.97      1  < 2.2e-16 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

