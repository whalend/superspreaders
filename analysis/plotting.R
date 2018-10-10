#' A lot of the plotting and code in this was written by Sarah Haas. I have made some additions and modifications while exploring ways to quantify what may constitute a superspreader individual.

library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

umca_stem <- read_csv("SH_data/UMCA_stems_live_2004-11.csv")

head(umca_stem)
summary(umca_stem %>% group_by(year, tag))
summary(umca_stem %>% group_by(year))

#' Get specific quantile values for SLC for each year ####
#' If the same stems are consistently in the upper quantile then they may be superspreaders because they support more infections compared to the population. 
#' The `plyr` method with sample code from <http://stackoverflow.com/questions/5473537/how-to-calculate-95th-percentile-of-values-with-grouping-variable-in-r-or-excel>

#Random seed
set.seed(42)
#Sample data
dat <- data.frame(Watershed = sample(letters[1:2], 100, TRUE), WQ = rnorm(100))
#plyr call
ddply(dat, "Watershed", summarise, WQ95 = quantile(WQ, .95))

dat %>% 
      group_by(Watershed) %>% 
      summarise(WQ95 = quantile(WQ, 0.95))

ddply(umca_stem, "year", summarise, slc95 = quantile(slc, .95))

# The `dplyr` method
umca_stem %>% group_by(year) %>% 
      summarise(slc90 = quantile(slc, 0.90), slc95 = quantile(slc, 0.95),
                slc85 = quantile(slc, 0.85), slc80 = quantile(slc, 0.80),
                slc99 = quantile(slc, 0.99))

data90 <- umca_stem %>% group_by(year) %>% 
      filter(year == 2004 & slc >= 82 | year == 2005 & slc >= 63 | 
                   year == 2006 & slc >= 112 | year == 2007 & slc >= 72 | 
                   year == 2008 & slc >= 62 | year == 2009 & slc >= 65 | 
                   year == 2010 & slc >= 104 | year == 2011 & slc >= 166)

data95 <- umca_stem %>% group_by(year) %>% 
      filter(year == 2004 & slc >= 114 | year == 2005 & slc >= 88 | 
                   year == 2006 & slc >= 130 | year == 2007 & slc >= 86 | 
                   year == 2008 & slc >= 76 | year == 2009 & slc >= 83 | 
                   year == 2010 & slc >= 139 | year == 2011 & slc >= 190)
data95 <- data95 %>% arrange(tag) %>% select(plot.id, tag, year, slc)
d <- subset(data95, duplicated(data95$tag))
qplot(year, slc, data = d) + geom_line(aes(color = tag))

p1 <- ggplot(data95, aes(x = year, y = slc, tag = tag))

p1 + geom_point(aes(color = tag), position = "jitter", alpha = 0.3) +
      #geom_line(aes(color = tag)) +
      theme(legend.position = "none") +
      geom_text(aes(label = tag))

p2 <- ggplot(data = umca_stem %>% group_by(tag, year, plot.id) %>% 
                   filter(slc >= 150), 
             aes(x = year, y = slc))

p2 <- p2 + geom_point(aes(color = plot.id), position = "jitter") +
      theme(legend.position = "top") +
      geom_line(aes(color = tag))

p2 + geom_text(aes(label = tag, color = ))


#' "Seasonality" and symptoms ####
#' Note that this has to be done within a sample season, so it's a bad assessment of true seasonality of symptoms. Jennifer Davidson's work addressed this pretty well, so it would be a proper reference.

ggplot(data = stems %>% select(Date, slc, species, status, year) %>% filter(species == "UMCA" & status == "Alive") %>% group_by(Date)) + 
      geom_point(aes(x=Date,y=slc)) +
      facet_grid(year ~ .)


#' Sarah's code below ####
annual.slc = tapply(slc, list(plot.id, year), sum, na.rm=TRUE) 
#is.matrix(annual.slc), #[1] TRUE
#annual.slc.export = write.csv(annual.slc, file="D:\\SOCO\\Plot_descriptors\\annual.slc.export.csv")
annual.slc.df = as.data.frame(annual.slc)

#' Matrix algebra to get average slc/year (for line graph):
annual.slc.avg = as.integer(colMeans(annual.slc, na.rm=TRUE))
#annual.slc.avg
#     2004  2005  2006  2007  2008  2009  2010  2011 
#[1]  361    293   681   458   430   409   606   1060
# I'm thinking that I DON'T want to plot avgs b/c they are not as informative as my boxplot.
# Now, I want the 'sd' for each year:
annual.slc.sd = as.integer(sd(annual.slc, na.rm=TRUE))
#    2004 2005  2006 2007  2008  2009  2010  2011 
#[1]  477  394  907  447   389   407   659   1153
# I'm having a hell of a time getting 'stderr' values in R, so I did it myself in excel. It matched up for 
  # 2011, which is the only year that didn't contain "na" values.
annual.slc.stderr = as.vector(c(35.6, 29.4, 67.7, 33.4, 29.1, 30.4, 49.1, 86.0))
#stderr <- function(x) sqrt(var(x)/length(x))


# BOXPLOTS ####

# Trying the "jitter plot" from library(ggplot2):
library(ggplot2)
#x1 = qplot(year, slc, geom="jitter", alpha=I(1/10), na.rm=TRUE)  #Remember, alpha changes teh transparency.
#  x2 + opts(panel.background = theme_rect(fill='ghostwhite'))
#ggsave(plot=x1, filename="D:\\SOCO\\SOCO_space_time\\Figures\\annual_slc_jitter.jpeg", height=4, width=5)

# Regular boxplot:
p1 <- ggplot(data = umca_stem %>% group_by(tag, year), 
             aes(x = year, y = slc, tag = tag))

p1 + geom_point(position = "jitter") +
      geom_text(aes(label = tag))

p2 <- ggplot(data = umca_stem %>% group_by(tag, year, plot.id) %>% 
                   filter(slc >= 150), 
             aes(x = year, y = slc))

p2 <- p2 + geom_point(aes(color = plot.id), position = "jitter") +
      theme(legend.position = "top") +
      geom_line(aes(color = tag))
      
p2 + geom_text(aes(label = tag))


p1 + geom_boxplot(na.rm=TRUE) +
  stat_summary(fun.y=mean, na.rm=TRUE, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, na.rm=TRUE,geom="point",shape=8) +
  ylab("Symptomatic Leaf Count")
x2b
#ggsave(plot=x2b, filename="C:\\Users\\Sarah\\Desktop\\SOCO_Nov23\\annual_slc_boxplot.jpeg", height=4, width=5)

# Log_base10:
x3 = ggplot(umca, aes(factor(year), log10(slc + 1)))
x3b=  x3 + geom_boxplot(na.rm=TRUE) +
  ylab(expression(paste(Log[10], " of Inoculum Load", sep=""))) +
  xlab("") +
  stat_summary(fun.y=mean, na.rm=TRUE, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, na.rm=TRUE,geom="point",shape=8)
x3b
#ggsave(plot=x3b, filename="D:\\SOCO\\SOCO_space_time\\Figures\\annual_log10slc_boxplot.jpeg", height=4, width=5)

################################################ Avec Errorbar:
library(plyr)
cdata= ddply(umca, .(year), summarise,
             N = sum(!is.na(slc)),
             mean = mean(slc, na.rm=TRUE),
             sd = sd(slc, na.rm=TRUE),
             se = sd(slc, na.rm=TRUE)/sqrt(sum(!is.na(slc))))
#year    N   mean       sd        se
#1 2004 2428 26.18369 40.23642 0.8165729
#2 2005 2448 21.22631 29.82954 0.6028939
#3 2006 2505 48.14890 44.30052 0.8851256
#4 2007 2420 33.32603 26.94906 0.5478176
#5 2008 2563 29.93094 23.71083 0.4683520
#6 2009 2475 29.42101 26.58503 0.5343792
#7 2010 2550 42.31569 43.29772 0.8574226
#8 2011 2583 73.90825 59.46067 1.1699508

cdata.umca.alive= ddply(umca.alive, .(year), summarise,
             N = length(log10_slc),
             mean = mean(log10_slc),
             sd = sd(log10_slc),
             se = sd(log10_slc)/sqrt(length(log10_slc)))
cdata.umca.alive
#year    N      mean        sd          se
#1 2004 2428 0.9135364 0.7399910 0.015017655
#2 2005 2448 0.8764221 0.7118471 0.014387357
#3 2006 2505 1.3329451 0.7287718 0.014560882
#4 2007 2420 1.3379578 0.5027549 0.010219947
#5 2008 2563 1.3450055 0.3964459 0.007830864
#6 2009 2475 1.2636763 0.5207769 0.010468010
#7 2010 2550 1.4156478 0.4787359 0.009480383
#8 2011 2583 1.6758250 0.4913556 0.009667935

#cdatalog= ddply(umca, .(year), summarise,  # Identical to "cdata.umca.alive":
#             N = sum(!is.na(slc)),
#             mean = mean(log10(slc+1), na.rm=TRUE),
#             sd = sd(log10(slc+1), na.rm=TRUE),
#             se = sd(log10(slc+1), na.rm=TRUE)/sqrt(sum(!is.na(slc))))
#cdatalog
#year    N      mean        sd          se
#1 2004 2428 0.9135364 0.7399910 0.015017655
#2 2005 2448 0.8764221 0.7118471 0.014387357
#3 2006 2505 1.3329451 0.7287718 0.014560882
#4 2007 2420 1.3379578 0.5027549 0.010219947
#5 2008 2563 1.3450055 0.3964459 0.007830864
#6 2009 2475 1.2636763 0.5207769 0.010468010
#7 2010 2550 1.4156478 0.4787359 0.009480383
#8 2011 2583 1.6758250 0.4913556 0.009667935


x5c = ggplot(cdata, aes(factor(year), mean, group=1)) +
  geom_boxplot(data=umca, aes(factor(year), slc, group=year),na.rm=TRUE) +
  geom_line() + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, colour="red") +geom_point(size=3) +
  ylab("Figure x5c")
x5c #Nothing logged.

x5c1 = ggplot(cdata.umca.alive, aes(factor(year), mean, group=1)) + theme_gray(base_size = 9) +
  geom_boxplot(data=umca.alive, aes(factor(year), log10_slc, group=year)) +
  geom_line(colour="red", lwd=.3) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, colour="red",lwd=.3) +geom_point(size=1, colour="red") +
  ylab(expression(paste(Log[10], " of Inoculum Load", sep=""))) +
  xlab("") 
x5c1 
ggsave(plot=x5c1, filename="D:\\SOCO\\SOCO_space_time\\Figures\\annual_log10slc_boxplot.png", height=4, width=5)
#Put simply, standard error is an estimate of how close to the population mean your sample mean is likely to be, 
  #whereas standard deviation is the degree to which individuals within the sample differ from the sample mean. 
  #Standard error should decrease with larger sample sizes, as the estimate of the population mean improves. 
  #Standard deviation will be unaffected by sample size.

# Tried using standard deviation:
x5c1b = ggplot(cdata.umca.alive, aes(factor(year), mean, group=1)) + theme_gray(base_size = 9) +
  geom_boxplot(data=umca.alive, aes(factor(year), log10_slc, group=year)) +
  geom_line(colour="red", lwd=.3) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, colour="red",lwd=.3) +geom_point(size=2, colour="red") +
  ylab(expression(paste(Log[10], " of Inoculum Load", sep=""))) +
  xlab("") 
x5c1b 
ggsave(plot=x5c1b, filename="D:\\SOCO\\SOCO_space_time\\Figures\\annual_log10slc_boxplot_sd.png", height=4, width=5)

################################################ Errorbar Completed:



# What Anya/Whalen and I did; this is just the 'se' line:
g5c = ggplot(cdata, aes(factor(year), mean, group=1)) +
 geom_line() + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, colour="red") +geom_point(size=3)
g5c

g5d = ggplot(cdata, aes(factor(year), log10(mean), group=1)) +
  geom_line() + geom_errorbar(aes(ymin=log10(mean)-log10(se), ymax=log10(mean)+log10(se)), width=.1, colour="red") +geom_point(size=3)
g5d


#*****************************************************>
#*****************************************************>
#Unless there's an easier within ggplot approach, what I'd do is write a function to calculate 
# the upper and lower ends of the error bar. Here, I've got two that calculate 95% confidence intervals. 

errorUpper <- function(x){ 
  x.mean <- mean(x) 
  x.sd <- sd(x) 
  SEM <- x.sd / (sqrt(length(x)))  
  return(x.mean + (SEM*1.96)) 
} 

errorLower <- function(x){ 
  x.mean <- mean(x) 
  x.sd <- sd(x)  
  SEM <- x.sd / (sqrt(length(x))) 
  return(x.mean - (SEM*1.96)) 
} 

#Now, the ggplot code would look like this: 
#ggplot(mtcars, aes(x=gear, y=mpg, group=as.factor(cyl), 
#colour=as.factor(cyl))) + 
# stat_summary(fun.y=mean, geom="point") + 
# stat_summary(fun.y=mean, geom="line")+ 
# stat_summary(fun.ymax = errorUpper, fun.ymin = errorLower, geom = 
#"errorbar") 

x4 = ggplot(umca.alive, aes(factor(year), log10_slc)) + theme_gray(base_size = 9) 
x4=  x4 + geom_boxplot() +
  ylab(expression(paste(Log[10], " of Inoculum Load", sep=""))) +
  xlab("") +
  stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red", lwd=.4)  + 
  stat_summary(fun.ymax=errorUpper, fun.ymin=errorLower, geom="errorbar", width=0.2, colour="red", lwd=.3) +
  stat_summary(fun.y=mean, geom="point", colour="red",size=1.5)
  x4
ggsave(plot=x4, filename="D:\\SOCO\\SOCO_space_time\\Figures\\annual_log10slc_boxplot_95CI.png", height=4, width=5)


#A ribbon might look nicer:
#ggplot(mtcars, aes(x=gear, y=mpg, group=as.factor(cyl), 
#colour=as.factor(cyl))) + 
# stat_summary(fun.ymax = errorUpper, fun.ymin = errorLower, geom = 
#"ribbon", alpha = 0.6)+ 
# stat_summary(fun.y=mean, geom="point") + 
# stat_summary(fun.y=mean, geom="line") 
