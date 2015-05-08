0umca.alive = read.csv("D:\\SOCO\\SOCO_space_time\\SatScan_SOCO\\UMCA\\UMCA_alive_2004-11.csv")
attach(umca.alive)
head(umca.alive)
# plot.id tag  spp dbh year stem.status slc log10_slc
#1   ANN01 807 UMCA  21 2004       alive   0 0.0000000
#2   ANN01 807 UMCA  21 2005       alive   2 0.4771213
#3   ANN01 807 UMCA  21 2006       alive  12 1.1139434
#4   ANN01 807 UMCA  21 2007       alive  11 1.0791812
#5   ANN01 807 UMCA  21 2008       alive  27 1.4471580
#6   ANN01 807 UMCA  21 2009       alive  40 1.6127839
#umca = read.csv("D:\\SOCO\\Plot_descriptors\\UMCA_stems_2004-11.csv")
#umca = read.csv("C:\\Users\\Sarah\\Desktop\\standarderrorbars\\UMCA_stems_2004-11.csv")


############## Tagged UMCA:
umca.stems = unique(tag)
#length(umca.stems)
#2748
plot.id.uniq = unique(plot.id)
#length(plot.id.uniq)
#180

################ SLC:
annual.slc = tapply(slc, list(plot.id, year), sum, na.rm=TRUE)
#umca.annual.slc;  Good, these numbers match my "umca_cumsum_slc" file in the satscan folder.
#ANN01     505  445 1012 1069 1328 1604 1260 2700
#ANN02     184  133  472  389  480  657  566  394
#ANN03     950 1274 1273 1058  757  946 1468 4684
#ANN04      66   89  138  542  576  876  542  272
#ANN05     685  455 1128  622  472  506  269  787
#ANN06     101 1167    2  649  626  506  437  540
#ANN07       4  131   24   69   55   73   80  112
#is.matrix(annual.slc), #[1] TRUE
#annual.slc.export = write.csv(annual.slc, file="D:\\SOCO\\Plot_descriptors\\annual.slc.export.csv")
annual.slc.df = as.data.frame(annual.slc)

# Matrix algegra to get average slc/year (for line graph):
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

####################################################################################################################
#**************************************************** BOXPLOTS:

# Trying the "jitter plot" from library(ggplot2):
library(ggplot2)
#x1 = qplot(year, slc, geom="jitter", alpha=I(1/10), na.rm=TRUE)  #Remember, alpha changes teh transparency.
#  x2 + opts(panel.background = theme_rect(fill='ghostwhite'))
#ggsave(plot=x1, filename="D:\\SOCO\\SOCO_space_time\\Figures\\annual_slc_jitter.jpeg", height=4, width=5)

# Regular boxplot:
x2 = ggplot(umca, aes(factor(year), slc))
x2b=  x2 + geom_boxplot(na.rm=TRUE) +
  stat_summary(fun.y=mean, na.rm=TRUE, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, na.rm=TRUE,geom="point",shape=8) +
  ylab("x2")
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
