### Exploring UMCA stem infection through time

setwd("~/Dropbox/DISSERTATION/superspreaders/analysis")

library(foreign)
plot_dat1<-read.dbf("data/UMCA_slc_dbh_04-12.dbf")
str(plot_dat1)
summary(plot_dat1)
head(plot_dat1)

plot_dat1$avg_slc_04_12<-with(plot_dat1, rowMeans(cbind(avgslc04,avgslc05,avgslc06,avgslc07,avgslc08,avgslc09,avgslc10,avgslc11,avgslc12)))

years<-c("2004","2005","2006","2007","2008","2009","2010","2011","2012")
with(plot_dat1, boxplot(avgslc04,avgslc05,avgslc06,avgslc07,avgslc08,
        avgslc09,avgslc10,avgslc11,avgslc12,names=years,xlab="Year",ylab="Mean SLC"))



library(reshape)
plot_dat1<-rename(plot_dat1,c(PLOT_ID="plotid"))
library(dplyr)
plot_dat1<-left_join(plot_dat1,plot_202_rain,by="plotid")

avg_slc_year<-with(plot_dat1, colMeans(cbind(avgslc04,avgslc05,avgslc06,avgslc07,avgslc08,
                         avgslc09,avgslc10,avgslc11,avgslc12)))


plot(years,avg_slc_year,xlab="Year",ylab="Mean SLC")

avg_rain_year<-with(plot_202_rain, colMeans(cbind(amou2004,amou2005,amou2006,amou2007,amou2008,
                                                  amou2009,amou2010,amou2011,amou2012)))
rainy_days<-with(plot_202_rain, colMeans(cbind(freq2004,freq2005,freq2006,freq2007,freq2008,
                                               freq2009,freq2010,freq2011,freq2012)))

df_rain_slc<-as.data.frame(rbind(avg_rain_year,avg_slc_year))

matplot(cbind(avg_rain_year,avg_slc_year),type="b")
matplot(cbind(rainy_days,avg_slc_year),type="b")



library(plotrix)
twoord.stackplot(lx=years,rx=years,ldata=avg_slc_year,rdata=avg_rain_year,lcol="red",rcol="blue",
                 ltype="b",rtype="b",rylab="Mean SLC",lylab="Mean Rainfall",xlab="Year")

library(ggplot2)
qplot(years,avg_slc_year,geom=c("point","line"))

ggplot() +
  geom_line(data=bay_dat1,aes(x=years,y=colMeans(cbind(amou2004,amou2005,amou2006,amou2007,amou2008,
                                                        amou2009,amou2010,amou2011,amou2012)),
                              colour="blue")) +
  geom_line(data=bay_dat1,aes(x=years,y=colMeans(cbind(cbind(avgslc04,avgslc05,avgslc06,avgslc07,avgslc08,
                                                              avgslc09,avgslc10,avgslc11,avgslc12))),
                              colour="red"))
