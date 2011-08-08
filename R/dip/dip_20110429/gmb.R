setwd("/home/lluc/Documents/R/dip/gamba")

source("/home/lluc/Documents/R/functions.R")

data<-read.csv("data.csv")
cond<-read.csv("cond.csv")
Fecha<-as.Date(Fecha)
cond.date<-as.Date(cond$Fecha.LF)

Tmin<-1
Tmax<-1
Tmed<-1
#RD<-1
DD<-1

for(i in 2:length(cond.date))
{
DD[i-1]<-sum(Dry.Days[Fecha<cond.date[i] & Fecha>cond.date[i-1] & is.na(Fecha)!=T])
#RD[i]<-sum(Lluvia[Fecha<cond.date[i] & Fecha>cond.date[i-1] & is.na(Fecha)!=T])
Tmin[i-1]<-mean(T.min[Fecha<cond.date[i] & Fecha>cond.date[i-1] & is.na(Fecha)!=T])
Tmax[i-1]<-mean(T.max[Fecha<cond.date[i] & Fecha>cond.date[i-1] & is.na(Fecha)!=T])
Tmed[i-1]<-mean(T.med[Fecha<cond.date[i] & Fecha>cond.date[i-1] & is.na(Fecha)!=T])
}

write.csv(data.frame(Tmin, Tmax, Tmed), "gmb.csv")

plot(cond.date[2:length(cond.date)], Tmed[2:length(cond.date)], type="n", ylim=c(20, max(Tmax[is.na(Tmax)!=T])))
lines(cond.date[2:length(cond.date)], Tmin[2:length(cond.date)], col="black")
lines(cond.date[2:length(cond.date)], Tmax[2:length(cond.date)], col="grey")
lines(cond.date[2:length(cond.date)], Tmed[2:length(cond.date)], col="darkgrey")
