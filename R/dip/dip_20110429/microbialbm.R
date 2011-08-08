pdf("microbialbiomass.pdf", width=12, height=12)

par(mfcol=c(3,2))

var<-alldata$C_mic
cond<-is.na(var)!=T
write.csv(timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], xlab="incubation (days)", ylab="microbial C (mg C g-1 d.w., SE)",
pch=pch, col=1, cex=2, errcol=1, lwd=.5, type="o", bg="blue", endsig=T, letters=c(0.05,0.00)), "export/timeseriesmicbm1.csv")


var<-alldata$N_.mic
cond<-is.na(var)!=T
write.csv(timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], xlab="incubation (days)", ylab="microbial N (mg N g-1 d.w., SE)",
pch=pch, col=1, cex=2, errcol=1, lwd=.5, type="o", bg="blue", endsig=T, letters=c(0.05,0.00)), "export/timeseriesmicbm2.csv")

var<-alldata$P_mic
cond<-is.na(var)!=T
write.csv(timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], xlab="incubation (days)", ylab="microbial P (mg P g-1 d.w., SE)",
pch=pch, col=1, cex=2, errcol=1, lwd=.5, type="o", bg="blue", endsig=T, letters=c(0.05,0.00)), "export/timeseriesmicbm3.csv")

var<-alldata$C.N_mic
cond<-is.na(var)!=T
write.csv(timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], xlab="incubation (days)", ylab="microbial C:N (mass ratio, SE)",
pch=pch, col=1, cex=2, errcol=1, lwd=.5, type="o", bg="blue", endsig=T, letters=c(0.05,0.00)), "export/timeseriesmicbm4.csv")

var<-alldata$C.Pmic
cond<-is.na(var)!=T
write.csv(timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], xlab="incubation (days)", ylab="microbial C:P ratio (mass ratio, SE)",
pch=pch, col=1, cex=2, errcol=1, lwd=.5, type="o", bg="blue", endsig=T, letters=c(0.05,0.00)), "export/timeseriesmicbm5.csv")

var<-alldata$N.Pmic
cond<-is.na(var)!=T
write.csv(timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], xlab="incubation (days)", ylab="microbial N:P ratio (mass ratio, SE)",
pch=pch, col=1, cex=2, errcol=1, lwd=.5, type="o", bg="blue", endsig=T, letters=c(0.05,0.00)), "export/timeseriesmicbm6.csv")

dev.off()