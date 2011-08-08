alldata.ylabs<-alldata[1,T]
alldata.titles<-alldata[1,T]

alldata.ylabs[38:49]<-c(rep("% dry weight, CI (95%)", 3), rep("w/w, CI (95%)", 3), rep("% dry weight, CI (95%)", 3), rep("ng g-1 dry weight, CI (95%)", 3))
alldata.titles[38:49]<-c("Litter C content", "Litter N content", "Litter P content", "Litter C:N ratio", "Litter C:P ratio", "Litter N:P ratio", "Litter K content", "Litter Ca content", "Litter Mg content", "Litter Fe content", "Litter Mn content", "Litter Zn content")


pdf("output/controls_h1.pdf", width=18, height=10)

par(mfrow=c(2,6), tck=0.02)
for (i in c(38:49))
{
var<-alldata[,i]
cond1<-alldata$Harvest=="I"
sep<-alldata$Litter
name<-"DOC H1"
ylab="µg C g-1 d.w."

m<-tapply(var, list(sep, alldata$Harvest), mean)+tapply(var, list(sep, alldata$Harvest), function(x) CI(x, 0.04))
ylim<-max(m[,1:3])

      mean<-tapply(var[cond1], sep[cond1], mean)
      sds<-tapply(var[cond1], sep[cond1], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), cex.axis=1.5, cex.main=2, cex.lab=1.7, plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=alldata.titles[i], ylab=paste(ylab, "95% CI"), col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)
}
dev.off()

pdf("output/controls_h2.pdf", width=18, height=10)

par(mfrow=c(2,6), tck=0.02)
for (i in c(38:49))
{
var<-alldata[,i]
cond1<-alldata$Harvest=="II"
sep<-alldata$Litter
ylab="µg C g-1 d.w."

m<-tapply(var, list(sep, alldata$Harvest), mean)+tapply(var, list(sep, alldata$Harvest), function(x) CI(x, 0.04))
ylim<-max(m[,1:3])

      mean<-tapply(var[cond1], sep[cond1], mean)
      sds<-tapply(var[cond1], sep[cond1], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), cex.axis=1.5, cex.main=2, cex.lab=1.7, plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=alldata.titles[i], ylab=paste(ylab, "95% CI"), col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)
}
dev.off()

