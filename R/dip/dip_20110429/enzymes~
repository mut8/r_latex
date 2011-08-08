pdf("doc_barplots.pdf", width=12, height=12)

par(mfrow=c(4,4), tck=0.02)

var<-alldata$Cellulase
cond1<-alldata$Harvest=="I"
sep<-alldata$Litter
cond2<-alldata$Harvest=="I"
name<-"Cellulase H1"
ylab="nn, CI (95%)"

m<-tapply(var, list(sep, alldata$Harvest), mean)+tapply(var, list(sep, alldata$Harvest), function(x) CI(x, 0.04))
ylim<-max(m[,1:3])


      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab=paste(ylab, "95% CI"), col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Cellulase
cond1<-alldata$Harvest=="II"
sep<-alldata$Litter
cond2<-alldata$Harvest=="II"
name<-"Cellulase H2"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Cellulase
cond1<-alldata$Harvest=="III"
sep<-alldata$Litter
cond2<-alldata$Harvest=="III"
name<-"Cellulase H3"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Cellulase
cond1<-alldata$Harvest=="IV"
sep<-alldata$Litter
cond2<-alldata$Harvest=="IV"
name<-"Cellulase H4"

ylim<-max(m[,4])

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Chitinase
cond1<-alldata$Harvest=="I"
sep<-alldata$Litter
cond2<-alldata$Harvest=="I"
name<-"Chitinase H1"
ylab="nn, CI (95%)"

m<-tapply(var, list(sep, alldata$Harvest), mean)+tapply(var, list(sep, alldata$Harvest), function(x) CI(x, 0.04))
ylim<-max(m[,1:3])


      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab=paste(ylab, "95% CI"), col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Chitinase
cond1<-alldata$Harvest=="II"
sep<-alldata$Litter
cond2<-alldata$Harvest=="II"
name<-"Chitinase H2"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Chitinase
cond1<-alldata$Harvest=="III"
sep<-alldata$Litter
cond2<-alldata$Harvest=="III"
name<-"Chitinase H3"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Chitinase
cond1<-alldata$Harvest=="IV"
sep<-alldata$Litter
cond2<-alldata$Harvest=="IV"
name<-"Chitinase H4"

ylim<-max(m[,4])

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Phenoloxidase
cond1<-alldata$Harvest=="I"
sep<-alldata$Litter
cond2<-alldata$Harvest=="I"
name<-"Phenoloxidase H1"
ylab="nn, CI (95%)"

m<-tapply(var, list(sep, alldata$Harvest), mean)+tapply(var, list(sep, alldata$Harvest), function(x) CI(x, 0.04))
ylim<-max(m[,1:3])


      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab=paste(ylab, "95% CI"), col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Phenoloxidase
cond1<-alldata$Harvest=="II"
sep<-alldata$Litter
cond2<-alldata$Harvest=="II"
name<-"Phenoloxidase H2"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Phenoloxidase
cond1<-alldata$Harvest=="III"
sep<-alldata$Litter
cond2<-alldata$Harvest=="III"
name<-"Phenoloxidase H3"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Phenoloxidase
cond1<-alldata$Harvest=="IV"
sep<-alldata$Litter
cond2<-alldata$Harvest=="IV"
name<-"Phenoloxidase H4"

ylim<-max(m[,4])

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Peroxydase
cond1<-alldata$Harvest=="I"
sep<-alldata$Litter
cond2<-alldata$Harvest=="I"
name<-"Peroxydase H1"
ylab="nn, CI (95%)"

m<-tapply(var, list(sep, alldata$Harvest), mean)+tapply(var, list(sep, alldata$Harvest), function(x) CI(x, 0.04))
ylim<-max(m[,1:3])


      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab=paste(ylab, "95% CI"), col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Peroxydase
cond1<-alldata$Harvest=="II"
sep<-alldata$Litter
cond2<-alldata$Harvest=="II"
name<-"Peroxydase H2"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Peroxydase
cond1<-alldata$Harvest=="III"
sep<-alldata$Litter
cond2<-alldata$Harvest=="III"
name<-"Peroxydase H3"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Peroxydase
cond1<-alldata$Harvest=="IV"
sep<-alldata$Litter
cond2<-alldata$Harvest=="IV"
name<-"Peroxydase H4"

ylim<-max(m[,4])

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

var<-alldata$Peroxydase + alldata$Phenoloxidase
cond1<-alldata$Harvest=="I" 
sep<-alldata$Litter
cond2<-alldata$Harvest=="I"
name<-"Peroxydase + Phenoloxidase H1"
ylab="nn, CI (95%)"
zz

m<-tapply(var, list(sep, alldata$Harvest), mean)+tapply(var, list(sep, alldata$Harvest), function(x) CI(x, 0.04))
ylim<-max(m[,1:3])


      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab=paste(ylab, "95% CI"), col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)


cond1<-alldata$Harvest=="II"
sep<-alldata$Litter
cond2<-alldata$Harvest=="II"
name<-"Peroxydase + Phenoloxidase H2"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05)  

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

cond1<-alldata$Harvest=="III"
sep<-alldata$Litter
cond2<-alldata$Harvest=="III"
name<-"Peroxydase + Phenoloxidase H3"

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)

cond1<-alldata$Harvest=="IV"
sep<-alldata$Litter
cond2<-alldata$Harvest=="IV"
name<-"Peroxydase + Phenoloxidase H4"

ylim<-max(m[,4])

      mean<-tapply(var[cond1], sep[cond2], mean)
      sds<-tapply(var[cond1], sep[cond2], sd)
      ci2<-CI2(sds, 5, 0.05) 

      height<-mean
      error<-ci2

      barplot2(height, ylim=c(0, ylim*1.2), plot.ci=TRUE, ci.u=height+error,ci.l=height-error, names.arg=typlev, 
      main=name, ylab="", col=colscale, tck=0.01)
      aov<-aov(lm(var[cond1] ~ sep[cond2]))
      hsd<-HSD.test(aov, "sep[cond2]")
      gr<-hsd$M[order(hsd$trt)]
      text(((0:3)+.6)*1.2, ylim*1.15, gr)




dev.off()