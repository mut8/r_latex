pdf("output/timeseries_inbalance.pdf")

par(mfcol=c(2,2))

var<-CN_inbal
cond<-is.na(var)!=T
ylab<-"resource C:N / microbial C:N, +/- 1SE)"
write.csv(timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], xlab="incubation time (days)", ylab=ylab, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05)), "export/inbal1.csv")
legend("topright", pch=pch, typlev, cex=1.2)

var<-CP_inbal
cond<-is.na(var)!=T
ylab<-"resource C:P / microbial C:P, +/- 1SE)"
write.csv(timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], xlab="incubation time (days)", ylab=ylab, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05)), "export/inbal2.csv")

var<-NP_inbal
cond<-is.na(var)!=T
ylab<-"resource N:P / microbial N:P, +/- 1SE)"
write.csv(timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], xlab="incubation time (days)", ylab=ylab, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05)), "export/inbal3.csv")

dev.off()


CrcTICdif<-rowSums(difinit.rcTIC[,peaks$origin=="C"])
LrcTICdif<-rowSums(difinit.rcTIC[,peaks$origin=="L"])
timeseries(LrcTICdif, days, type, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05))
timeseries(CrcTICdif, days, type, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05))


timeseries(CrcTICdif, days, type, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05))


exclude<-c(rep(T, 29), F, rep(T, 49), F)
alldata$Litter[alldata$days==97&exclude]
cor.test(-LrcTICdif[harvest==2],CN_inbal[alldata$days==97&exclude])
plot(LrcTICdif[harvest==2],CN_inbal[alldata$days==97&exclude]) 

cor.test(-CrcTICdif[harvest==2],CN_inbal[alldata$days==97&exclude])
plot(LrcTICdif[harvest==2],CN_inbal[alldata$days==97&exclude]) 


cor.test(-LrcTICdif[harvest==6],alldata$C.N_mic[alldata$days==181])

pdf("output/substrateuse_cninbalance.pdf")
cor.test(-LrcTICdif[harvest==6],CN_inbal[alldata$days==181])
plot(CN_inbal[alldata$days==181],-LrcTICdif[harvest==6], ylab="Lignin decomposition (% cTIC)", xlab="C:N inbalance (C:N litter / C:N micr.)", pch=16)  
abline(lm(-LrcTICdif[harvest==6] ~ CN_inbal[alldata$days==181]))

cor.test(-CrcTICdif[harvest==6],CN_inbal[alldata$days==181])
plot(CN_inbal[alldata$days==181], -CrcTICdif[harvest==6], ylab="Carbohydrate decomposition (% cTIC)", xlab="C:N inbalance", pch=16) 
abline(lm(-CrcTICdif[harvest==6] ~ CN_inbal[alldata$days==181]))

cor.test(LrcTICdif[harvest==6]/CrcTICdif[harvest==6],CN_inbal[alldata$days==181])

dev.off()

cor.test(-LrcTICdif[harvest==15],CN_inbal[alldata$days==375&exclude])
plot(LrcTICdif[harvest==15],CN_inbal[alldata$days==375&exclude]) 


cor.test(-CrcTICdif[harvest==15],CN_inbal[alldata$days==375&exclude])
plot(CrcTICdif[harvest==15],CN_inbal[alldata$days==375&exclude]) 


alldata$Litter[alldata$days==181]


  pdf("output/timeseries_lci.pdf", width=12, height=6)

  par(mfrow=c(1,2), oma=c(2,.5,0,.5))

  par(mar=c(2.1,5.1,2.1,0.1))

write.csv(timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05)), "export/timeserieslci01.csv")
  timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05))

  mtext("Lignin/(Lignin+Carbohydrates) \n peak area (%TIC, SE)", side=2, padj=-1.5)
  axis(1, tck=0.01)
  axis(4, tck=0.01, labels=F)
  mtext("litter incubation (days)", side=1, padj=3)


  par(mar=c(2.1,0.1,2.1,5.1))
write.csv(timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, masslossSE=masslossSE, massloss=massloss, normalize=1,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, col=1, type="o", yaxt="n", xaxt="n"), "export/timeserieslci02.csv")
  timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, masslossSE=masslossSE, massloss=massloss, normalize=1,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, col=1, type="o", yaxt="n", xaxt="n")
  abline(h=100, col="grey", lty=2 )
# curve(100/(1-x), add=T, col="grey", lty=3)
  #timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, nam="Lignin/(Lignin+Carbohydrates)", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
  #abline(h=100, col="grey")s
  axis(4, tck=0.01)
  axis(1, tck=0.01)
  axis(2, tck=0.01, labels=F)

  mtext("accumulated respiration (g CO2-C g-1 litter-C, SE)", side=1, padj=3)
  mtext("Lignin/(Lignin+Carbohydrates) \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

#   par(mar=c(2.1,5.1,0.6,1.1))
#   timeseries(orig_cTIC$L/orig_cTIC$N, days, type,  errcol=errcol, pch=pch, lwd=lwd, endsig=T, add=F, cex=cex,  type="o", letters=c(0.05,0.05))
#   mtext("Lignin/Nitrogen compounds \n peak area (%TIC, SE)", side=2, padj=-1.5)
#   axis(4, tck=0.01, labels=F)
#   mtext("litter incubation (days)", side=1, padj=3)

#   par(mar=c(2.1,1.1,0.6,5.1))
#   timeseries(orig_cTIC$L/orig_cTIC$N, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n")
#   abline(h=100, col="grey", endsig=T, letters=c(0.05,0.05))
#   #timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#   #abline(h=100, col="grey")
#   axis(2, tck=0.01, labels=F)
#   axis(4, tck=0.01)
#   mtext("Lignin/Nitrogen compounds, \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)
#   mtext("accumulated respiration (g CO2-C g-1 litter-C, SE)", side=1, padj=3)
  legend("topleft", pch=pch, typlev, cex=1.2)

dev.off()

  pdf("output/timeseries_lphci.pdf", width=12, height=6)

  par(mfrow=c(1,2), oma=c(2,.5,0,.5))

  par(mar=c(2.1,5.1,2.1,0.1))

write.csv(timeseries((orig_cTIC$L+orig_cTIC$Ph)/(orig_cTIC$L+orig_cTIC$C+orig_cTIC$Ph), days, type, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05)), "export/timeseries03.csv")
  timeseries((orig_cTIC$L+orig_cTIC$Ph)/(orig_cTIC$L+orig_cTIC$C+orig_cTIC$Ph), days, type, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05))
  mtext("(Lignin+Phenoles)/(Lignin+Phenoles+Carbohydrates) \n peak area (%TIC, SE)", side=2, padj=-1.5)
  axis(1, tck=0.01)
  axis(4, tck=0.01, labels=F)
  mtext("litter incubation (days)", side=1, padj=3)


  par(mar=c(2.1,0.1,2.1,5.1))
write.csv(timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, errcol=errcol, pch=pch, lwd=lwd, endsig=T, cex=cex, type="o", xaxt="n", letters=c(0.05,0.05)), "export/timeserieslci01.csv")
  timeseries((orig_cTIC$L+orig_cTIC$Ph)/(orig_cTIC$L+orig_cTIC$C+orig_cTIC$Ph), days, type, masslossSE=masslossSE, massloss=massloss, normalize=1,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, col=1, type="o", yaxt="n", xaxt="n")
  abline(h=100, col="grey", lty=2 )
# curve(100/(1-x), add=T, col="grey", lty=3)
  #timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, nam="Lignin/(Lignin+Carbohydrates)", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
  #abline(h=100, col="grey")s
  axis(4, tck=0.01)
  axis(1, tck=0.01)
  axis(2, tck=0.01, labels=F)

  mtext("accumulated respiration (g CO2-C g-1 litter-C, SE)", side=1, padj=3)
  mtext("(Lignin+Phenoles)/(Lignin+Phenoles+Carbohydrates) \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

#   par(mar=c(2.1,5.1,0.6,1.1))
#   timeseries(orig_cTIC$L/orig_cTIC$N, days, type,  errcol=errcol, pch=pch, lwd=lwd, endsig=T, add=F, cex=cex,  type="o", letters=c(0.05,0.05))
#   mtext("Lignin/Nitrogen compounds \n peak area (%TIC, SE)", side=2, padj=-1.5)
#   axis(4, tck=0.01, labels=F)
#   mtext("litter incubation (days)", side=1, padj=3)

#   par(mar=c(2.1,1.1,0.6,5.1))
#   timeseries(orig_cTIC$L/orig_cTIC$N, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n")
#   abline(h=100, col="grey", endsig=T, letters=c(0.05,0.05))
#   #timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#   #abline(h=100, col="grey")
#   axis(2, tck=0.01, labels=F)
#   axis(4, tck=0.01)
#   mtext("Lignin/Nitrogen compounds, \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)
#   mtext("accumulated respiration (g CO2-C g-1 litter-C, SE)", side=1, padj=3)
  legend("topleft", pch=pch, typlev, cex=1.2)


  dev.off()


pdf("output/timeseries_waxes.pdf", width=12, height=16)

par(mfrow=c(4,2), oma=c(3,1,1,2))

par(mar=c(0.1,5.1,0.6,0.1))

write.csv(timeseries(class_cTIC$cut0, days, type, errcol=errcol, pch=pch, lwd=lwd, legsig=F, cex=cex, type="o", xaxt="n", endsig=T, letters=c(0.05,0.05)), "export/timeserieswaxes01.csv")
timeseries(class_cTIC$cut0, days, type, errcol=errcol, pch=pch, lwd=lwd, legsig=F, cex=cex, type="o", xaxt="n", endsig=T, letters=c(0.05,0.05))
mtext("Alcanes (C25, C27, C29) \n peak area (%TIC, SE)", side=2, padj=-1)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(4, tck=0.01, labels=F)
legend("topleft", pch=pch, typlev, cex=1.5)

par(mar=c(0.1,0.1,0.6,5.1))
write.csv(timeseries(class_cTIC$cut0, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", xaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeserieswaxes2.csv")
timeseries(class_cTIC$cut0, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", xaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
abline(h=100, col="darkgrey")
#timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, nam="Lignin/(Lignin+Carbohydrates)", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")s
axis(4, tck=0.01)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(2, tck=0.01, labels=F)
mtext("Alcanes (C25, C27, C29) \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(0.1,5.1,0.1,0.1))
write.csv(timeseries(class_cTIC$cut1, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05)), "export/timeserieswaxes3.csv")
timeseries(class_cTIC$cut1, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05))
mtext("Alcenes (C25, C27, C29) \n peak area (%TIC, SE)", side=2, padj=-1)
axis(4, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)

par(mar=c(0.1,0.1,0.1,5.1))
write.csv(timeseries(class_cTIC$cut1, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, xaxt="n", pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeserieswaxes4.csv")
timeseries(class_cTIC$cut1, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, xaxt="n", pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01)
axis(1, tck=0.01, labels=F)

mtext("Alcenes (C25, C27, C29) \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(0.1,5.1,0.1,0.1))
write.csv(timeseries(orig_cTIC$al, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05)), "export/timeserieswaxes5.csv")
timeseries(orig_cTIC$al, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05))
mtext("Unknown terpenoid (C20H40O, RT=20.00) \n peak area (%TIC, SE)", side=2, padj=-1)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)

par(mar=c(0.1,0.1,0.1,5.1))
write.csv(timeseries(orig_cTIC$al, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, xaxt="n", errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeserieswaxes6.csv")
timeseries(orig_cTIC$al, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, xaxt="n", errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(4, tck=0.01)
axis(3, tck=0.01, labels=F)
mtext("Unknown terpenoid (C20H40O, RT=20.00) \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(2.1,5.1,0.1,0.1))
write.csv(timeseries(class_cTIC$fa, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", endsig=T, letters=c(0.05,0.05)), "export/timeserieswaxes7.csv")
timeseries(class_cTIC$fa, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", endsig=T, letters=c(0.05,0.05))
mtext("Fatty acids (14:0, 16:0, 18:0) \n peak area (%TIC, SE)", side=2, padj=-1)
axis(4, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
mtext("litter incubation \n (days)", side=1, padj=1.5)

par(mar=c(2.1,0.1,0.1,5.1))
write.csv(timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeserieswaxes8.csv")
timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01)
mtext("Fatty acids (14:0, 16:0, 18:0) \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)
mtext("accumulated respiration \n (g CO2-C g-1 litter-C, SE)", side=1, padj=1.5)
legend("bottomright", pch=pch, typlev, cex=1.2)

dev.off()

pdf("output/timeseries_lignin.pdf", width=12, height=16)

par(mfrow=c(4,2), oma=c(3,1,1,2))

par(mar=c(0.1,5.1,0.6,0.1))

write.csv(timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig1.csv")
timeseries(class_cTIC$g, days, type, errcol=errcol, pch=pch, lwd=lwd, legsig=F, cex=cex, type="o", xaxt="n", endsig=T, letters=c(0.05,0.05))
mtext("Guaiacyl derrivatives \n peak area (%TIC, SE)", side=2, padj=-1)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(4, tck=0.01, labels=F)
legend("topleft", pch=pch, typlev, cex=1.5)

par(mar=c(0.1,0.1,0.6,5.1))
write.csv(timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig1.csv")
timeseries(class_cTIC$g, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", xaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
abline(h=100, col="darkgrey")
#timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, nam="Lignin/(Lignin+Carbohydrates)", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")s
axis(4, tck=0.01)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(2, tck=0.01, labels=F)
mtext("Guaiacyl derrivatives \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(0.1,5.1,0.1,0.1))
write.csv(timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig1.csv")
timeseries(class_cTIC$sy, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05))
mtext("Syringol derrivatives \n peak area (%TIC, SE)", side=2, padj=-1)
axis(4, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)

par(mar=c(0.1,0.1,0.1,5.1))
write.csv(timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig1.csv")
timeseries(class_cTIC$sy, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, xaxt="n", pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01)
axis(1, tck=0.01, labels=F)
mtext("Syringol derrivatives \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(0.1,5.1,0.1,0.1))
write.csv(timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig1.csv")
timeseries(rcTIC[,peaks$name=="G1"]/rcTIC[,peaks$name=="G0"], days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05))
mtext("Methylguaiacol:Guaiacol ratio  \n peak area (%TIC, SE)", side=2, padj=-1)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)

par(mar=c(0.1,0.1,0.1,5.1))
write.csv(timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig1.csv")
timeseries(rcTIC[,peaks$name=="G1"]/rcTIC[,peaks$name=="G0"], days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, xaxt="n", errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(4, tck=0.01)
axis(3, tck=0.01, labels=F)
mtext("Methylguaiacol:Guaiacol ratio \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(2.1,5.1,0.1,0.1))
write.csv(timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig1.csv")
timeseries(rsim[,peaks$name=="S1"]/rsim[,peaks$name=="S0"], days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", endsig=T, letters=c(0.05,0.05))
mtext("Methylsyringol:Syringol ratio \n peak area (%TIC, SE)", side=2, padj=-1)
axis(4, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
mtext("litter incubation \n (days)", side=1, padj=1.5)

par(mar=c(2.1,0.1,0.1,5.1))
write.csv(timeseries(class_cTIC$fa, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig1.csv")
timeseries(rsim[,peaks$name=="S1"]/rsim[,peaks$name=="S0"], days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01)
mtext("Methylsyringol:Syringol ratio \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)
mtext("accumulated respiration \n (g CO2-C g-1 litter-C, SE)", side=1, padj=1.5)
legend("bottomright", pch=pch, typlev, cex=1.2)

dev.off()

pdf("output/timeseries_carbohydrates.pdf", width=12, height=16)

par(mfrow=c(4,2), oma=c(3,1,1,2))

par(mar=c(0.1,5.1,0.6,0.1))

timeseries(class_cTIC$f, days, type, errcol=errcol, pch=pch, lwd=lwd, legsig=F, cex=cex, type="o", xaxt="n", endsig=T, letters=c(0.05,0.05))
mtext("Furan derrivatives \n peak area (%TIC, SE)", side=2, padj=-1)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(4, tck=0.01, labels=F)
legend("topleft", pch=pch, typlev, cex=1.5)

par(mar=c(0.1,0.1,0.6,5.1))
timeseries(class_cTIC$f, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", xaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
abline(h=100, col="darkgrey")
#timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, nam="Lignin/(Lignin+Carbohydrates)", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")s
axis(4, tck=0.01)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(2, tck=0.01, labels=F)
mtext("Furan derrivatives \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(0.1,5.1,0.1,0.1))
timeseries(class_cTIC$cp, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05))
mtext("Cyclopentenone derrivatives \n peak area (%TIC, SE)", side=2, padj=-1)
axis(4, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)

par(mar=c(0.1,0.1,0.1,5.1))
timeseries(class_cTIC$cp, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, xaxt="n", pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01)
axis(1, tck=0.01, labels=F)
mtext("Cyclopentenone derrivatives \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(0.1,5.1,0.1,0.1))
timeseries(class_cTIC$f/class_cTIC$cp, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05))
mtext("Furan:Cyclopentenone ratio  \n peak area (%TIC, SE)", side=2, padj=-1)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)

par(mar=c(0.1,0.1,0.1,5.1))
timeseries(class_cTIC$f/class_cTIC$cp, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, xaxt="n", errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(4, tck=0.01)
axis(3, tck=0.01, labels=F)
mtext("Furan:Cyclopentenone ratio \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(2.1,5.1,0.1,0.1))
timeseries(rsim[,peaks$name=="S1"]/rsim[,peaks$name=="S0"], days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", endsig=T, letters=c(0.05,0.05))
mtext("Methylsyringol:Syringol ratio \n peak area (%TIC, SE)", side=2, padj=-1)
axis(4, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
mtext("litter incubation \n (days)", side=1, padj=1.5)

par(mar=c(2.1,0.1,0.1,5.1))
timeseries(rsim[,peaks$name=="S1"]/rsim[,peaks$name=="S0"], days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01)
mtext("Methylsyringol:Syringol ratio \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)
mtext("accumulated respiration \n (g CO2-C g-1 litter-C, SE)", side=1, padj=1.5)
legend("bottomright", pch=pch, typlev, cex=1.2)

dev.off()

 
pdf("output/timeseries_orig.pdf", width=12, height=16)

par(mfrow=c(4,2), oma=c(3,1,1,2))

par(mar=c(0.1,5.1,0.6,0.1))
write.csv(timeseries(orig_cTIC$L, days, type, errcol=errcol, pch=pch, lwd=lwd, legsig=F, cex=cex, type="o", xaxt="n", endsig=T, letters=c(0.05,0.05)), "export/timeseriesorig1.csv")
timeseries(orig_cTIC$L, days, type, errcol=errcol, pch=pch, lwd=lwd, legsig=F, cex=cex, type="o", xaxt="n", endsig=T, letters=c(0.05,0.05))
mtext("Lignin \n peak area (%TIC, SE)", side=2, padj=-1)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(4, tck=0.01, labels=F)
legend("topleft", pch=pch, typlev, cex=1.5)

par(mar=c(0.1,0.1,0.6,5.1))
write.csv(timeseries(orig_cTIC$L, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", xaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig2.csv")
timeseries(orig_cTIC$L, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", xaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
abline(h=100, col="darkgrey")
#timeseries(orig_cTIC$L/(orig_cTIC$L+orig_cTIC$C), days, type, nam="Lignin/(Lignin+Carbohydrates)", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")s
axis(4, tck=0.01)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(2, tck=0.01, labels=F)
mtext("Lignin \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(0.1,5.1,0.1,0.1))
write.csv(timeseries(orig_cTIC$C, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05)), "export/timeseriesorig3.csv")
timeseries(orig_cTIC$C, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05))
mtext("Carbohydrates \n peak area (%TIC, SE)", side=2, padj=-1)
axis(4, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)

par(mar=c(0.1,0.1,0.1,5.1))
write.csv(timeseries(orig_cTIC$C, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, xaxt="n", pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig4.csv")
timeseries(orig_cTIC$C, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, xaxt="n", pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01)
axis(1, tck=0.01, labels=F)

mtext("Carbohydrates \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(0.1,5.1,0.1,0.1))
write.csv(timeseries(orig_cTIC$Ph, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05)), "export/timeseriesorig5.csv")
timeseries(orig_cTIC$Ph, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex, xaxt="n", type="o", endsig=T, letters=c(0.05,0.05))
mtext("Phenolic compounds \n peak area (%TIC, SE)", side=2, padj=-1)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)

par(mar=c(0.1,0.1,0.1,5.1))
write.csv(timeseries(orig_cTIC$Ph, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, xaxt="n", errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig6.csv")
timeseries(orig_cTIC$Ph, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, xaxt="n", errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(4, tck=0.01)
axis(3, tck=0.01, labels=F)
mtext("Phenolic compounds \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)

par(mar=c(2.1,5.1,0.1,0.1))
write.csv(timeseries(orig_cTIC$N, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", endsig=T, letters=c(0.05,0.05)), "export/timeseriesorig7.csv")
timeseries(orig_cTIC$N, days, type,  errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", endsig=T, letters=c(0.05,0.05))
mtext("Nitrogen compounds \n peak area (%TIC, SE)", side=2, padj=-1)
axis(4, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
mtext("litter incubation \n (days)", side=1, padj=1.5)

par(mar=c(2.1,0.1,0.1,5.1))
write.csv(timeseries(orig_cTIC$N, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05)), "export/timeseriesorig8.csv")
timeseries(orig_cTIC$N, days, type, masslossSE=masslossSE, massloss=massloss, normalize=1, errcol=errcol, pch=pch, lwd=lwd, legsig=F, add=F, cex=cex,  type="o", yaxt="n", endsig=F, letters=c(0.05,0.05))
abline(h=100, col="grey", lty=2 )
curve(100/(1-x), add=T, col="grey", lty=3)
#timeseries(orig_cTIC$L/orig_cTIC$N, days, type, nam="Lignin/Nitrogen compounds", masslossSE=masslossSE, massloss=massloss, normalize=2, xlab="accumulated respiration (g CO2-C g-1 litter-C, SE)", ylab="peak area (%sim)", col=c(grey(0), grey(.3), grey(.5), grey(.7)), allpoints=T, pch=19, lwd=2, legsig=F, add=F)
#abline(h=100, col="grey")
axis(2, tck=0.01, labels=F)
axis(1, tck=0.01, labels=F)
axis(3, tck=0.01, labels=F)
axis(4, tck=0.01)
mtext("Nitrogen compounds \n peak area / initial peak area (%TIC, SE)", side=4, padj=2)
mtext("accumulated respiration \n (g CO2-C g-1 litter-C, SE)", side=1, padj=1.5)
legend("bottomright", pch=pch, typlev, cex=1.2)

dev.off()

 

difinit.rsim
