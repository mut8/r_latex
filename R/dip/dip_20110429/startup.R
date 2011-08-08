errcol<-"black"
lwd<-.3
cex<-2
lty=1
pch<-c(1,16,2,17)
pch.all<-c(rep(pch[1],19), rep(pch[2],18), rep(pch[3],19), rep(pch[4],18))

addpeaks<-function(source, sep=",") {

for (i in 1:length(source))
{
print(i)
tmp<-read.csv(source[i], sep=sep, header=F)
pe<-data.frame(t(tmp[1:10,2:ncol(tmp)]))
colnames(pe)<-tmp[1:10, 1]
rownames(pe)<-paste(pe$class, 1:nrow(pe), sep="")

si<-tmp[12:85, 2:ncol(tmp)]
colnames(si)<-paste(rownames(pe), "sim", sep="")
rownames(si)<-rownames(samples) 

ti<-tmp[102:175, 2:ncol(tmp)]
colnames(ti)<-paste(rownames(pe), "tic", sep="")
rownames(ti)<-rownames(samples) 

if (i==1) 
{
peaks<-pe
sim<-si
tic<-ti
} else 
{
peaks<-rbind(peaks, pe)
sim<-data.frame(sim, si)
tic<-data.frame(tic, ti)
}
}
write.csv(peaks, "raw_data/peaks.csv")
write.csv(sim, "raw_data/sim.csv")
write.csv(tic, "raw_data/tic.csv")
print("files added")
}

samples<-read.csv("raw_data/samples.csv", sep=";", header=TRUE)
addpeaks(c("raw_data/furane.csv","raw_data/aromates.csv", "raw_data/phenole.csv", "raw_data/carbos_neu.csv", "raw_data/cyclopentenone.csv", "raw_data/n_compounds_neu.csv", "raw_data/alkan_en_cleaned.csv", "raw_data/the_rest.csv", "raw_data/lignin2.csv", "raw_data/keto-alkyl_alcohole.csv", "raw_data/fa.csv"), sep=",")

minerals<-read.csv("raw_data/minerals.csv", sep=",", header=TRUE)
alldata<-read.csv("raw_data/micdif_alldata.csv", sep=",", header=TRUE)
alldata$days<-rep(c(rep(14,5), rep(97,5), rep(181,5), rep(375,5)),4)



sim<-read.csv("raw_data/sim.csv", sep=",", header=TRUE)
rownames(sim)<-sim[,1]
sim[,1]<-NULL
peaks<-read.csv("raw_data/peaks.csv", sep=",", header=TRUE)

corr<-read.csv("raw_data/corr_data.csv")
corr$pH.CaCl2<-NULL

#remove excluded
sim<-sim[,peaks$exclude==F]
peaks<-peaks[peaks$exclude==F,T]

colnames(peaks[10])<-"orig2"

peaks$origin[peaks$orig2=="C"]<-"Carbohydrates"
peaks$origin[peaks$orig2=="L"]<-"Lignin"
peaks$origin[peaks$orig2=="N"]<-"Nitrogen"
peaks$origin[peaks$orig2=="Cut"]<-"Waxes"
peaks$origin[peaks$orig2=="non"]<-"non-specific"
peaks$origin[peaks$orig2=="unk"]<-"unknown"

peaks$origin<-factor(peaks$origin)

days<-samples$day
accresp<-samples$accresp
harvest<-samples$harvest
type<-samples$type
harlev<-as.numeric(levels(as.factor(harvest)))
typlev<-levels(type)

harvest

# carbos<-peaks[,2]=="c"|peaks[,2]=="ch"|peaks[,2]=="f"
# ligs<-peaks[,2]=="sy"|peaks[,2]=="g"|peaks[,2]=="h"
# #ligs<-peaks[,2]=="sy"|peaks[,2]=="g"|peaks[,2]=="bf"|peaks[,2]=="h"
# nitr<-peaks[,2]=="ind"|peaks[,2]=="pr"
# phen<-peaks[,2]=="ph"
# al<-peaks[,2]=="al"
# ka<-peaks[,2]=="ka"
# unk<-peaks[,2]=="u"

peaknr<-length(colnames(sim))

      # code2<-vector(length=peaknr)
      # code2[ligs]<-"L"
      # code2[nitr]<-"Pr"
      # code2[carbos]<-"C"
      # code2[phen]<-"Ph"
      # code2[peaks$origin=="a"]<-"a"
      # code2[al]<-"al"
      # code2[ka]<-"ka"
      # code2[peaks$origin=="s"]<-"s"
      # code2[peaks$origin=="t"]<-"t"
      # code2[peaks$origin=="bf"]<-"bf"
      # code2[unk]<-"u"
      # 
      # code3<-vector(length=peaknr)
      # code3[ligs]<-paste("L", 1:sum(ligs), sep="")
      # code3[nitr]<-paste("Pr", 1:sum(nitr), sep="")
      # code3[carbos]<-paste("C", 1:sum(carbos), sep="")
      # code3[phen]<-paste("Ph", 1:sum(phen), sep="")
      # code3[peaks$origin=="a"]<-paste("a", 1:sum(peaks$origin=="a"), sep="")
      # code3[al]<-paste("al", 1:sum(al), sep="")
      # code3[ka]<-paste("ka", 1:sum(ka), sep="")
      # code3[peaks$origin=="s"]<-paste("s", 1:sum(peaks$origin=="s"), sep="")
      # code3[peaks$origin=="t"]<-paste("t", 1:sum(peaks$origin=="t"), sep="")
      # code3[peaks$origin=="bf"]<-paste("t", 1:sum(peaks$origin=="bf"), sep="")
      # code3[unk]<-paste("u", 1:sum(unk), sep="")

      # peaks<-data.frame(peaks, code2, code3)

colnames(peaks)

# calculate cTIC
    cTIC<-sim
    for (i in 1:peaknr)
	 cTIC[,i]<-sim[,i] * peaks$tic[i]
 cTIC

#subs by categories
 colnames(peaks)

      rsim<-100*sim/rowSums(sim)
      rcTIC<-100*cTIC/rowSums(cTIC)
      class_rsim<-sumif(rsim,peaks$class)
      orig_rsim<-sumif(rsim,peaks$origin)
      class_cTIC<-sumif(rcTIC,peaks$class)
      orig_cTIC<-sumif(rcTIC,peaks$origin)



samples$massloss->massloss
samples$masslossSE->masslossSE

processes<-data.frame(corr[,38], corr[,3:30], corr[,51:56], corr[,61:69])
controls<-data.frame(class_cTIC[harvest==2|harvest==6, T], orig_cTIC[harvest==2|harvest==6, T],  corr[,31:34], corr[,39:44], corr[,35:37], corr[,57:59], corr[,45:50], corr[,60])

contnames<-colnames(controls)
procnames<-colnames(processes)

h2<-corr$harvest==2
h3<-corr$harvest==6
h23<-harvest==2|harvest==6
h234<-harvest==2|harvest==6|harvest==15
h34<-harvest==6|harvest==15
h4<-harvest==15

#h2.mds<-metaMDS(rsim[harvest==2,1:ncol(rsim)])
#h3.mds<-metaMDS(rsim[harvest==6,1:ncol(rsim)])
#h23.mds<-metaMDS(rsim[harvest==6|harvest==2,1:ncol(rsim)])
#all.mds<-metaMDS(rsim)

#h2.pca<-rda(rsim[harvest==2,1:ncol(rsim)], scale=TRUE)
#h3.pca<-rda(rsim[harvest==6,1:ncol(rsim)], scale=TRUE)
#h23.pca<-rda(rsim[harvest==6|harvest==2,1:ncol(rsim)], scale=TRUE)
#all.pca<-rda(rsim, scale=TRUE)
#h23.wo.n.pca<-rda(rsim[h23, 1:ncol(rsim)],,controls[,52], scale=TRUE)
#h23.wo.p.pca<-rda(rsim[h23, 1:ncol(rsim)],,controls[,53], scale=TRUE)
#h23.wo.n.p.pca<-rda(rsim[h23, 1:ncol(rsim)],,controls[,52:53], scale=TRUE)
#h2.wo.n.pca<-rda(rsim[harvest==2, 1:ncol(rsim)],,controls[h2,52], scale=TRUE)
#h2.wo.p.pca<-rda(rsim[harvest==2, 1:ncol(rsim)],,controls[h2,53], scale=TRUE)
#h2.wo.n.p.pca<-rda(rsim[harvest==2, 1:ncol(rsim)],,controls[h2,52:53], scale=TRUE)
#h3.wo.n.pca<-rda(rsim[harvest==6, 1:ncol(rsim)],,controls[h3,52], scale=TRUE)
#h3.wo.p.pca<-rda(rsim[harvest==6, 1:ncol(rsim)],,controls[h3,53], scale=TRUE)
#h3.wo.n.p.pca<-rda(rsim[harvest==6, 1:ncol(rsim)],,controls[h3,52:53], scale=TRUE)


 	  rec<-matrix(ncol=5, nrow=peaknr) 
 	  colnames(rec)<-c(typlev, "mean")
 	  rownames(rec)=rownames(peaks)
 	  for (i in 1:peaknr)
 	  {
 	  for (j in 1:length(typlev))
 	  {

 	  tmp<-rsim[type==typlev[j]&harvest!=0,i]/mean(rsim[type==typlev[j]&harvest!=0,i])
 	  mod<-lm(tmp ~ samples$massloss[type==typlev[j]&harvest!=0])
 	  rec[i,j]<-mod$coefficients[2]
 	  }
 	  rec[i,5]<-mean(rec[i, 1:4])
 	  }
 	  rec

# 	  rec2<-matrix(ncol=5, nrow=peaknr) 
# 	  colnames(rec2)<-c(typlev, "mean")
# 	  rownames(rec2)=rownames(peaks)
# 	  for (i in 1:peaknr)
# 	  {
# 	  for (j in 1:length(typlev))
# 	  {
# 	  rsim[type==typlev[j]&harvest!=0,i]
# 	  tmp<-rsim[type==typlev[j]&harvest!=0,i]*(100-samples$massloss[type==typlev[j]&harvest!=0])/mean(rsim[type==typlev[j]&harvest==0,i])
# 	  mod<-lm(tmp ~ samples$massloss[type==typlev[j]&harvest!=0])
# 	  rec2[i,j]<-mod$coefficients[2]
# 	  }
# 	  rec2[i,5]<-mean(rec2[i, 1:4])
# 	  }
# 	  rec2
# 
# 	  rec.ctic<-matrix(ncol=5, nrow=peaknr) 
# 	  colnames(rec.ctic)<-c(typlev, "mean")
# 	  rownames(rec.ctic)=rownames(peaks)
# 	  for (i in 1:peaknr)
# 	  {
# 	  for (j in 1:length(typlev))
# 	  {
# 	  tmp<-rcTIC[type==typlev[j]&harvest!=0,i]/mean(rcTIC[type==typlev[j]&harvest==0,i])
# 	  mod<-lm(tmp ~ samples$massloss[type==typlev[j]&harvest!=0])
# 	  rec.ctic[i,j]<-mod$coefficients[2]
# 	  }
# 	  rec.ctic[i,5]<-mean(rec.ctic[i, 1:4])
# 	  }
# 	  rec.ctic
# 
# 	  h234<-harvest==2|harvest==6|harvest==15
# 
# 	  rec2.ctic<-matrix(ncol=5, nrow=peaknr) 
# 	  colnames(rec2.ctic)<-c(typlev, "mean")
# 	  rownames(rec2.ctic)=rownames(peaks)
# 	  for (i in 1:peaknr)
# 	  {
# 	  for (j in 1:length(typlev))
# 	  {
# 	  rsim[type==typlev[j]&harvest!=0,i]
# 	  tmp<-rsim[type==typlev[j]&harvest!=0,i]*(100-samples$massloss[type==typlev[j]&harvest!=0])/mean(rcTIC[type==typlev[j]&harvest==0,i])
# 	  mod<-lm(tmp ~ samples$massloss[type==typlev[j]&harvest!=0])
# 	  rec2.ctic[i,j]<-mod$coefficients[2]
# 	  }
# 	  rec2.ctic[i,5]<-mean(rec2.ctic[i, 1:4])
# 	  }
# 
# 	  rsim.colsums<-colSums(rsim)
# 
 	  recidx<-matrix(ncol=1, nrow=nrow(samples)) 
 	  for (i in 1:nrow(samples))
 	  recidx[i]<-sum(t(rcTIC[i, 1:peaknr])*rec[,5])
	 
# 	  recidx.wo.n<-matrix(ncol=1, nrow=nrow(samples)) 
# 	  for (i in 1:nrow(samples))
# 	  recidx.wo.n[i]<-sum(t(rsim[i, peaks$origin!="p"&peaks$origin!="ind"])*rec[peaks$origin!="p"&peaks$origin!="ind",5])
# 
# 	  recidx2<-matrix(ncol=1, nrow=nrow(samples)) 
# 	  for (i in 1:nrow(samples))
# 	  recidx2[i]<-sum(t(rsim[i, 1:peaknr])*rec2[,5])/1000
# 
# 	  recidx2.wo.n<-matrix(ncol=1, nrow=nrow(samples)) 
# 	  for (i in 1:nrow(samples))
# 	  recidx2.wo.n[i]<-sum(t(rsim[i, peaks$origin!="p"&peaks$origin!="ind"])*rec2[peaks$origin!="p"&peaks$origin!="ind",5])/1000
# 
# 	  recidx.ctic<-matrix(ncol=1, nrow=nrow(samples)) 
# 	  for (i in 1:nrow(samples))
# 	  recidx.ctic[i]<-sum(t(rcTIC[i, 1:peaknr])*rec.ctic[,5])
# 
# 	  recidx.wo.n.ctic<-matrix(ncol=1, nrow=nrow(samples)) 
# 	  for (i in 1:nrow(samples))
# 	  recidx.wo.n.ctic[i]<-sum(t(rcTIC[i, peaks$origin!="p"&peaks$origin!="ind"])*rec.ctic[peaks$origin!="p"&peaks$origin!="ind",5])
# 
# 	  recidx2.ctic<-matrix(ncol=1, nrow=nrow(samples)) 
# 	  for (i in 1:nrow(samples))
# 	  recidx2.ctic[i]<-sum(t(rsim[i, 1:peaknr])*rec2[,5])/1000
# 
# 	  recidx2.wo.n.ctic<-matrix(ncol=1, nrow=nrow(samples)) 
# 	  for (i in 1:nrow(samples))
# 	  recidx2.wo.n.ctic[i]<-sum(t(rcTIC[i, peaks$origin!="p"&peaks$origin!="ind"])*rec2.ctic[peaks$origin!="p"&peaks$origin!="ind",5])/1000


minerals.stat<-data.frame(matrix(nrow=length(typlev), ncol=2*(ncol(minerals[,2:ncol(minerals)]))))
colnames(minerals.stat)<-paste(c(rep("means", ncol(minerals[,2:ncol(minerals)])), rep("se", ncol(minerals[,2:ncol(minerals)]))), rep(colnames(minerals[,2:ncol(minerals)]),2))
rownames(minerals.stat)<-typlev

for (i in 1:length(typlev))
for (j in 1:(ncol(minerals)-1))
{
minerals.stat[i,j]<-mean(minerals[minerals[,1]==typlev[i], j+1])
minerals.stat[i,j+(ncol(minerals)-1)]<-stderr(minerals[minerals[,1]==typlev[i], j+1])
}


alldata.stat<-data.frame(matrix(nrow=length(typlev)*length(levels(alldata$Harvest)), ncol=2+ncol(alldata[3:ncol(alldata)])*2))
colnames(alldata.stat)<-c("type", "harvest", paste(c(rep("means", ncol(alldata[,3:ncol(alldata)])), rep("se", ncol(alldata[,3:ncol(alldata)]))), rep(colnames(alldata[,3:ncol(alldata)]),2)))
alldata.stat$type<-sort(rep(typlev, length(harlev)))
alldata.stat$harvest<-rep(levels(alldata$Harvest), length(typlev))


for (i in 1:nrow(alldata.stat))
for (k in 3:(ncol(alldata)))
{
alldata.stat[i,k] <- 
mean(alldata[alldata$Litter==alldata.stat$type[i] & alldata$Harvest==alldata.stat$harvest[i], k])
alldata.stat[i,k-1+(ncol(alldata.stat)/2)] <- 
stderr(alldata[alldata$Litter==alldata.stat$type[i] & alldata$Harvest==alldata.stat$harvest[i], k])
}

colscale<-c(grey(0), grey(.3), grey(.5), grey(.7))
controls$recidx<-recidx[harvest==2 | harvest==6]
controls$L_PH<-controls$L+controls$Ph

pch<-c(1,16,2,17)

data<-rsim
cond<-harvest==0
by<-type




# 
# for (i in 1:peaknr)
# 	 cTIC[,i]<-sim[,i] * peaks$tic[i]
# 
#       rsim<-100*sim/rowSums(sim)
#       rcTIC<-100*cTIC/rowSums(cTIC)
#       class_rsim<-sumif(rsim,peaks$class)
#       orig_rsim<-sumif(rsim,peaks$origin)
#       class_cTIC<-sumif(rcTIC,peaks$class)
#       orig_cTIC<-sumif(rcTIC,peaks$origin)

alldata$CN_inbal<-alldata$C.N_lit/alldata$C.N_mic
alldata$CP_inbal<-alldata$C.P_lit/alldata$C.Pmic
alldata$NP_inbal<-alldata$N.P_lit/alldata$C.Pmic
