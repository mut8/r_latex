
ord<-rda(-1*difinit.rsim, scale=T)
name<-"diff.pca"
xvar<-eigenvals(ord)/sum(eigenvals(ord))
plot(ord)
plot(xvar)


pdf("output/allpeaks_diffs_PCA12.pdf") 
      plot(ord, choices=c(1,2), type="n", tck=.01, 
      xlab=paste("PCA1", formatC(xvar[1]*100, digits=3), "% variance"),
      ylab=paste("PCA2", formatC(xvar[2]*100, digits=3), "% variance"))

      scores.sites<-(scores(ord, display="sites", choices=1:2))
#        for (i in 1:length(typlev))
#        for (j in 1:length(harlev))
#        points(ord, choices=c(1,2),  display="sites", select=harvest==harlev[j] & type==typlev[i], col=colscale[j], pch=pch[i], cex=.6)
      mat<-matrix(nrow=length(typlev)*length(harlev), ncol=6)
      colnames(mat)<-c("type", "harvest", "pca1", "pca1.error", "pca2", "pca2.error")
      for (i in 1:length(typlev))
      for (j in 1:length(harlev))
      {
      x<-mean(scores.sites[harvest==harlev[j] & type==typlev[i],1])
      y<-mean(scores.sites[harvest==harlev[j] & type==typlev[i],2])
      x.err<-stderr(scores.sites[harvest==harlev[j] & type==typlev[i],1])
      y.err<-stderr(scores.sites[harvest==harlev[j] & type==typlev[i],2])
      plotCI(x, y, uiw=y.err, liw=y.err, col=colscale[j], pch=pch[i], cex=1, add=T, gap=0)
      plotCI(x, y, uiw=x.err, liw=x.err, err="x", col=colscale[j], pch=pch[i], cex=1, add=T, gap=0)
      print(paste(x, x.err,y, y.err))
      mat[length(harlev)*(i-1)+j, T]<-c(typlev[i], harlev[j], x,x.err,y,y.err)
      }
      mat
      write.csv(mat, "export/dif.sites.csv")
      scores<-(scores(ord, display="species", choices=1:2))
      text(scores[,1]*2.5, scores[,2]*2.5, labels=peaks$orig, cex=.4)
      write.csv(data.frame(scores,peaks$orig), "export/dif.species.csv")
      title(name)

dev.off()

ord<-rda(rsim, scale=T)
name<-"allpeaks.pca"
xvar<-eigenvals(ord)/sum(eigenvals(ord))

plot(ord)
plot(xvar)


pdf("output/allpeaks_PCA12.pdf") 
      plot(ord, choices=c(1,2), type="n", tck=.01, 
      xlab=paste("PCA1", formatC(xvar[1]*100, digits=3), "% variance"),
      ylab=paste("PCA2", formatC(xvar[2]*100, digits=3), "% variance"))

      scores.sites<-(scores(ord, display="sites", choices=1:2))
#        for (i in 1:length(typlev))
#        for (j in 1:length(harlev))
#        points(ord, choices=c(1,2),  display="sites", select=harvest==harlev[j] & type==typlev[i], col=colscale[j], pch=pch[i], cex=.6)
      mat<-matrix(nrow=length(typlev)*length(harlev), ncol=6)
      colnames(mat)<-c("type", "harvest", "pca1", "pca1.error", "pca2", "pca2.error")
      for (i in 1:length(typlev))
      for (j in 1:length(harlev))
      {
      x<-mean(scores.sites[harvest==harlev[j] & type==typlev[i],1])
      y<-mean(scores.sites[harvest==harlev[j] & type==typlev[i],2])
      x.err<-stderr(scores.sites[harvest==harlev[j] & type==typlev[i],1])
      y.err<-stderr(scores.sites[harvest==harlev[j] & type==typlev[i],2])
      plotCI(x, y, uiw=y.err, liw=y.err, col=colscale[j], pch=pch[i], cex=1, add=T, gap=0)
      plotCI(x, y, uiw=x.err, liw=x.err, err="x", col=colscale[j], pch=pch[i], cex=1, add=T, gap=0)
      print(paste(x,y,x.err, y.err))
      mat[length(harlev)*(i-1)+j, T]<-c(typlev[i], harlev[j], x,x.err,y,y.err)
      }
      write.csv(mat, "export/all.sites.csv")
      scores<-(scores(ord, display="species", choices=1:2))
      text(scores[,1]*2.5, scores[,2]*2.5, labels=peaks$orig, cex=.4)
      write.csv(data.frame(scores,peaks$orig), "export/all.species.csv")
      title(name)

dev.off()


ord<-cca(rsim, scale=T)
name<-"allpeaks.cca"
xvar<-eigenvals(ord)/sum(eigenvals(ord))

plot(ord)
plot(xvar)


pdf("output/allpeaks_CCA12.pdf") 
      plot(ord, choices=c(1,2), type="n", tck=.01, 
      xlab=paste("CCA1", formatC(xvar[1]*100, digits=3), "% variance"),
      ylab=paste("CCA2", formatC(xvar[2]*100, digits=3), "% variance"))

      scores.sites<-(scores(ord, display="sites", choices=1:2))
#        for (i in 1:length(typlev))
#        for (j in 1:length(harlev))
#        points(ord, choices=c(1,2),  display="sites", select=harvest==harlev[j] & type==typlev[i], col=colscale[j], pch=pch[i], cex=.6)

      for (i in 1:length(typlev))
      for (j in 1:length(harlev))
      {
      x<-mean(scores.sites[harvest==harlev[j] & type==typlev[i],1])
      y<-mean(scores.sites[harvest==harlev[j] & type==typlev[i],2])
      x.err<-stderr(scores.sites[harvest==harlev[j] & type==typlev[i],1])
      y.err<-stderr(scores.sites[harvest==harlev[j] & type==typlev[i],2])
      plotCI(x, y, uiw=y.err, liw=y.err, col=colscale[j], pch=pch[i], cex=1, add=T, gap=0)
      plotCI(x, y, uiw=x.err, liw=x.err, err="x", col=colscale[j], pch=pch[i], cex=1, add=T, gap=0)
      print(paste(x,y,x.err, y.err))
      }
      scores<-(scores(ord, display="species", choices=1:2))
      text(scores[,1]*2.5, scores[,2]*2.5, labels=peaks$class, cex=.4)
      title(name)

dev.off()

trans.difinit.rsim<-difinit.rsim
for (i in 1:peaknr)
trans.difinit.rsim[,i]<-(difinit.rsim[,i]-(min(difinit.rsim[,i])))/max(difinit.rsim[,i])

ord<-cca(trans.difinit.rsim)
name<-"diff.cca"
xvar<-eigenvals(ord)/sum(eigenvals(ord))
plot(ord)
plot(xvar)


pdf("output/allpeaks_diffs_CCA12.pdf") 
      plot(ord, choices=c(1,2), type="n", tck=.01, 
      xlab=paste("CCA1", formatC(xvar[1]*100, digits=3), "% variance"),
      ylab=paste("CCA2", formatC(xvar[2]*100, digits=3), "% variance"))

      scores.sites<-(scores(ord, display="sites", choices=1:2))
#        for (i in 1:length(typlev))
#        for (j in 1:length(harlev))
#        points(ord, choices=c(1,2),  display="sites", select=harvest==harlev[j] & type==typlev[i], col=colscale[j], pch=pch[i], cex=.6)

      for (i in 1:length(typlev))
      for (j in 1:length(harlev))
      {
      x<-mean(scores.sites[harvest==harlev[j] & type==typlev[i],1])
      y<-mean(scores.sites[harvest==harlev[j] & type==typlev[i],2])
      x.err<-stderr(scores.sites[harvest==harlev[j] & type==typlev[i],1])
      y.err<-stderr(scores.sites[harvest==harlev[j] & type==typlev[i],2])
      plotCI(x, y, uiw=y.err, liw=y.err, col=colscale[j], pch=pch[i], cex=1, add=T, gap=0)
      plotCI(x, y, uiw=x.err, liw=x.err, err="x", col=colscale[j], pch=pch[i], cex=1, add=T, gap=0)
      print(paste(x,y,x.err, y.err))
      }
      scores<-(scores(ord, display="species", choices=1:2))
      text(scores[,1]*2.5, scores[,2]*2.5, labels=peaks$class, cex=.4)
      title(name)
dev.off()
