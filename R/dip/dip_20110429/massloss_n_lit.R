pdf("massloss_n_lit.pdf")

cor<-cor.test(samples$h4.massloss[h4], samples$N_lit[h4])
cor
plot(samples$N_lit[h4], samples$h4.massloss[h4], tck=.02, xlab="litter N content", ylab="mass loss (15 month)", pch=19)
if(as.numeric(cor$p.value)<0.05)
abline(lm(samples$h4.massloss ~ samples$N_lit))

dev.off()
