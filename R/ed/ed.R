setwd("/home/lluc/Documents/R/ed")
tea<-read.csv("tea_peakareas.csv")

sink("output.txt")
names<-colnames(tea)
for (i in 4:ncol(tea))
{
print(names[i])
aov<-anova(tea[,i] ~ tea[,2])
print(summary(aov))
names(aov)

aov$qr
#print(TukeyHSD(aov))

}
sink()