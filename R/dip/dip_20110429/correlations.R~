h2.pro<-harvest[harvest==2 | harvest==6]==2
h3.pro<-harvest[harvest==2 | harvest==6]==6

tmp<-data.frame(processes, controls)

cond<-h2.pro
write.csv(corr.ab(tmp[cond], tmp[cond]),"proc_cont_cor_h2.csv")
cond<-h3.pro
write.csv(corr.ab(tmp[cond], tmp[cond]),"proc_cont_cor_h3.csv")

