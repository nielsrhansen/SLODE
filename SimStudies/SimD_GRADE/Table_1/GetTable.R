source('./help_func.R')
fns<-c("Ecoli1","Ecoli2","Yeast1","Yeast2","Yeast3")

NRep=100
ind=1
us.auc<-list()

p=10
si=1
for(fni in fns){
  gold <- read.csv( paste('./GNW/Graph_NeRDS/InSilicoSize10-',fni,'_goldstandard.tsv' ,sep=""),sep='\t',header=F)
  dt <- read.csv(paste('./GNW/Data/InSilicoSize10-',fni,'_noexpnoise_multifactorial_timeseries.tsv',sep=""),sep='\t',header=T)
  gene_names<-colnames(dt)[-1]
  A <- GNWgraph(gold,gs=gene_names) 
  diag(A)<-0
  us.auc[[ind]]<-numeric(NRep)
  for(iRandom in 1:NRep){
    load(file=paste("./RData/fig4",fni,iRandom,"US",p,"S",si, "PSFinal_Apr.RData",sep=""))  
    AUCus<-auc_us(tps=CD.PS$TP_sum,nps=CD.PS$NP_sum,p=p,ntrue=sum(A))

    us.auc[[ind]][iRandom]=AUCus$aucroc
  }
  ind<-ind+1
}

ind=1
us100.auc<-list()
p=100
si=1
for(fni in fns){
  gold <- read.csv( paste('./GNW/Graph_NeRDS/InSilicoSize100-',fni,'_goldstandard.tsv' ,sep=""),sep='\t',header=F)
  dt <- read.csv(paste('./GNW/Data/InSilicoSize100-',fni,'_noexpnoise_multifactorial_timeseries.tsv',sep=""),sep='\t',header=T)
  gene_names<-colnames(dt)[-1]
  A <- GNWgraph(gold,gs=gene_names) 
  diag(A)<-0
  
  us100.auc[[ind]]<-numeric(NRep)
  for(iRandom in 1:NRep){
    load(file=paste("./RData/fig4",fni,iRandom,"US",p,"S",si, "PSFinal.RData",sep=""))
    
    # Add the last one:
    CD.PS$TP_sum<-c(CD.PS$TP_sum)
    CD.PS$NP_sum<-c(CD.PS$NP_sum)
    
    AUCus<-auc_us(tps=CD.PS$TP_sum,nps=CD.PS$NP_sum,p=p,ntrue=sum(A))
    us100.auc[[ind]][iRandom]=AUCus$aucroc
    
  }
  ind<-ind+1
}


us.tab<-us100.tab<-matrix(0,nrow=length(fns),ncol=3)
for(i in 1:length(fns)){
  us.tab[i,]<-makeCI(us.auc[[i]])
  us100.tab[i,]<-makeCI(us100.auc[[i]])
}
pastingCIs<-function(onetab){
  paste(round(onetab[1],3), " (", round(onetab[2],3), " ", round(onetab[3],3), ")",sep="")
}
outtable<-numeric(0)
for(i in 1:length(fns)){
  temp<-c(pastingCIs(us.tab[i,]),  pastingCIs(us100.tab[i,]))
  outtable<-rbind(outtable,temp)
}

library(xtable)

rownames(outtable)<-c("Ecoli1", "Ecoli2", "Yeast1", "Yeast2", "Yeast3")
colnames(outtable)<-c("GRADE 10","GRADE 100")
xtable(outtable)
