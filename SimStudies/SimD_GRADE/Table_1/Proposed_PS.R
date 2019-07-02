# Loading packages and functrions
source("./help_func.R")
source("../Sources/Smoother.R")
source("../Sources/Main.R")
source("../Sources/Estimator.R")
source("../Sources/GroupCorD.R")
source("../Sources/Recoverer.R")
source("../Sources/Regressor.R")

library(fda)

############################
# For simulation
#args <- commandArgs(TRUE)
#RandomSeed <- as.numeric(args[[1]])


# For debugging:
RandomSeed <- 1



s=0.025 # Standard deviation of noise
si = 1

fns<-c("Ecoli1","Ecoli2","Yeast1","Yeast2","Yeast3")

for(fni in fns){
  
  # Read the GNW data sets
  set.seed(RandomSeed)
  gold <- read.csv( paste('./GNW/Graph_NeRDS/InSilicoSize10-',fni,'_goldstandard.tsv' ,sep=""),sep='\t',header=F)
  #dtold <- read.csv(paste('./GNW/Data_NeRDS/InSilicoSize10-',fni,'_noexpnoise_multifactorial_timeseries.tsv',sep=""),sep='\t',header=T)
  dt <- read.csv(paste('./GNW/Data/InSilicoSize10-',fni,'_noexpnoise_multifactorial_timeseries.tsv',sep=""),sep='\t',header=T)
  
  R <- 10;d<-p <- 10;n <- 21
  gene_names<-colnames(dt)[-1]
  A <- GNWgraph(gold,gs=gene_names) 
  X <- list()
  for(r in 1:R){
    X[[r]] <- dt[seq(n*(r-1)+1,n*r,1),2:(d+1)]
  }
  
  ## Add noise and normalize ##
  Y<-list()
  for(r in 1:R){
    Y[[r]] <- X[[r]] + rnorm(n*d,0,s)
    Y[[r]][Y[[r]]<0] <- 0
    Y[[r]] <- Y[[r]] / max(Y[[r]])
  }
  
  # Obtain smooth estimates
  times<-tm<-seq(0,1,1/{n-1})
  type_data="Perturbations"
  type_reg="Y";
  type_smooth="smoothspline";
  maxk=5000;
  n.knots=4;nspecs=c(4,n.knots-2)
  
  times_e<-seq(from=0,to=1,by=0.001)
  lambda_N<- 100
  observations<-list()
  for(r in 1:R){
    observations[[r]]<- cbind(tm,Y[[r]])
  }
  
  smthed<-smoothX(observations=observations,times_e=times_e,deg=NULL,h=NULL,type_smooth=type_smooth,type_data=type_data)
  xhat<-smthed$Xhat
  
  # Fit the penalized regression
  lambda_range<-c(-5.5,-1)
  lambdas<-exp(seq(from=lambda_range[2],to=lambda_range[1],length.out=lambda_N)) 
  lambdasCD<-lambdas*R
  L1s<-0.1
  
  fits.fun.PS<-GRADE(observations=observations, times_e=times_e,type_reg=type_reg, type_data=type_data,type_basis="penalized spline", type_smooth=type_smooth, xhat=xhat, xprime=NULL,lambdas=lambdasCD, lambdas.fun=L1s,nspecs=nspecs)
  CD.PS<-countGraph(graph_est=fits.fun.PS$neighbour.path,graph_true=t(A),self=F)
  print(fni)
  save(CD.PS, file=paste("./RData/fig4",fni,RandomSeed,"US",p,"S", si, "PSFinal_Apr.RData",sep=""))
}

