rm(list = ls())

source("Sources/Smoother.R")
source("Sources/Main.R")
source("Sources/Estimator.R")
source("Sources/GroupCorD.R")
source("Sources/Recoverer.R")
source("Sources/Regressor.R")
source("Sources/Lineargen.R")
source("Sources/Generator.R")
source("Sources/Wrapper.R")
source("Sources/GroupGrp.R")
library(splines,MASS)
library(grplasso) #install.packages('grplasso');
library(locfit)#install.packages('locfit');
library(fda)#install.packages('fda');


RandomSeed=1;


#---------------------------------------------#
# Generate a simple graph
source("Example/Graph_generation.R")
image(tgraph)

#---------------------------------------------#
# Generation noiseless and noisy trajectories 
source("Example/Data_generation.R")
# The simulated data is stored in a list Y with R elements.
# Each of the R elements is an n by p matrix. 

# Noiseless trajectories 
plot(0,0,xlim=c(0,upp),ylim=range(solutions$mean_curve[[1]]))
for(j in 1:p){
  for(r in 1:R){
    lines(solutions$mean_curve[[r]][,j]~tgrid,col=j)
  }
}


# Noisy observations 
plot(0,0,xlim=c(0,upp),ylim=range(Y[[1]]))
for(j in 1:p){
  for(r in 1:R){
    lines(Y[[r]][,j]~times,col=j)
  }
}


#--------------------------------------------#
# Edge recovery performance 
source("Example/Analysis.R")

plot(0,0,pch=16,col="white",type="l",lwd=3,xlim=c(0,p^2),ylim=c(0,p), cex.lab=1,cex.main=1,ylab='True positives',xlab='Total positives',main='Edge recovery')
lines(recint$TP_sum~recint$NP_sum,col='black',lwd=2)






