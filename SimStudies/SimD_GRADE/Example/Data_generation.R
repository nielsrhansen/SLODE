##---------------------------------------------------------------##
#### Data generation #####
R<-5; type_data="Perturbations";sigma<-sqrt(1); # standard deviation of the noise 
n<- 40;upp<-1 # number of samples, number of variables,range of time 

#R<-1; type_data="Replicates";
tgrid<-seq(from=0,to=upp,by=0.001)

# Find solutions of the ODE system given different initial values 
solutions<-list()
solutions$mean_curve<-list()
set.seed(RandomSeed) # Set seed for reproducibility 
for(r in 1:R){
  inis<-rnorm(p/2)
  initialvalue<-c(rbind(sin(inis), cos(inis)))
  Xs<-EulerM(initialvalue,range=c(0,upp),gap=0.001,f=LinearODE,theta=thetas)
  solutions$mean_curve[[r]]<-Xs 
}
solutions$times_fine<-tgrid
solutions$gap<-0.001
solutions$graph=tgraph
solutions$theta<-thetas

# Add noise to the noiseless trajectories 
times<-seq(from=upp/n,to=upp,length.out=n) # Observed time points 
times_e<-seq(from=0.001,to=upp,by=0.001)# dense time grid for estimating the integrals 
observations<-getObservations(times=times,solutions=solutions,sigma=sigma,R=R,type_data=type_data)

# Remove the initial values from the noisy observations 
Y<-list();
for(r in 1:R){
  Y[[r]] <- observations[[r]][,-1]
  
}
  
