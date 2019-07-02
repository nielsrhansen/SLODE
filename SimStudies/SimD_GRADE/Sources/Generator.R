
###----------------------------------------###
## The Euler-Maruyama method
## Input: 
##       Initial.value
##       range: the range of time to consider
##       gap: the time gap (default is 0.01)
##       f: the right-hand side of the ODE
## Output: 
##       A matrix of n times p, where p is the length of Initial.value and n is the length of Time.points
## Note: Length.N-1 should be the multiplicative of N-1
EulerM<-function(Initial_value=c(-1,1), range=c(0,20), gap=0.001, f, theta,... ){
  T=(range[2]-range[1])/gap + 1
  P<-length(Initial_value)
  out_full<-matrix(0,T,P)
  out_full[1,]<-Initial_value
  
  for(i in 1:(T-1)){
    out_full[i+1,]<-out_full[i,]+f(out_full[i,],theta)*gap
  }
  return(out_full)
}

FitzHugh_Nagumo<-function(X,theta){
  X1Prime<-theta[3]*(X[1]-X[1]^3/3+X[2] )
  X2Prime<- -(X[1]-theta[1]+theta[2]*X[2])/theta[3]
  return( c(X1Prime, X2Prime))
}

# Here theta is a 7 by 2 vector

FitzHugh_Nagumo_extend<-function(X,theta){
  X1Prime<-theta[1,1]+theta[1,2]*X[1]+theta[1,3]*X[1]^2+
    theta[1,4]*X[1]^3+theta[1,5]*X[2]+theta[1,6]*X[2]^2+theta[1,7]*X[2]^3
  X2Prime<-theta[2,1]+theta[2,2]*X[1]+theta[2,3]*X[1]^2+
    theta[2,4]*X[1]^3+theta[2,5]*X[2]+theta[2,6]*X[2]^2+theta[2,7]*X[2]^3
  return( c(X1Prime, X2Prime))
}

###
## Adding noise to the mean curves
## as in Wu et al. 2013
## Input: 
##      X: the data matrix
##      times:: sequence of time to measure 
##      sigma: standard deviation of noise

corrupt<-function(X,times,sigma,gap=0.001){
  Y<-X[times/gap+1,]
  total_samples<-dim(Y)[1]*dim(Y)[2]
  noise<-matrix(rnorm(total_samples, mean=0,sd=sigma), nrow=dim(Y)[1], ncol =dim(Y)[2] )
  Y<-Y+noise
  return(Y)
}


# Let the noise variables be smoothed random walks.
# initial: initial values
# variance: the variance coefficient, default is 1
# time_grid: the time grid to consider
randomWalk<-function(initial,variance=1,time_grid){
  
  out<-numeric(length(time_grid))
  out[1]<-initial
  for(i in 2:length(time_grid)){
    out[i]<-out[i-1]+rnorm(1,mean=0,sd=variance*(time_grid[i]-time_grid[i-1]) )
  }
  return(out)
}

# Adding simple ODE noise variables
simple_ode<-function(time_grid,c,iv){
  out<-time_grid*c+iv
  return(out)
}


## Scaler: to scale the functions to have a certain range
## vec: target vector
## range: the target range (max-min)
scaler<-function(vec,range){
  out<-vec*range/(max(vec)-min(vec))
  out<-out-mean(out)
  return(out)
}
