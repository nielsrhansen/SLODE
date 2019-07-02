
#######################
# Data generation code
EulerM<-function(Initial_value=NULL, range=NULL, gap=0.001, f, theta,... ){
  T=(range[2]-range[1])/gap + 1
  P<-length(Initial_value)
  out_full<-matrix(0,T,P)
  out_full[1,]<-Initial_value
  
  for(i in 1:(T-1)){
    out_full[i+1,]<-out_full[i,]+f(out_full[i,],theta)*gap
  }
  return(out_full)
}

LinearODE<-function(X,theta){
  #theta=(alpha,beta,gamma, dealta)
  XPrime<- theta[,-1]%*%X+theta[,1]
  return( XPrime)
  }



#######################

###############
# Example data:
# 
# upp<-2*pi
# n=250
# tgrid<-seq(from=0,to=upp,by=0.001)
# 
# set.seed(1)
# p<-10
# 
# initialvalue<- {runif(p)-0.5}*2 
# ks<-seq(from=1,to= p/2,length.out=p/2)
# 
# thetas<-matrix(0,p,p)
# 
# for(j in 1:{p/2}){
#   thetas[ 2*j-1  ,2*j]<- ks[j]
#   thetas[2*j, 2*j-1]<- {-ks[j]}
# }
# thetas<-cbind(rep(0,p),thetas)
# 
# # initial values are pertubed around SP and bounded away from zero
# Xs<-EulerM(initialvalue,range=c(0,upp),gap=0.001,f=LinearODE,theta=thetas)
# #Xprime<-t(apply(Xs,1,Lotka_Volterra,theta=c(1,0.2,0.2,1))
# plot(Xs[,1]~tgrid,type="l",xlab="t",ylab="X",ylim=c(min(Xs),max(Xs)),col='white')
# for(j in 1:p){
#   lines(Xs[,j]~tgrid, col="red")  
#   }
# 
# Xprime = Xs
# for(i in 1:{dim(Xs)[1]}){
#   Xprime[i,]<-LinearODE(X=Xs[i,], theta=thetas )  
# }
# 
# plot(Xprime[,1]~tgrid,type="l",lty=2,xlab="t", ylab="X'",col='red',ylim=range(Xprime))
# for(j in 1:p){
#   lines(Xprime[,j]~tgrid, col="red",lty=2)  
# }
# 
# # Note: need to sample the observations
# 
# 
# corrupt<-function(X,times,sigma,gap=0.001){
#   Y<-X[times/gap+1,]
#   total_samples<-dim(Y)[1]*dim(Y)[2]
#   noise<-matrix(rnorm(total_samples, mean=0,sd=sigma), nrow=dim(Y)[1], ncol =dim(Y)[2] )
#   Y<-Y+noise
#   return(Y)
# }
# 
# tobs<-seq(from=0,to=1,length.out=n)
# observations<-corrupt(X=Xs,times=tobs,sigma=0.5,gap=0.001)
# 
# plot(Xs[,1]~tgrid,type="l",xlab="t",ylab="X",ylim=c(min(Xs),max(Xs)),col='red')
# for(j in 1:p){
#   lines(Xs[,j]~tgrid, col="red")  
# }
# 
# for(j in 1:p){
#   points(observations[,j]~tobs, col="red",pch=16)  
# }

###################
## Apply our methods on this data set:

grpNeighbours.Linear.integral<-function(Ys, Xhat,times,times_e,lambdas){
  # Evaluate the integrals at times
  dels=times_e[2]-times_e[1]
  n=length(times)
  p=dim(Ys)[2]
  Psi=matrix(0, nrow=n,ncol=p )
  for(i in 1:length(times)){
    Psi[i,] = apply(Xhat[ times_e<times[i],],2,sum)*dels  
  }
  response <- Ys
  intercepts<-rep(1,n)
  group<-c(1,1+{1:p})
  group_t<-group
  group_t[1]<-NA;cInd<-T
  
  covariate<-cbind(intercepts,Psi)
  scaleInd<-F;
  std_reg<-T
  out<-list()
  for(j in 1:p){
    
    out[[j]]<-grplasso(x=covariate, y=scale(response[,j],center=F,scale=scaleInd),
                       center=cInd,index=group_t,lambda=(lambdas*2*n), 
                       penscale=sqrt,model=LinReg(),standardize=std_reg)
    out[[j]]$beta<-out[[j]]$coef
  }  
  
  neighs<-array(0,dim=c(p,length(lambdas),p) )  
  beta.path<-list()
  for(j in 1:p){
    # Get the estimated neighbours
    neighs[,,j]<-Neighbours(est=out[[j]]$beta,grp=group,deg=dim(B[[1]])[2],p=p)
    beta.path[[j]]<-out[[j]]$beta
  }
  return(list(neighbour.path=neighs,coefficients.path=beta.path,covariate=covariate,group=group))
  
}

grpNeighbours.Linear.derivative<-function(Xdot, Xhat,times,lambdas){
  
  response <- Xdot
  intercepts<-rep(1,dim(response)[1])
  n<-dim(response)[1]
  group<-c(1,1+{1:p})
  group_t<-group
  group_t[1]<-NA;cInd<-T
  
  covariate<-Xhat
  covariate<-cbind(intercepts,covariate)
  scaleInd<-F;
  std_reg<-T
  out<-list()
  for(j in 1:p){
    
    out[[j]]<-grplasso(x=covariate, y=scale(response[,j],center=F,scale=scaleInd),
                       center=cInd,index=group_t,lambda=(lambdas*2*n), 
                       penscale=sqrt,model=LinReg(),standardize=std_reg)
    out[[j]]$beta<-out[[j]]$coef
  }  
  
  neighs<-array(0,dim=c(p,length(lambdas),p) )  
  beta.path<-list()
  for(j in 1:p){
    # Get the estimated neighbours
    neighs[,,j]<-Neighbours(est=out[[j]]$beta,grp=group,deg=dim(B[[1]])[2],p=p)
    beta.path[[j]]<-out[[j]]$beta
  }
  return(list(neighbour.path=neighs,coefficients.path=beta.path,covariate=covariate,group=group))
  
}