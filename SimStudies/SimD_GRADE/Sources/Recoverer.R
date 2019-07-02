
getGraph<-function(estimates,type="Y"){
  #observations is useless here)
  
  R_sm<-  estimates$R_sm
  deg_basis<-estimates$deg_basis
  if(type=="Derivative"){
    p<-(max(estimates$fits[[1]]$group)-R_sm)
  } else{p<-(max(estimates$fits[[1]]$group)-R_sm-1) }
  neighs<-array(0,dim=c(p,length(estimates$fits[[1]]$lambda),p) )  
  for(j in 1:p){
    fits<-estimates$fits[[j]]
    # Get the estimated neighbours
    neighs[,,j]<-Neighbours(est=fits$beta,grp=fits$group,deg=deg_basis,p=p)
  }
  return(neighbours=neighs)
}

### Neigh_Recover
## Turn the estimated coefficients to a neighbourhood (regulators) of a node
##
## 
CoeftoNeigh<-function(theta,group,deg,p,pre=0.001){
  neigh<-numeric(p)
  # Ignore the first two components of vectors theta and group
  NumInt<-(max(group)-p)
  for(i in (NumInt+1):length(group)){
    neigh[group[i]-NumInt]<-neigh[group[i]-NumInt]+theta[i]^2
  }
  return((sqrt(neigh)>pre))
}
Neighbours<-function(est,grp,pre=0.001,deg,p){
  neigh<-matrix(0,ncol=dim(est)[2],nrow=p)
  for(i in 1:dim(est)[2]){
    neigh[,i]<-CoeftoNeigh(theta=est[,i],group=grp,pre=pre,deg=deg,p=p)    
  }  
    return(neigh)
}

Betas<-function(est,grp,pre=0.001,deg,p,m){
  betas<-array(0,c(p,dim(est)[2],m ))
  for(i in 1:dim(est)[2]){
    betas[,i,]<-CoeftoBeta(theta=est[,i],group=grp,pre=pre,deg=deg,p=p,m=m)    
  }  
    return(betas)
}

CoeftoBeta<-function(theta,group,deg,p,pre=0.001,m){
  beta<-array(0,c(p,m))
  # Ignore the first two components of vectors theta and group
  NumInt<-(max(group)-p)
  grpcount<-rep(1,p)
  for(i in (NumInt+1):length(group)){
    itemp<-group[i]-NumInt
	beta[itemp,grpcount[itemp]]<- beta[itemp]+theta[i]^2
	grpcount[itemp]<-grpcount[itemp]+1
  }
  return(beta)
}
countGraph<-function(graph_est,graph_true,self=T){
  dims<-dim(graph_est)
  NP<-TP<-matrix(0,nrow=dims[1],ncol=dims[2])
  if(self==T){
    if(length(dims)==3) {
    for(j in 1:dims[3]){
      for(i in 1:dims[2]){
        NP[j,i]<-sum(graph_est[,i,j]) #NP
        TP[j,i]<-sum( (graph_est[,i,j]==graph_true[,j])&(graph_true[,j]==1) ) #TP      
      }
    }
  } else{
    for(j in 1:dims[2]){
      NP[j,1]<-sum(graph_est[,j]) #NP
      TP[j,1]<-sum( (graph_est[,j]==graph_true[,j])&(graph_true[,j]==1) ) #TP      
    }
    
  }
    
  } else{
    
    if(length(dims)==3) {
    for(j in 1:dims[3]){
      for(i in 1:dims[2]){
        NP[j,i]<-sum(graph_est[-j,i,j]) #NP minus diagonal
        TP[j,i]<-sum( (graph_est[-j,i,j]==graph_true[-j,j])&(graph_true[-j,j]==1) ) #TP  minus diagonal    
      }
    }
  } else{
    for(j in 1:dims[2]){
      NP[j,1]<-sum(graph_est[-j,j]) #NP
      TP[j,1]<-sum( (graph_est[-j,j]==graph_true[-j,j])&(graph_true[-j,j]==1) ) #TP      
    }
    
  }
    
  }
  
  NP_sum<-apply(NP,2,sum)
  TP_sum<-apply(TP,2,sum)  
  return(list(NP=NP,TP=TP,NP_sum=NP_sum,TP_sum=TP_sum))
}
