

grpNeighbours.Y<-function(observations,bases,times_e,lambdas,std_reg,type_data="Replicates"){
  R<-length(observations)
  p<-dim(observations[[1]])[2]-1
  B<-list() 
  Bobs<-list()
    for(r in 1:R){
      bases$Psi[[r]]<-bases$Psi[[1]]
    }
  for(r in 1:R){
    indexes<-match(round(observations[[r]][,1],5), round(times_e,5) )
    Bobs[[r]]<-list()
    for(j in 1:p){
      Bobs[[r]][[j]]<-scale(bases$Psi[[r]][[j]][indexes,],center=T,scale=F)
    }
    
  }
  B<-Bobs[[1]]
  if(R > 1){
    for(r in 2:R){
      for(j in 1:p)
        B[[j]]<-rbind(B[[j]],Bobs[[r]][[j]])
    }
  }
  
  m<-dim(B[[1]])[2]
  n<-dim(B[[1]])[1]
  
  # create intercepts if data is from perturbation experiments
    response<-do.call(rbind,observations)[,-1]
    response<-scale(response,T,F)
    intercepts<-rep(1,dim(response)[1])
    times<-do.call(rbind,observations)[,1]
    group<-c(1,2,2+rep(1:p,each=m))
    group_t<-group
    group_t[c(1,2)]<-NA;cInd<-T
    
    covariate<-do.call(cbind,B)
    covariate<-cbind(intercepts,times, covariate)
    scaleInd<-T;
    
  # check if the basis expansion is nonfullrank
  minimunvalue<-1
  for(j in 1:p){
    minimunvalue<-min(minimunvalue,min(svd(B[[j]])$d))
  }
  if(minimunvalue<1e-5){
    gcvscore<-bicscore<-array(Inf,dim=c(length(lambdas),p ))  
    
    
    return(list(neighbour=NULL, neighbour.path=NULL,coefficients.path=NULL,BIC=bicscore,GCV=gcvscore, covariate=covariate,psi=bases$psi,group=group))
  } else{
    
  out<-list()
  gcvscore<-bicscore<-array(0,dim=c(length(lambdas),p ))  
  for(j in 1:p){
    
    out[[j]]<-grplasso(x=covariate, y=scale(response[,j],center=F,scale=scaleInd),
                       center=cInd,index=group_t,lambda=(lambdas*2*n), 
                       penscale=sqrt,model=LinReg(),standardize=std_reg)
    out[[j]]$beta<-out[[j]]$coef
    
    for(l in 1:length(lambdas)){
      betaest<-out[[j]]$beta[,l]
      MSE<-mean((scale(response[,j],F,scaleInd)-covariate%*%betaest)^2)
      nonezero<-round(abs(betaest),8)
      df<-sum(nonezero>0)
      bicscore[l,j]<- n*log(MSE)+log(n)*df
      gcvTemp<-MSE/(1-df/n)^2
      if(df<n){
        gcvscore[l,j]<-gcvTemp  
      } else{
        gcvscore[l,j]<- Inf
      }
      
    }
    
  }  
  
  neighs<-array(0,dim=c(p,length(lambdas),p) )  
  BICneighs<-array(0,dim=c(p,p) )  
  temp<-which.min( apply(bicscore,1,sum))
  beta.path<-list()
  for(j in 1:p){
    # Get the estimated neighbours
    neighs[,,j]<-Neighbours(est=out[[j]]$beta,grp=group,deg=dim(B[[1]])[2],p=p)
    BICneighs[,j]<-neighs[,temp,j]
	beta.path[[j]]<-out[[j]]$beta
  }
  
  return(list(neighbour=BICneighs, neighbour.path=neighs,coefficients.path=beta.path,BIC=bicscore,GCV=gcvscore, covariate=covariate,psi=bases$psi,group=group))
  }
  }


grpNeighbours.Y.P<-function(observations,bases,times_e,lambdas,std_reg){
  R<-length(observations)
  p<-dim(observations[[1]])[2]-1
  B<-list() 
  Bobs<-list()
  for(r in 1:R){
    indexes<-match(round(observations[[r]][,1],5), round(times_e,5) )
    Bobs[[r]]<-list()
    for(j in 1:p){
      Bobs[[r]][[j]]<-scale(bases$Psi[[r]][[j]][indexes,],center=T,scale=F)
    }
  }
  B<-Bobs[[1]]
  if(R > 1){
    for(r in 2:R){
      for(j in 1:p)
        B[[j]]<-rbind(B[[j]],Bobs[[r]][[j]])
    }
  }
  
  m<-dim(B[[1]])[2]
  n<-dim(B[[1]])[1]
  
    obs<-list()
  intercepts<-numeric(0)
    for(r in 1:R){
      obs[[r]]<-scale(observations[[r]],T,T)
      temp<-numeric(R)
      temp[r]<-1   
      intercepts<-rbind(intercepts,do.call(rbind,rep(list(temp),dim(obs[[r]])[1] )))
      
    }
    response<-do.call(rbind,obs)[,-1]
    
    times<-do.call(rbind,observations)[,1]
    times<-scale(times,T,F)
    group<-c(1:R,R+1,R+1+rep(1:p,each=m))
    group_t<-group
    group_t[1:(R+1)]<-NA;
    cInd<-T
    
    covariate<-do.call(cbind,B)
    covariate<-cbind(intercepts,times, covariate)
    scaleInd<-F;
  
  out<-list()
    gcvscore<-bicscore<-array(0,dim=c(length(lambdas),p ))  
  for(j in 1:p){
    out[[j]]<-grplasso(x=covariate, y=scale(response[,j],center=T,scale=scaleInd),
                       center=cInd,index=group_t,lambda=(lambdas*2*n), 
                       penscale=sqrt,model=LinReg(),standardize=std_reg)
    out[[j]]$beta<-out[[j]]$coef
    for(l in 1:length(lambdas)){
      betaest<-out[[j]]$beta[,l]
      MSE<-mean((scale(response[,j],T,scaleInd)-covariate%*%betaest)^2)
      nonezero<-round(abs(betaest),8)
      df<-sum(nonezero>0)
      bicscore[l,j]<- n*log(MSE)+log(n)*df
      gcvTemp<-MSE/(1-df/n)^2
      if(df<n){
        gcvscore[l,j]<-gcvTemp  
      } else{
        gcvscore[l,j]<- Inf
      }
    }
  }  
  neighs<-array(0,dim=c(p,length(lambdas),p) )  
  BICneighs<-array(0,dim=c(p,p) )  
  temp<-which.min( apply(bicscore,1,sum))
  beta.path<-list()  
  for(j in 1:p){
    # Get the estimated neighbours
    neighs[,,j]<-Neighbours(est=out[[j]]$beta,grp=group,deg=dim(B[[1]])[2],p=p)
    BICneighs[,j]<-neighs[,temp,j]
	beta.path[[j]]<-out[[j]]$beta
	
  }
  return(list(neighbour.path=neighs,neighbour=BICneighs,coefficients.path=beta.path,BIC=bicscore,GCV=gcvscore,covariate=covariate,psi=bases$psi,group=group)) 
}
# Need edits
grpNeighbours.Xprime.P<-function(observations,bases,times_e,lambdas,std_reg){
  R<-length(observations)
  p<-dim(observations[[1]])[2]-1
  B<-list() 
  Bobs<-list()
  for(r in 1:R){
    indexes<-match(round(observations[[r]][,1],5), round(times_e,5) )
    Bobs[[r]]<-list()
    for(j in 1:p){
      Bobs[[r]][[j]]<-scale(bases$psi[[r]][[j]][indexes,],center=T,scale=F)
    }
  }
  B<-Bobs[[1]]
  if(R > 1){
    for(r in 2:R){
      for(j in 1:p)
        B[[j]]<-rbind(B[[j]],Bobs[[r]][[j]])
    }
  }
  
  m<-dim(B[[1]])[2]
  n<-dim(B[[1]])[1]
  
    obs<-list()
  intercepts<-numeric(0)
    for(r in 1:R){
      obs[[r]]<-scale(observations[[r]],T,T)
      temp<-numeric(R)
      temp[r]<-1   
      intercepts<-rbind(intercepts,do.call(rbind,rep(list(temp),dim(obs[[r]])[1] )))
      
    }
    response<-do.call(rbind,obs)[,-1]
    
    times<-do.call(rbind,observations)[,1]
    times<-scale(times,T,F)
    group<-c(1,1+rep(1:p,each=m))
    group_t<-group
    group_t[1:(R+1)]<-NA;
    cInd<-T
    
    covariate<-do.call(cbind,B)
    covariate<-cbind(intercepts,times, covariate)
    scaleInd<-F;
  
  out<-list()
    gcvscore<-bicscore<-array(0,dim=c(length(lambdas),p ))  
  for(j in 1:p){
    out[[j]]<-grplasso(x=covariate, y=scale(response[,j],center=T,scale=scaleInd),
                       center=cInd,index=group_t,lambda=(lambdas*2*n), 
                       penscale=sqrt,model=LinReg(),standardize=std_reg)
    out[[j]]$beta<-out[[j]]$coef
    for(l in 1:length(lambdas)){
      betaest<-out[[j]]$beta[,l]
      MSE<-mean((scale(response[,j],T,scaleInd)-covariate%*%betaest)^2)
      nonezero<-round(abs(betaest),8)
      df<-sum(nonezero>0)
      bicscore[l,j]<- n*log(MSE)+log(n)*df
      gcvTemp<-MSE/(1-df/n)^2
      if(df<n){
        gcvscore[l,j]<-gcvTemp  
      } else{
        gcvscore[l,j]<- Inf
      }
    }
  }  
  neighs<-array(0,dim=c(p,length(lambdas),p) )  
  BICneighs<-array(0,dim=c(p,p) )  
  temp<-which.min( apply(bicscore,1,sum))
  beta.path<-list()  
  for(j in 1:p){
    # Get the estimated neighbours
    neighs[,,j]<-Neighbours(est=out[[j]]$beta,grp=group,deg=dim(B[[1]])[2],p=p)
    BICneighs[,j]<-neighs[,temp,j]
	beta.path[[j]]<-out[[j]]$beta
	
  }
  return(list(neighbour.path=neighs,neighbour=BICneighs,coefficients.path=beta.path,BIC=bicscore,GCV=gcvscore,covariate=covariate,psi=bases$psi,group=group)) 
}
