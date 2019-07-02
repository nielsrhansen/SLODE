

neighbourhood.Y<-function(observations,bases,times_e,lambda1s.M,lambda2s,l1l2.M){
  R<-length(observations)
  p<-dim(observations[[1]])[2]-1
  B<-list()
  K<-list()
  K<-bases$K
    for(r in 1:R){
      bases$Psi[[r]]<-bases$Psi[[1]]
    }
  # matching is done for each replicate in case the observe time is different, but the smooth estimates are the same
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
    for(r in 1:R){
      for(j in 1:p)
        B[[j]]<-rbind(B[[j]],Bobs[[r]][[j]])
    }
  }
  
  m<-dim(B[[1]])[2]
  n<-dim(B[[1]])[1]
  
  response<-do.call(rbind,observations)[,-1]
  
  intercepts<-rep(1,dim(response)[1])
  times<-do.call(rbind,observations)[,1]
  
  # using the given matrix of lambda1  
  CDcoef<-array(0,dim=c(p*dim(B[[1]])[2]+2,length(lambda2s),p ))
  for(j in 1:p){
    for(l in 1:length(lambda2s)){
      L2<-lambda2s[l]
      L1<- lambda1s.M[l,j]
      LHSpen<-list()
      for(tt in 1:p){
        LHSpen[[tt]]<-t(B[[tt]])%*%B[[tt]]+L1*K[[tt]]
      }
      
      CDfits <- groupCD(Y=scale(response[,j],T,T), B, LHSpen, Bup=times,lambda2=L2,
                        thresh = 1e-05, max.iter = 1000) 
      CDcoef[,l,j]<- c(CDfits$means,CDfits$gamma0,do.call(c,CDfits$gammas))
    }
  }
  
  group<-c(1,1,1+1+rep(1:p,each=dim(B[[1]])[2]))
  
  neighs<-array(0,dim=c(p,length(lambda2s),p) )  
  for(j in 1:p){
    # Get the estimated neighbours
    neighs[,,j]<-Neighbours(est=CDcoef[,,j],grp=group,deg=dim(B[[1]])[2],p=p)
  }
BICneighs<-array(0,dim=c(p,p) )  
for(j in 1:p){
  temp<-which(lambda2s==l1l2.M[2,j])
  # Get the estimated neighbours
  BICneighs[,j]<-neighs[,temp,j]
}

return(list(neighbour.path=neighs,neighbour=BICneighs,coefficients.path=CDcoef, psi=bases$psi))

}




neighbourhood.Y.P<-function(observations,bases,times_e,lambda1s.M,lambda2s,l1l2.M){
  R<-length(observations)
  p<-dim(observations[[1]])[2]-1
  B<-list()
  K<-list()
  K<-bases$K
  
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
  for(r in 1:R){
    obs[[r]]<-scale(observations[[r]],T,F)
  }
  response<-do.call(rbind,obs)[,-1]
  
  # create intercepts if data is from perturbation experiments
  intercepts<-rep(1,dim(response)[1])
  times<-do.call(rbind,observations)[,1]
  if(sd(lambda1s.M)==0){
	LHSpen<-list()
	sLHS<-list()
	for(tt in 1:p){
			LHSpen[[tt]]<-t(B[[tt]])%*%B[[tt]]+lambda1s.M[1,1]*K[[tt]]
			sLHS[[tt]]<-solve(LHSpen[[tt]])
	}
	LHSinvBmat<-sLHS[[1]]%*%t(B[[1]])
	LHSinvB.array<-array(0,dim=c(p,dim(LHSinvBmat)))
	for(tt in 1:p){
		LHSinvB.array[tt,,]<-sLHS[[tt]]%*%t(B[[tt]])
	}
	} else{
	# do nothing
	}
  # using the given matrix of lambda1  
  CDcoef<-array(0,dim=c(p*dim(B[[1]])[2]+1+1,length(lambda2s),p ))
  for(j in 1:p){
    for(l in 1:length(lambda2s)){
      L2<-lambda2s[l]
      L1<- lambda1s.M[l,j]
	  if(sd(lambda1s.M)!=0){
		  LHSpen<-list()
		  for(tt in 1:p){
			LHSpen[[tt]]<-t(B[[tt]])%*%B[[tt]]+L1*K[[tt]]
		  }
		  sLHS<-NULL
		  LHSinvB.array<-NULL
	  } else{
	   # do nothing
	   
	  }
      if(l==1){
	  CDfits <- groupCD(Y=scale(response[,j],T,T), B, LHSpen, Bup=times,lambda2=L2,
                        thresh = 1e-03, max.iter = 2e2,sLHS=sLHS,LHSinvB.array=LHSinvB.array) 
		gamma0.pre<- CDfits$gamma0;
		gammas.pre<-CDfits$gammas
	  } else{
	  CDfits <- groupCD(Y=scale(response[,j],T,T), B, LHSpen, Bup=times,lambda2=L2,
                        thresh = 1e-03, max.iter = 2e2,gamma0.ini=gamma0.pre,gammas.ini=gammas.pre,sLHS=sLHS,LHSinvB.array=LHSinvB.array) 
		gamma0.pre<- CDfits$gamma0;
		gammas.pre<-CDfits$gammas
	  }
	  if(sum(do.call(c,CDfits$gammas)!=0)=={m*p}){
		print(c(j, "finished"))
		for(li in l:length(lambda2s)){
		CDcoef[,li,j]<- c(CDfits$means,CDfits$gamma0,do.call(c,CDfits$gammas))
		}
		break;
	  } else{
	  print(c(j, "lambdas", l, sum(do.call(c,CDfits$gammas)!=0)))
      CDcoef[,l,j]<- c(CDfits$means,CDfits$gamma0,do.call(c,CDfits$gammas))
	  }
      
    }
  }
  
  
  group<-c(1,1,2+rep(1:p,each=dim(B[[1]])[2]))
  
  neighs<-array(0,dim=c(p,length(lambda2s),p) )  
  for(j in 1:p){
    # Get the estimated neighbours
    neighs[,,j]<-Neighbours(est=CDcoef[,,j],grp=group,deg=dim(B[[1]])[2],p=p)
  }
  
  #BICneighs<-array(0,dim=c(p,p) )  
  #for(j in 1:p){
  #  temp<-which(lambda2s==l1l2.M[2,j])
    # Get the estimated neighbours
  #  BICneighs[,j]<-neighs[,temp,j]
  #}
  
  return(list(neighbour.path=neighs,coefficients.path=CDcoef,psi=bases$psi))
  
}

