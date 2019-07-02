# The main function

# Need to scale the observations before any regression!


GRADE<-function(observations, times_e,type_reg, lambdas, lambdas.fun=NULL,
                type_data,type_basis, type_smooth, xhat, xprime=NULL,...){
    
  p<-dim(observations[[1]])[2]-1
  if(type_basis=="penalized spline"){
   # Need to use the groupCD appraoch
    bases<-getPsi(xhat=xhat,times_e=times_e,type_basis=type_basis,  ...)
    lambda2s<-lambdas
    lambda1s<-lambdas.fun
      if(type_data=="Replicates"){
	  scores<-groupSearch.Y(observations=observations,bases=bases,times_e=times_e,
                            lambda1s=lambda1s,lambda2s=lambda2s)
      lambda1s.M<-array(0,dim=c(length(lambda2s),p )) 
      l1l2.M<-array(0,dim=c(2,p ))
      bicl2<-gcvl2<-array(0,dim=c(length(lambda2s),p ))
      for(j in 1:p){
        for(l in 1:length(lambda2s)){
          lambda1s.M[l,j]<-lambda1s[which.min(scores$gcv[,l,j])]
          bicl2[l,j]<-scores$bic[which.min(scores$gcv[,l,j]),l,j]
        }
      }
      l1l2.M[2,]<-lambda2s[which.min(apply(bicl2,1,sum))]
        fits<-neighbourhood.Y(observations=observations,bases=bases,
                            lambda1s.M=lambda1s.M,times_e=times_e,
                            lambda2s=lambda2s,l1l2.M=l1l2.M)
    
      } else if(type_data=="Perturbations"){
	  	  if(length(lambdas.fun)==1){
		  lambda1s.M<-array(lambdas.fun,dim=c(length(lambda2s),p )) 
		  l1l2.M<-NULL
	  } else{
		scores<-groupSearch.Y.P(observations=observations,bases=bases,times_e=times_e,
                                  lambda1s=lambda1s,lambda2s=lambda2s)
            lambda1s.M<-array(0,dim=c(length(lambda2s),p )) 
            l1l2.M<-array(0,dim=c(2,p ))
            bicl2<-gcvl2<-array(0,dim=c(length(lambda2s),p ))
            for(j in 1:p){
              for(l in 1:length(lambda2s)){
                lambda1s.M[l,j]<-lambda1s[which.min(scores$gcv[,l,j])]
                bicl2[l,j]<-scores$bic[which.min(scores$gcv[,l,j]),l,j]
              }
            }
            l1l2.M[2,]<-lambda2s[which.min(apply(bicl2,1,sum))]
	  }
             
            fits<-neighbourhood.Y.P(observations=observations,bases=bases,
                                  lambda1s.M=lambda1s.M,times_e=times_e,
                                  lambda2s=lambda2s,l1l2.M=l1l2.M)    
      }
    
   
 } else{
   
   # Prepare the basis
   bases<-getPsi(xhat=xhat,times_e=times_e,type_basis=type_basis,  ...)
     # check if the basis expansion is nonfullrank

	
   # Call the regression with the basis (keep things simple)
   if(type_reg=="Y"){
     # Call regression that comes with BICs
     if(type_data=="Replicates"){
       fits<-grpNeighbours.Y(observations=observations,bases=bases,
                             times_e=times_e,lambdas=lambdas,std_reg=T) 
     } else if(type_data=="Perturbations"){
       fits<-grpNeighbours.Y.P(observations=observations,bases=bases,
                             times_e=times_e,lambdas=lambdas,std_reg=T)        
     }
      
   }
	
     # Return the estimated neighbours 
   
    
 } 
 
# fits should be the return value of the intermediate function  
 return(fits)
  
}