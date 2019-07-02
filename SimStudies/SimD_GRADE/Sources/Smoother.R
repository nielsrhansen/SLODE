
## Smoothing: smooth the data

smoothX<-function(observations,times_e,deg,h=NULL,maxk=5000,
                       type_smooth=c("local polynomial","smoothspline" ,"linear interpolation"),
                       type_data="Replicates"){
  Xprime<-Xhat<-list()
  gcvscore<-list()
  if(type_data=="Replicates"){    
      obs<-do.call(rbind,observations)
      p<-dim(obs)[2]-1

      Xprime[[1]]<-Xhat[[1]]<-matrix(0,nrow=length(times_e),ncol= p)
      
      if( type_smooth=="linear interpolation") {
        for(i in 1:p){Xhat[[1]][,i]<-approx(x=obs[,1],y=obs[,i+1] ,xout=times_e,rule=2)$y
          gcvscore[[1]]<-NULL
          Xprime[[1]]<-NULL
          deg<-NULL
        }
      } else if(type_smooth=="local polynomial"){
        gcvscore[[1]]<-matrix(0,nrow=length(h),ncol= p)
        
        # Get the LP estimates 
        for(i in 1:p){
          # For the lp estimates, we pick the bandwidth that minimize the generalized cross-validation
          # for each node.
          for(j in 1:length(h)){
            h_temp<-h[j]
            temp <- locfit(obs[,i+1]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)
            gcvscore[[1]][j,i]<-gcv(obs[,i+1]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)[4]
          }
          h_temp<-h[which.min(gcvscore[[1]][,i])]
          temp <- locfit(obs[,i]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)
          Xhat[[1]][,i]<-predict(temp,newdata=times_e)
          Xprime[[1]][,i]<-  NULL
        }
        
      } else if (type_smooth == "smoothspline"){
        for(i in 1:p){
          gcvscore[[1]]<-numeric(p)
          
          temp <- smooth.spline(obs[,1], obs[,i+1], all.knots = T)
          gcvscore[[1]][i]<-temp$cv.crit
          Xhat[[1]][,i]<-predict(temp,times_e)$y
          Xprime[[1]][,i]<-  predict(temp,times_e,deriv=1)$y
        }
        
      }   
  } else if (type_data=="Perturbations"){
    # Similar to the previous function, but create a list for every experiment
    R<-length(observations)
    for(r in 1:R){
      obs<-observations[[r]]
      p<-dim(obs)[2]-1
      
      Xprime[[r]]<-Xhat[[r]]<-matrix(0,nrow=length(times_e),ncol= p)
      
      if( type_smooth=="linear interpolation") {
        for(i in 1:p){
          Xhat[[r]][,i]<-approx(x=obs[,1],y=obs[,i+1] ,xout=times_e,rule=2)$y
          gcvscore[[r]]<-NULL
          Xprime[[r]]<-NULL
          deg<-NULL
        }
      } else if(type_smooth=="local polynomial"){
        # Get the LP estimates 
        for(i in 1:p){
          gcvscore[[r]]<-matrix(0,nrow=length(h),ncol= p)
          for(j in 1:length(h)){
            h_temp<-h[j]
            temp <- locfit(obs[,i+1]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)
            gcvscore[[r]][j,i]<-gcv(obs[,i+1]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)[4]
          }
          h_temp<-h[which.min(gcvscore[[r]][,i])]
          temp <- locfit(obs[,i]~lp(obs[,1],deg=deg,h= h_temp),maxk=maxk)
          Xhat[[r]][,i]<-predict(temp,newdata=times_e)
          Xprime[[r]][,i]<-  NULL
        }
        
      } else if (type_smooth == "smoothspline"){
        for(i in 1:p){
          gcvscore[[r]]<-numeric(p)
          
          temp <- smooth.spline(obs[,1], obs[,i+1], all.knots = T)
          gcvscore[[r]][i]<-temp$cv.crit
          Xhat[[r]][,i]<-predict(temp,times_e)$y
          Xprime[[r]][,i]<-  predict(temp,times_e,deriv=1)$y
        }
        
      }
      
    
    }
    
    }
  
  return(list(Xhat=Xhat,Xprime=Xprime,times_e=times_e,gcvscore=gcvscore,type_smooth=type_smooth))
  
}
