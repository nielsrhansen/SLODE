
###############
# Start of function: getSolution
# Input: theta1, theta2, theta3
#        Nnoise: number of noise solutions
#        times_fine 
# Output: a list: mean_curve, a matrix of mean curves
#                 times_fine, a vector of times points that are measured
#######################


getSolution<-function(theta1,theta2,theta3,iv1,iv2,iv3,Nnoise,gap=0.001,range=c(0,20)){
  coef<-rnorm(Nnoise)
  initial<-rnorm(Nnoise)  
  times_fine<-seq(from=range[1],to=range[2],by=0.001) # The time points used to solve the ODE systems
  noise<-matrix(0,length(times_fine),Nnoise)
  for(i in 1:Nnoise){
    noise[,i]<-simple_ode(time_grid=times_fine,c=coef[i],iv=initial[i])
  }  
  Out1<-EulerM(f=FitzHugh_Nagumo_extend, theta=theta1,Initial_value=iv1,gap=gap)  
  Out2<-EulerM(f=FitzHugh_Nagumo_extend, theta=theta2,Initial_value=iv2,gap=gap) 
  Out3<-EulerM(f=FitzHugh_Nagumo_extend, theta=theta3,Initial_value=iv3,gap=gap) 
  
  mean_curve<-cbind(Out1,Out2,Out3)
#   for(i in 1:6){
#     mean_curve[,i]<-scaler(mean_curve[,i],range=4)
#   }
# #   
  mean_curve<-cbind(mean_curve,noise)
  # Create a matrix to store the graph information
  # the second index indices the respond variable
  Nvar<-6+Nnoise
  graph<-matrix(0,nrow=Nvar,ncol=Nvar)
  
  #several elements are determined
  graph[c(1:2),1]<-1
  graph[c(1:2),2]<-1
  graph[c(4),3]<-1
  graph[3,4]<-1
  graph[6,5]<-1
  graph[5,6]<-1
  
  theta<-matrix(0,nrow=Nvar,ncol=3*Nvar+1)
  theta[c(1,2),1:7]<-theta1
  theta[c(3,4),c(1,8:13)]<-theta2
  theta[c(5,6),c(1,14:19)]<-theta3
  theta[7,1]<-coef[1]
  theta[8,1]<-coef[2]
  theta[9,1]<-coef[3]
  theta[10,1]<-coef[4]
  mc<-list()
  mc[[1]]<-mean_curve
  return(list(mean_curve=mc, times_fine=times_fine, gap=gap,graph=graph,theta=theta))  
}




##########
# Start of function: observations
# Input: 
#        times, the time points where the values are observed.
#        solutions, the true solutions (a list from getSolution)
# Output: 
#        observations
#        times
##########################################

# R is the number of systems(replicates)

getObservations<-function(times,solutions,sigma=1,R=1,type_data){
 out<-list()
 if(type_data=="Replicates"){
   for(r in 1:R){
     out[[r]]<-cbind(times,corrupt(X=solutions$mean_curve[[1]],times=times,sigma=sigma,gap=solutions$gap))
   }
   
 } else if( type_data=="Perturbations") {
   for(r in 1:R){
     out[[r]]<-cbind(times,corrupt(X=solutions$mean_curve[[r]],times=times,sigma=sigma,gap=solutions$gap))
   }
 }
    return(out)    
}


wrapper_bandwidth<-function(...){
  res_hs<-sapply(hs,wrapper, simplify=F, USE.NAMES=F,...)
  
  gcvscore<-numeric(length(hs))
  bic<-rss<-matrix(0,nrow=length(hs),ncol=lambda_N)
  out<-array(0,dim=c(dim(res_hs[[1]]$out ), length(hs)) ) 
  for(l in 1:length(hs)){
    gcvscore[l]<-res_hs[[l]]$gcvscore
    bic[l,]<-res_hs[[l]]$bic
    rss[l,]<-res_hs[[l]]$rss
    out[,,l]<-res_hs[[l]]$out  
  }
  return(list(ROC=out,bic=bic,rss=rss,gcvscore=gcvscore) )  
}

wrapper<-function(h,...){
  
  if(type_reg=="Derivative"){
    smoothed<-smooth_curve(observations=observations,times_e=times,deg=deg,h=h,
                           maxk=maxk,type=type_smooth,type_data=type_data,R=R)  
    
  } else{
    smoothed<-smooth_curve(observations=observations,times_e=times_e,deg=deg,h=h,
                           maxk=maxk,type=type_smooth,type_data=type_data,R=R)  
    
  }
  
  Psihat<-getPsi(smoothed=smoothed,observations=observations,
                 type=type_reg,fun_basis=fun_basis,std=std,...)
  
  est<-grp_regression(smoothed=smoothed, Psi=Psihat,upp=lambda_range[2],low=lambda_range[1],
                      Nlambda=lambda_N,stdus=stdus, std_reg=std_reg,observations=observations,
                      type=type_reg,type_data=type_data,R=R,type_pack=type_pack)    
  
  graph<-getGraph(estimates=est,type=type_reg)
  
  counts<-countGraph(graph_est=graph,graph_true=solutions$graph,self=self)  
  
  criteria<-getCriteria(neighbours=graph,observations=observations,type=type_reg,
                        type_data=type_data, smoothed=smoothed,Psi=Psihat,estimates=est,R=R)
  
  gcvscore<-sum(smoothed$gcvscore)
  out<-cbind(counts$NP_sum,counts$TP_sum)
  bic<-apply(criteria$bic,MARGIN=2,sum)
  rss<-apply(criteria$rss,MARGIN=2,sum)
  print(paste("h =", h, "done."))
  return(list(out=out,criteria=criteria,rss=rss,bic=bic,gcvscore=gcvscore))
}

