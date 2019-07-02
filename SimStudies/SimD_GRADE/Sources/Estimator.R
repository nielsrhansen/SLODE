
###########
# Start of function: getPsi
# Output should be a list of Basis matrices and penalty matrices (if any)
# For each system, create one set of estimates

getPsi<-function(xhat,times_e,type_basis, std=T, ...){
  
  R<-length(xhat)  # For replicates, R_sm=1; For Pertubations, R_sm>1; 
  Psi<-psi<-list()
  K<-list()
  if(type_basis=="penalized spline"){
    bsb<-list()
    p<-dim(xhat[[1]])[2]
    N<-length(times_e)
    xall<-do.call(rbind,xhat)
    for(j in 1:p){
      x_temp<-xall[,j]
      specs<-getPercentile(x=x_temp,...)
      brs<-specs$percentile    
      bsb[[j]]<-create.bspline.basis(norder=specs$norder,breaks=brs )
      temp2<-getbasispenalty(bsb[[j]])
      K[[j]]<-temp2
    }
    for(r in 1:R){    
      psi[[r]]<-list()
      Psi[[r]]<-list()
      for(j in 1:p){
        temp<-getbasismatrix(xhat[[r]][,j], bsb[[j]])
        deg_basis<-dim(temp)[2]
        psi[[r]][[j]]<-temp
        tempPsi<-temp
        for(k in 1:deg_basis){
          tempPsi[,k]<-cumsum(temp[,k])*max(times_e)/N
        }
        Psi[[r]][[j]]<-tempPsi
      }
    }
  } else {
    
    
    for(r in 1:R){    
      p<-dim(xhat[[r]])[2]
      N<-length(times_e)
      psi[[r]]<-list()
      Psi[[r]]<-list()
      
      for(j in 1:p){
        temp2<-numeric(0)
        if(type_basis== "bspline"){
          specs<-getPercentile(x=xhat[[r]][,j],...)
          knts<-specs$percentile[-c(1,length(specs$percentile))]
          temp<-bs(xhat[[r]][,j],degree=specs$norder,knots=knts )
        } else  if(type_basis=="SAODEspline") {
          
          specs<-getPercentile(x=xhat[[r]][,j],...)
          brs<-c(rep(min(xhat[[r]][,j]),3), specs$percentile, rep(max(xhat[[r]][,j]),3))
          bsb = SplineBasis(brs)
          temp<-evaluate(bsb,xhat[[r]][,j])
          
        }  else if(type_basis=="fourier"){
          rangevalue<-range(xhat[[r]][,j])+c(-0.1,0.1)
          bsb <- create.fourier.basis(rangeval=rangevalue, ...)
          temp<-eval.basis(xhat[[r]][,j],bsb)
          
        } else if(type_basis=="monomial") {
          rangevalue<-range(xhat[[r]][,j])+c(-0.1,0.1)
          bsb <- create.monomial.basis(rangeval=rangevalue, ...)
          temp<-eval.basis(xhat[[r]][,j],bsb)
        }
        deg_basis<-dim(temp)[2]
        psi[[r]][[j]]<-temp
        tempPsi<-temp
        for(k in 1:deg_basis){
          tempPsi[,k]<-cumsum(temp[,k])*max(times_e)/N
        }
        Psi[[r]][[j]]<-tempPsi
      }
    }   
    
  }
  return(list(Psi=Psi,psi=psi,deg_basis=deg_basis,K=K))
}

getPercentile<-function(x,nspecs){
  nknots<-nspecs[2]
  norder<-nspecs[1]
  if(nknots>0){
  out<-seq(from=0, to=1, length.out=(nknots+2) )
  out<-quantile(x,probs=out )
  } else{out=NULL }
  return(list(percentile=out,norder=norder))
}


# use the Gram-Schimit process to orthogonalized the basis 

#       
#       if(std==T){
#         temp<-scale(temp,center=T,scale=F)
#         temp[,1]<-temp[,1]/sqrt(sum(temp[,1]^2))
#         if( deg_basis>1){
#           for(k in 2:deg_basis){
#             temp[,k]<-residuals(lsfit(x=temp[,1:(k-1)],y=temp[,k],intercept=F))
#             temp[,k]<-temp[,k]/sqrt(sum(temp[,k]^2))
#           }    
#         }
#       }


### Extract values at measured time points
#     if(type=="Y"){
#       d<-1    
#       times<-observations[[r]][,1]
#       covariate[[r]]<-matrix(0,nrow=length(times),ncol=dim(Psi[[r]])[2])
#       for(k in 1:length(smoothed$times_e)){
#         if(abs(times[d]-smoothed$times_e[k])<1e-5){
#           covariate[[r]][d,]<-c(Psi[[r]][k,])
#           d<-d+1
#         }  
#         if(d > length(times)) break;          
#       }
#     }
