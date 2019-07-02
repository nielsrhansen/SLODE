
######################
# B are the penalized blocks, and Bup are unpenalized blocks (not including the intercept)

groupCD<-function (Y, B, LHS,Bup,lambda2, gamma0.ini=NULL,gammas.ini=NULL,
                      thresh = 1e-03, max.iter = 50,sLHS=NULL,LHSinvB.array=NULL){
  p<-length(B)
  if(is.null(sLHS)){
  sLHS<-list()
	for(j in 1:p){
    sLHS[[j]]<-solve(LHS[[j]])
  }
  } else{
	# Do nothing if sLHS is feeded
	
  }
  
  resid <- as.matrix(Y- mean(Y))
  crit <- 1
  crits <- c()
  iter <- 0
  sup <- rep(0, p)
  
  fs <- matrix(0, nrow = dim(resid)[1], ncol = p)
  f0<-numeric(dim(resid)[1])
  
  m <- dim(B[[1]])[2]
  # Use Initial values:
  if(is.null(gammas.ini)){
    gammas <- list()
  for (j in 1:p) {
    gammas[[j]] <-numeric(m)
  }
  } else{
  gammas=gammas.ini
  }
if(is.null(gamma0.ini) ){
	  gamma0<-0
} else{ 
	gamma0=gamma0.ini
}
  
  mus <- mean(Y)
  Yresid <- Y-mean(Y)
  if(is.null(Bup)){
    sup0<-0
  } else{
    Bup<-scale(Bup,center=T,scale=F)
    sBup<-solve(t(Bup)%*%Bup)
  } 
	
	# Put parameters into arrays
	B.array<-array(0,dim=c(p,dim(B[[1]])))
	for(j in 1:p){
		B.array[j,,]<-B[[j]]
	}
	sLHS.array<-array(0,dim=c(p,dim(sLHS[[1]])))
	for(j in 1:p){
		sLHS.array[j,,]<-sLHS[[j]]
	}
	gammas.array<-array(0,dim=c(p,length(gammas[[1]])))
	for(j in 1:p){
		gammas.array[j,]<-gammas[[j]]
	}
	if( is.null(LHSinvB.array)){
	LHSinvBmat<-sLHS.array[1,,]%*%t(B.array[1,,])
	LHSinvB.array<-array(0,dim=c(p,dim(LHSinvBmat)))
	for(j in 1:p){
		LHSinvB.array[j,,]<-sLHS.array[j,,]%*%t(B.array[j,,])
	}
	} else{
	# Do nothing since it's provided
	} 
	
	while (iter < max.iter & crit > thresh) {
    if(is.null(Bup)){
      gamma0<-NULL
    } else {
      resid<-Yresid-apply(as.matrix(fs), 1, sum)
      gamma0 <- sBup%*%t(Bup)%*%resid
      f<-f0
      f0<-Bup%*%gamma0
      sup0<- max(abs(f-f0))
    }
    for (j in 1:p) {
      resid <- Yresid - apply(as.matrix(fs[, -j]), 1, sum)-f0
      f <- fs[, j]
      gammas.array[j,] <- LHSinvB.array[j,,]%*%resid
      fs[, j] <- B.array[j,,] %*% gammas.array[j,]
      si <- sqrt(sum(fs[, j]^2))
      soft.thresh <- max(1 - lambda2/si, 0)
      gammas.array[j,] <- soft.thresh * gammas.array[j,]
      fs[, j] <- soft.thresh * fs[, j]
      fs[, j] <- fs[,j]-mean(fs[,j])
      sup[j] <- max(abs(fs[, j] - f))
    }
    crit <- max(c(sup,sup0))
    iter <- iter + 1
    crits <- c(crits, crit)
  }
  fit <- list()

  fit$gammas <- list()
  for(j in 1:p){
  fit$gammas[[j]]<-gammas.array[j,]
  }
  
  fit$f.hat <- fs
  cnv <- ifelse(crit < thresh, 0, 1)
  if (cnv == 1) 
    warning(cat("Model did not converge \n"))
  fit$means <- mus
  fit$gamma0 <- gamma0
  fit$cnv <- cnv
  return(fit)
}



