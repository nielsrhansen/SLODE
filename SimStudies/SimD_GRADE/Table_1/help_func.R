#
makeCI <- function(x,alpha=.95){
  xbar <- mean(x)
  s  <- sd(x)
  n <- length(x)
  z <- qnorm(alpha,0,1)
  return(c(xbar,xbar - z*s/sqrt(n),xbar + z*s/sqrt(n)))
}


# A function to reformat the GNW graph
# The first column is the regulator
# The second column is the regulatee
GNWgraph<-function(raw,gs){
  
  graph<-matrix(0,nrow=length(gs),ncol=length(gs))  
  for(i in 1:dim(raw)[1]){
    graph[gxTox(gx=raw[i,2],set=gs),
          gxTox(gx=raw[i,1],set=gs)]<-raw[i,3]
    if(raw[i,3]==0) break;
  }  
  return(graph)
}
# A mapping for GX and X
gxTox<-function(gx,set){ 
  return((1:length(set))[set==gx])
}

# A function to estimate the AUC 
# Adapted from Jamse Henderson's code. 

auc_us<-function(tps=NULL,nps=NULL,p=NULL,ntrue=NULL,approx=T,n.mc=1e3){
  ## Compute AUCs on ranked portion ##
  precision <- recall <- tpr <- fpr <- c()
  # Find out when the first edge comes in:
  istart<-min({1:length(nps)}[nps!=0])
  unranked<- ( max(nps) < {p^2-p}  )
  for(i in 1:{length(tps)-istart+1}){
    isub=istart+i-1;
    
    precision[i] <- tps[isub]/nps[isub]
    
    tpr[i]       <- tps[isub]/ ntrue
    fpr[i]       <- {nps[isub]-tps[isub]} / { p^2-ntrue-p  }
    
    if(i == 1){
      aucroc <- tpr[i]*(1-fpr[i])
      aucpr  <- precision[i]*tpr[i] #Extend horizontally to 0
    } else {
      aucroc <- aucroc + (tpr[i]-tpr[i-1]) * (1-fpr[i])	
      aucpr  <- aucpr +
        .5*(precision[i-1] + precision[i])*(tpr[i]-tpr[i-1])
    }
  }
  
  # Unranked edges 
  # Used when the edge sets are not full:
  ## Approximate expected AUC on unranked portion ##
  if(unranked==T){
    if(approx==T){
      out <- list()
      out$aucroc <- aucroc+(1-fpr[length(fpr)])*(1-tpr[length(fpr)])/2 
      out$aucpr <- aucpr
      out$precision <- precision
      out$tpr <- tpr
      out$fpr <- fpr
      return(out)
    } else{
      n.sat <- p^2-p
      max.edges=max(nps)
      n.u <- n.sat - max.edges
      ## Compute classification matrix for ranked edges ##
      n.pos <- ntrue 
      n.neg <- n.sat - n.pos
      tp.r <- max(tps) 
      fp.r <- max.edges - tp.r
      
      q <- n.pos - tp.r
      
      ## Monte Carlo Samples to Evaluate ##
      roc <- rep(NA,n.mc)
      pr  <- rep(NA,n.mc)
      
      for(k in 1:n.mc){
        tp.u <- tp.r
        fp.u <- fp.r
        
        precision.u <- c(precision[length(precision)])
        tpr.u <- c(tpr[length(tpr)])
        fpr.u <- c(fpr[length(fpr)])
        aucroc.u <- aucpr.u <- 0    
        
        ## Locations of true positives ##
        tp.loc <- sample(1:(n.sat-max.edges),q)
        
        ## Compute the AUCs for these locations ##
        for(i in 1:(n.sat-max.edges)){
          if(sum(i == tp.loc)==1){
            tp.u <- tp.u + 1
          } else{
            fp.u <- fp.u + 1
          }
          #           
          precision.u[i+1] <- tp.u / (tp.u + fp.u)
          tpr.u[i+1] <- tp.u / n.pos
          fpr.u[i+1] <- fp.u / n.neg
          #        
          aucroc.u <- aucroc.u + (tpr.u[i+1]-tpr.u[i]) * (1-fpr.u[i+1])
          aucpr.u  <- aucpr.u + .5*(precision.u[i] + precision.u[i+1])*
            (tpr.u[i+1]-tpr.u[i])
        }
        
        pr[k] <- aucpr.u
        roc[k] <- aucroc.u
      }
      out <- list()
      out$aucroc <- aucroc + mean(roc)
      out$aucpr <- aucpr + mean(pr)
      out$precision <- precision
      out$tpr <- tpr
      out$fpr <- fpr
      return(out)
    }
  } 
    
   
    
    
  out <- list()
  out$aucroc <- aucroc
  out$aucpr <- aucpr
  out$precision <- precision
  out$tpr <- tpr
  out$fpr <- fpr
  return(out)
}

# To verify the validity of this code
# We check with the NeRDS function
# p=50
# i=1
# load(file=paste("./RData/fig4",i,"NeRDS",p,".RData",sep=""))  
# AUC
# TPS<-(AUC$tpr*50)
# NPS<-((AUC$fpr*( p^2-p-50))+AUC$tpr*50)
# res<-auc_us(tps=TPS,nps=NPS,p=p,ntrue=50)
# res$aucpr
# AUC$aucpr
# res$aucroc
# AUC$aucroc


