##---------------------------------------------------------------##
#### Graph generation #####
# Generate a simple graph 
p<-8   # number of variables 

set.seed(1)
ks<-seq(from=1,to=p/2,length.out=p/2)
ks<-ks*pi*2
thetas<-matrix(0,p,p)
for(j in 1:{p/2}){
  thetas[ 2*j-1  ,2*j]<- ks[j]
  thetas[2*j, 2*j-1]<- {-ks[j]}
}
intercepts<-c(rep(0,p) )
thetas<-cbind(intercepts,thetas)
tgraph<-t(thetas[,-1]!=0)

