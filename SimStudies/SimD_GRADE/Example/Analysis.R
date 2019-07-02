  ##------------------------------------------------------------------------##
  #### Method  ####
  
  type_reg="Y"; # use Y as the outcome 
  type_smooth = "smoothspline" # use smooth splines to fit the smooth trajectories 
  
  # Find the smooth trajectories
  smthed<-smoothX(observations=observations,times_e=times_e,deg=deg,h=h,type_smooth=type_smooth,type_data=type_data)

  
  xhat<-smthed$Xhat

  # Set the tuning parameter 
  lambda_N<-30 
  lambda_range<-c(-4,0)
  lambdas<-exp(seq(from=lambda_range[2],to=lambda_range[1],length.out=lambda_N)) 
  lambdas.int<-lambdas*3

  # Fit GRADE
  fits<-GRADE(observations=observations, times_e=times_e,type_reg=type_reg, type_data=type_data,type_basis="monomial", type_smooth=type_smooth, xhat=xhat, xprime=NULL,lambdas=lambdas,nbasis=2)

  # Evaluate the performance 
  recint<-countGraph(graph_est=fits$neighbour.path,graph_true=solutions$graph,self=T)  
  