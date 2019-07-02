##-------------------------##
##  Simulation Parameters  ##
##-------------------------##


if (simname == "sim_0.1.001.001") {
  sim_parameters <- expand.grid(
    n = c(30),  
    n_series = c(2, 4, 8, 16),
    sigma = c(.1, .5, 1),
    w_type = c("a"),
    stringsAsFactors = FALSE)
} else if (simname == "sim_0.1.001.002") {
  sim_parameters <- expand.grid(
    n = c(30),  
    n_series = c(2, 4, 8, 16),
    sigma = c(.1, .5, 1),
    w_type = c("a"),
    stringsAsFactors = FALSE) 
} else if (simname == "sim_0.1.001.003") {
  sim_parameters <- expand.grid(
    n = c(30),  
    n_series = c(2, 4, 8, 16),
    sigma = c(.1, .5, 1),
    w_type = c("a"),
    stringsAsFactors = FALSE) 
} else if (simname == "sim_0.1.001.004") {
  sim_parameters <- expand.grid(
    n = c(30),  
    n_series = c(5, 10, 15, 20),
    sigma = c(.1, .25, .5),
    w_type = c("a"),
    stringsAsFactors = FALSE) 
} else {
  stop("simname not recognised")
}

## _factory are used to pass a character value to a function that produces the intended parameter
w_factory <- function(value) {
  switch(value,
    a = function(d) {return(rep(1, d))},
    b = function(d) {
      d2    <- floor(d / 2)
      even  <- d2 == d / 2
      if (even) {
        return(c(rep(.5, d2), rep(2, d2)))
      } else {
        return(c(rep(.5, d2), 1, rep(2, d2)))
      }
    })
}




## get the true rate parameters and trajectory of system 
get_truth <- function(sim_param, makobj, x0s) {
  ## sim_param is a row in sim_parameters
  with(sim_param, {
    p <- nrow(makobj$A)
    
    ## rate parameter
    rate  <- numeric(p)
    nz    <- ceiling(alpha * d)
    rate[seq_len(nz)] <- k_factory(k_type)(nz) 
    
    ## Adapt time scale
    coarse_time <- seq(0, 10, by = .1)
    coarse_traj <- numsolve(makobj, time = rep(coarse_time, n_series), param = rate, x0 = x0s)
    drift <- coarse_traj %>% pmax(0) %>%
      extract(, -1, drop = FALSE) %>%
      apply(1, field, o = makobj, param = rate) %>% '^'(2) %>%
      apply(2, sum)
    drift <- split(drift, rep(seq_len(n_series), each = length(coarse_time)))
    stop_ind <- lapply(drift, {
      . %>% '/'(extract(., 1)) %>% '>'(1e-3) %>% {ifelse(all(.), length(.), which.min(.))}
    }) %>% unlist %>% max
    fine_time <- seq(0, coarse_time[stop_ind], length.out = n)
    if (stop_ind <= 1) { stop_ind <- 2 }
    fine_time <- c(0, sqrt(2)^((5 - n):3)) * coarse_time[stop_ind]
    fine_traj <- numsolve(makobj, time = rep(fine_time, n_series), param = rate, x0 = x0s)
    
    return(c(list(makobj = makobj, trajectory = fine_traj, k = rate, x0 = x0s, w = w_factory(w_type)(d)), sim_param))
  })
}



## returns y (trajectory with noise)
get_data <- function(truth, epsilon) {
  with(truth, {
    y <- trajectory
    y[, -1]  <- y[, -1] + epsilon
    return(y)
  })
}


## generate x0s from equilibrium
get_x0s <- function(equilibrium, d, knockups = c(2, 3), n_max = 16) {
  stopifnot(all(knockups < d) & all(knockups > 0))
  which_amplify <- sapply(knockups, combn, x = d, simplify = FALSE)
  if (Reduce("+", lapply(which_amplify, ncol)) < n_max) stop("Number of knockups not enough to generate n_max x0s.")
  x0s <- lapply(which_amplify, function(w) {
    apply(w, 2, function(amp) {
      namp  <- setdiff(seq_len(d), amp)
      mnamp <- sum(equilibrium[namp]) * .9 
      equilibrium[namp] <- equilibrium[namp] * .1
      equilibrium[amp]  <- equilibrium[amp] * (mnamp / sum(equilibrium[amp]) + 1)
      equilibrium
    })
  })
  x0s <- do.call(cbind, x0s)
  mdist <- rowMeans(as.matrix(dist(t(x0s), upper = TRUE, diag = TRUE)))
  mdist <- mdist + sqrt(colSums((x0s - equilibrium)^2)) / 4 # weighting between being far apart from other x0s and equilibrium
  x0s[, sample(rev(order(mdist))[seq_len(n_max)])]
}

# #  Truth parameters  
# #--------------------
# 
# ## Truth parameters are used to create the true parameters
# truth_parameters <- expand.grid(
#   alpha = c(.05, .1, .2),
#   d     = c(3, 5, 7),
#   k     = c("a", "b"),
#   stringsAsFactors = FALSE)
# truth_parameters$seed <- seq_along(truth_parameters$alpha)
# 
# ## _factory are used to pass a character value to a function that produces the intended parameter
# k_factory <- function(value) {
#   switch(value,
#     a = function(nz) {return(rep(1, nz))},
#     b = function(nz) {
#       nz2   <- floor(nz / 2)
#       even  <- nz2 == nz / 2 
#       if (even) {
#         return(c(rep(.5, nz2), rep(5, nz2)))
#       } else {
#         return(c(rep(.5, nz2), 1, rep(5, nz2)))
#       }
#     })
# }
# 
# ## get the true rate parameters, x0, epsilon-matrix, mak-obj and trajectory of system 
# get_truth <- function(truth) {
#   ## truth is a row in truth_parameters
#   with(truth, {
#     m     <- mak_enzyme(d, empty = FALSE)
#     p     <- nrow(m$A)
#     shuff <- sample(seq_len(p))
#     m$A   <- m$A[shuff, ]; m$B <- m$B[shuff, ];
#     
#     rate  <- numeric(p)
#     nz    <- ceiling(alpha * p)
#     rate[seq_len(nz)] <- k_factory(k)(nz) 
#     x0  <- runif(d, min = 1, max = 10)
#     
#     n_max   <- max(data_parameters$n)
#     epsilon <- matrix(rnorm(d * n_max), nrow = n_max, ncol = d)
#     
#     coarse_time <- seq(0, 10, by = .1)
#     coarse_traj <- solve(m, time = coarse_time, k = rate, x0 = x0)
#     relative_drift <- coarse_traj %>% pmax(0) %>%
#       extract(, -1, drop = FALSE) %>%
#       apply(1, f_mak, makobj = m, k = rate) %>% '^'(2) %>%
#       apply(2, sum) %>% '/'(extract(., 1)) 
#     stop_ind  <- relative_drift %>% '>'(1e-3) %>% {ifelse(all(.), length(.), which.min(.))}
#     fine_time <- seq(0, coarse_time[stop_ind], length.out = n_max)
#     fine_traj <- solve(m, time = fine_time, k = rate, x0 = x0)
#     
#     return(list(makobj = m, epsilon = epsilon, trajectory = fine_traj, k = rate, x0 = x0))
#   })
# }
# 
# 
# 
# 
# 
# #  Data parameters  
# #-------------------
# 
# ## Data parameters are used to create data from the true parameters
# data_parameters <- expand.grid(
#   n     = c(10, 20, 50, 100, 200),
#   sigma = c(.1, .5, 1),
#   w     = c("a", "b"),
#   stringsAsFactors = FALSE)
# 
# ## _factory are used to pass a character value to a function that produces the intended parameter
# w_factory <- function(value) {
#   switch(value,
#     a = function(d) {return(rep(1, d))},
#     b = function(d) {
#       d2    <- floor(d / 2)
#       even  <- d2 == d / 2
#       if (even) {
#         return(c(rep(.5, d2), rep(2, d2)))
#       } else {
#         return(c(rep(.5, d2), 1, rep(2, d2)))
#       }
#     })
# }
# 
# ## Get the desired data set
# get_data <- function(data, param) {
#   ## data is row of data_parameters
#   ## param is returned from get_truth
#   with(data, {
#     obs <- param$trajectory
#     tin <- obs[, 1]
#     tou <- seq(min(tin), max(tin), length.out = n)
#     obs <- obs[,-1, drop = FALSE]
#     
#     e <- param$epsilon
#     e <- t(t(e) * sigma * w_factory(w)(ncol(e)))
#     obs %<>% '+'(e)
#     
#     obs %<>% apply(2, approx, x = tin, xout = tou) %>% lapply(. %>% '$'(y)) %>% do.call(cbind, .)
#     
#     return(cbind(Time = tou, obs))
#   })
# }
# 
# 



#  Informals Tests
#-------------------
stopifnot(
  identical(w_factory("a")(5), c(1, 1, 1, 1, 1)),
  identical(w_factory("a")(6), c(1, 1, 1, 1, 1, 1)),
  identical(w_factory("b")(5), c(.5, .5, 1, 2, 2)),
  identical(w_factory("b")(6), c(.5, .5, .5, 2, 2, 2))
)
