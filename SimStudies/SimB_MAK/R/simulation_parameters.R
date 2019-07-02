##-------------------------##
##  Simulation Parameters  ##
##-------------------------##


if (simname == "sim_01") {
  sim_parameters <- expand.grid(
    alpha = c(1, 2),
    d = c(7, 9, 11),
    k_type = c("a"),
    n = c(10),  
    n_series = c(2, 4, 8),
    sigma = c(.1, .5, 1),
    w_type = c("a"),
    stringsAsFactors = FALSE)
} else {
  stop("simname not recognised")
}

## _factory are used to pass a character value to a function that produces the intended parameter
k_factory <- function(value) {
  switch(value,
    a = function(nz) {return(rep(1, nz))},
    b = function(nz) {
      nz2   <- floor(nz / 2)
      even  <- nz2 == nz / 2 
      if (even) {
        return(c(rep(.5, nz2), rep(5, nz2)))
      } else {
        return(c(rep(.5, nz2), 1, rep(5, nz2)))
      }
    })
}
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
    fine_time <- seq(0, 3, by = 0.01)
    coarse_time <- c(0, sqrt(2)^((5 - n):3))
    chosen <- sapply(coarse_time, function(ct) which.max(fine_time >= ct))
    fine_traj <- numsolve(makobj, time = rep(fine_time, n_series), param = list(k = rate), x0 = x0s)
    coarse_traj <- fine_traj[chosen + 301 * rep(0:(n_series - 1), each = length(chosen)), ]
    # coarse_traj <- numsolve(makobj, time = rep(coarse_time, n_series), param = list(k = rate), x0 = x0s)
    
    return(c(list(makobj = makobj, trajectory = coarse_traj, 
      k = rate, x0 = x0s, w = w_factory(w_type)(d)), sim_param))
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
      mnamp <- sum(equilibrium[namp]) * .5 
      equilibrium[namp] <- equilibrium[namp] * .5
      equilibrium[amp]  <- equilibrium[amp] * (mnamp / sum(equilibrium[amp]) + 1)
      equilibrium
    })
  })
  x0s <- do.call(cbind, x0s)
  mdist <- rowMeans(as.matrix(dist(t(x0s), upper = TRUE, diag = TRUE)))
  mdist <- mdist + sqrt(colSums((x0s - equilibrium)^2)) / 4 # weighting between being far apart from other x0s and equilibrium
  x0s[, sample(rev(order(mdist))[seq_len(n_max)])]
}

