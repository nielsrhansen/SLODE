rm(list = ls())


## Libraries ##
library(tsars)
library(episode)
library(magrittr)
library(foreach)
library(doParallel)
library(plyr)
# library(ggplot2); library(reshape)

## Sources ##
simname <- "sim_01"
source("R/simulation_parameters.R")
source("R/estimators.R")
source("R/assessment.R")


mak_creator <- function(d) {
  A <- sapply(seq_len(d), function(i) {
    t(sapply(setdiff(seq_len(d), i), function(j) {
      ret <- numeric(d)
      ret[c(i, j)] <- 1
      ret
    }))
  }, simplify = FALSE)
  A <- do.call(rbind, A)
 
  B <- sapply(seq_len(d), function(i) {
    t(sapply(setdiff(seq_len(d), i), function(j) {
      ret <- numeric(d)
      ret[i] <- 2
      ret
    }))
  }, simplify = FALSE)
  B <- do.call(rbind, B)
  
  mak(A = A, B = B, s = solver(step_max = 10000))
}



## Sample truths ##
makobjs_param <- expand.grid(alpha = unique(sim_parameters$alpha),
  d = unique(sim_parameters$d),
  n_series = unique(sim_parameters$n_series),
  k_type = unique(sim_parameters$k_type))
makobjs_param$B <- c(100, 100, 100)[log2(makobjs_param$n_series)]

# sim_param <- makobjs_param[10, ]
makobjs <- apply(makobjs_param, 1, function(sim_param) {
  d <- as.numeric(sim_param["d"])
  alpha <- as.numeric(sim_param["alpha"])
  n_series <- as.numeric(sim_param["n_series"])
  B <- as.numeric(sim_param["B"])
  k_type <- as.character(sim_param["k_type"])
  m_orig <- mak_creator(d)
  
  set.seed(d)
  replicate(B, {
    ## maks
    m     <- m_orig
    p     <- nrow(m$A)
    
    first <- with(m, {
      inds <- apply((B - A) != 0, 2, which)
      
      io <- sample(d)
      first <- integer(0)
      
      while (TRUE) {
        for (i in io) {
          if (length(first) > 0) {
            legal <- apply(A[inds[, i],,drop = FALSE], 1, function(a) {
              !any(apply(t(A[first,, drop = FALSE]) == a, 2, all))
            })
            sel <- inds[legal, i]
          } else {
            sel <- inds[, i]
          }
          
          if (length(sel) >= alpha) {
            proposal <- sample(sel, size = alpha)
            first <- c(first, proposal)
          }
        }
        if (length(first) == d * alpha) {
          break
        }
      }

      first
    })
    
    
    first <- as.vector(first)
    first <- first[sample(seq_along(first))]
    shuff <- c(first, setdiff(seq_len(p), first)) # make sure they are first (makes it easier to read off results)
    m$A   <- m$A[shuff, ]; m$B <- m$B[shuff, ];
    
    ## x0s (first find equilibrium)
    rate  <- numeric(p)
    nz    <- ceiling(alpha * d)
    rate[seq_len(nz)] <- k_factory(k_type)(nz) 
    coarse_time <- seq(0, 10, by = .1)
    coarse_traj <- numsolve(m, time = rep(coarse_time, n_series), param = list(rate = rate), 
      x0 = rep(5, d))
    eq  <- pmax(coarse_traj[nrow(coarse_traj), -1], 1e-1)
    x0s <- get_x0s(equilibrium = 5 * eq / mean(eq), d = d, knockups = c(1, 2, 3, 4), n_max = max(makobjs_param$n_series))
    
    list(makobj = m, x0s = x0s[, seq_len(n_series)])
  }, simplify = FALSE)
})

epsilons_param <- expand.grid(
  d = unique(sim_parameters$d),
  n = unique(sim_parameters$n),
  n_series = unique(sim_parameters$n_series),
  sigma = unique(sim_parameters$sigma),
  w_type = unique(sim_parameters$w_type)
)
epsilons_param$B <- c(100, 100, 100)[log2(epsilons_param$n_series)]
epsilons <- apply(epsilons_param, 1, function(sim_param) {
  set.seed(sim_param["d"])
  d <- as.numeric(sim_param["d"])
  n <- as.numeric(sim_param["n"])
  n_series <- as.numeric(sim_param["n_series"])
  B <- as.numeric(sim_param["B"])
  sigma <- as.numeric(sim_param["sigma"])
  w <- w_factory(sim_param["w_type"])(d)
  epsi <- array(rnorm(n * d * max(epsilons_param$n_series) * B) * 
      sigma * w,
    dim = c(d, n * max(epsilons_param$n_series), B))
  epsi <- epsi[, seq_len(n * n_series) ,]
  alply(epsi, 3, t)
})
# row = what to remove, col = series, slice = stab subset
stab_remove <- replicate(50, replicate(50, no_neighbours(n = max(sim_parameters$n), l = 5)))


## Comp. setups ##
cl <- makeCluster(min(3, nrow(sim_parameters)), outfile = "")
registerDoParallel(cl)

# i <- 1
## Actual simulation ##
res <- foreach(i = seq_len(nrow(sim_parameters)), .packages = c('episode', 'magrittr', 'tsars')) %dopar% {
  sim_param <- sim_parameters[i, ]

  log_dir <- paste0("Log/", simname, "/log_", i)
  cat(paste0("Started simulation setup ", i, "\n"), file = log_dir)
  
  # Extract all the randomness:
  m_i <- which(
    makobjs_param$alpha == sim_param$alpha & 
      makobjs_param$d == sim_param$d &
      makobjs_param$n_series == sim_param$n_series &
      makobjs_param$k_type == sim_param$k_type)
  ms  <- makobjs[[m_i]]
  x0  <- lapply(ms, getElement, "x0s")
  ms  <- lapply(ms, getElement, "makobj")
  e_i <- which(epsilons_param$d == sim_param$d & 
      epsilons_param$n_series == sim_param$n_series &
      epsilons_param$sigma == sim_param$sigma &
      epsilons_param$w_type == sim_param$w_type)
  eps <- epsilons[[e_i]]

  cat(paste0("Drew randomness\n"), file = log_dir, append = TRUE)
  
  
  ## Get truths
  trus  <- Map(function(makobj, x0s) get_truth(sim_param, makobj, x0s), makobj = ms, x0s = x0)
  cat(paste0("Got the truths\n"), file = log_dir, append = TRUE)
  datas <- Map(get_data, truth = trus, epsilon = eps)
  val_datas <- Map(get_data, truth = trus, epsilon = c(eps[-1], eps[1]))
  cat(paste0("Got the data\n"), file = log_dir, append = TRUE)
  
  ## Do the estimation and validate
  s <- 0
  ests    <- Map(function(y, makobj) {
    s <<- s + 1
    cat(paste0("Ran estimator no: ", s), "\n", file = log_dir, append = TRUE)
    try(estimators(y = y, makobj = makobj, stab_removes = stab_remove, stability = FALSE), silent = TRUE)
    }, y = datas, makobj = lapply(trus, getElement, "makobj"))
  tt <- table(unlist(lapply((lapply(ests, is)), paste, collapse = ",")))
  cat(paste0("Ran estimator. Breakdown: \n"), names(tt), "\n", tt, "\n",
    file = log_dir, append = TRUE)
  lapply(ests, function(l) {
    if (is(l, "try-error")) {
      cat(as.character(l), "\n", file = log_dir, append = TRUE)
    }
  })
  
  ## Assess
  # est_return <- ests[[1]]; truth = trus[[1]]; validation_y = val_datas[[1]]
  s <- 0
  to.res  <- Map(function(e, tr, val) {
    s <<- s + 1
    cat(paste0("Ran assessment no: ", s), "\n", file = log_dir, append = TRUE)
    try(assess(est_return = e, truth = tr, validation_y = val, stability = FALSE), silent = TRUE)
  }, e = ests, tr = trus, val = val_datas)
  tt <- table(unlist(lapply((lapply(to.res, is)), paste, collapse = ",")))
  cat(paste0("Ran assessment. Breakdown: \n"), names(tt), "\n", tt, "\n",
    file = log_dir, append = TRUE)
  lapply(to.res, function(l) {
    if (is(l, "try-error")) {
      cat(as.character(l), "\n", file = log_dir, append = TRUE)
    }
  })
  

  ## Save the result just in case
  save(to.res, file = paste0("Results/", simname, "/res_", i, ".RData"))
  
  ## Return
  to.res
}





