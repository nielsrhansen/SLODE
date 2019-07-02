rm(list = ls())


## Libraries ##
library(magrittr)
library(plyr)
library(reshape)
library(ggplot2)
library(glmnet)
library(foreach)
library(doParallel)
library(episode)

## Sources ##
simname <- "simprior_3"
source("R/Hynne_model.R")
source("simulate_data.R")

if (simname == "simprior_1") {
  sim_parameters <- expand.grid(
    sigma = c(0.1, 0.5, 1),
    a = c("correct", "extra"),
    stringsAsFactors = FALSE)
  sim_parameters <- do.call(rbind, Map(cbind, replicate(10, sim_parameters, simplify = FALSE), eps = as.list(1:10)))
  sim_parameters$eps <- as.factor(sim_parameters$eps)
} else if (simname == "simprior_2") {
  sim_parameters <- expand.grid(
    sigma = c(0.1, 0.25, 0.5),
    a = c("correct", "extra", "no_prior"),
    stringsAsFactors = FALSE)
  sim_parameters <- do.call(rbind, Map(cbind, replicate(10, sim_parameters, simplify = FALSE), eps = as.list(1:10)))
  sim_parameters$eps <- as.factor(sim_parameters$eps)
} else if (simname == "simprior_3") {
  sim_parameters <- expand.grid(
    sigma = c(0.1, 0.25, 0.5),
    a = c("correct", "extra", "no_prior"),
    stringsAsFactors = FALSE)
  sim_parameters <- do.call(rbind, Map(cbind, replicate(100, sim_parameters, simplify = FALSE), eps = as.list(1:100)))
  sim_parameters$eps <- as.factor(sim_parameters$eps)
} else {
  stop("simname not recognised")
}


# True network
rat <- ratmak(A = A_correct, C = CC, s = solver(step_max = 100000, h_init = 0.001))
netw <- field(o = rat, x = x0, param = list(theta1 = K1, theta2 = K2), differentials = TRUE)
netw <- netw$f_dx != 0

## From trajectory to discrete time points
trajs0 <- lapply(trajs, function(x) {
  tt <- exp(seq(0, log(nrow(x) - 10), length.out = 30))
  tt <- round(tt)
  tt <- tt + cumsum(c(0, diff(tt) <= 0)) # if two neigbours are identical, add one
  x[tt , ]
})

## Noise
set.seed(27 + 09 + 1990)
Neps <- length(levels(sim_parameters$eps))
eps <- array(rnorm(30 * rat$d * Neps), dim = c(30, rat$d, Neps))

## Comp. setups ##
cl <- makeCluster(min(60, nrow(sim_parameters)), outfile = "")
registerDoParallel(cl)

# i <- 7
## Actual simulation ##
res <- foreach(i = seq_len(nrow(sim_parameters)), .packages = c('episode', 'glmnet', 'magrittr', 'tsars')) %dopar% {
  sim_param <- sim_parameters[i, ]
  sigma <- sim_param$sigma
  if (sim_param$a == "correct") {
    A <- A_correct
  } else if (sim_param$a == "extra") {
    A <- A_extra
  } else if (sim_param$a == "no_prior") {
    dd <- rat$d
    AB <- list(
      A = kronecker(rbind(0, diag(dd)), c(0, rep(1, dd))),
      B = kronecker(c(0, rep(1, dd)), rbind(0, diag(dd)))
    )
  } else {
    stop("a not recognised")
  }
  if (sim_param$a == "no_prior") {
    # Use approximate fit
    odeobj <- mak(A = AB$A, B = AB$B, s = solver(step_max = 100000, h_init = 0.001))
  } else {
    # Use ratmak class for exact model class
    odeobj <- ratmak(A = A, C = CC, s = solver(step_max = 100000, h_init = 0.001))
  }
  
  # Log file
  log_dir <- paste0("Log/", simname, "/log_", i)
  cat(paste0("Started simulation setup ", i, "\n"), file = log_dir)
  
  # Get un-noised data and noise
  # y_raw <- do.call(rbind, ys0[seq_len(E)])
  ss <- 0
  yse <- lapply(trajs0, function(x) {
    ss <<- ss + 1
    x[, -1] <- x[, -1] + sigma * eps[,,as.numeric(sim_param$eps)]
    x
  })
  
  ss <- 0
  # yse <- yse[1:2]
  # y <- yse[[1]]
  # stab_remove <- replicate(50, replicate(50, no_neighbours(n = 30, l = 5)))
  to.res <- lapply(yse, function(y) {
    ss <<- ss + 1
    
    if (sim_param$a == "no_prior") {
      # Use approximate fit
      odeobj_ <- mak(A = AB$A[, -ss], B = AB$B[, -ss], s = solver(step_max = 100000, h_init = 0.001))
    } else {
      # Use ratmak class for exact model class
      odeobj_ <- ratmak(A = A[, -ss], C = CC[, -ss], s = solver(step_max = 100000, h_init = 0.001))
    }
    
    a <- episode::aim(odeobj_, opt(y = y[, -(ss+1)]), adapts = NULL)
     
    if (sim_param$a == "no_prior") {
      ret <- list(params=a$params[[1]], odeobj=odeobj)
      
    } else {
      # Run rodeo and collect loss value and network
      ret <- sapply(seq_len(ncol(a$params$theta1)), function(i) {
        # "exclude" not available for ratmak class, so simple work-around via box constraints and remove unused complexes #
        theta1 <- matrix(a$params$theta1[, i], ncol = nrow(A))
        theta2 <- matrix(a$params$theta2[, i], ncol = nrow(A))
        reduced <- apply(theta1 != 0 | theta2 != 0, 2, any)
        if (any(reduced)) {
          lower1 <- c(0, -Inf)[as.numeric(theta1[, reduced] != 0) + 1]
          upper1 <- c(0, Inf)[as.numeric(theta1[, reduced] != 0) + 1]
          upper2 <- c(0, Inf)[as.numeric(theta2[, reduced] != 0) + 1]
          rattrod <- ratmak(A = A[reduced, -ss, drop = FALSE], C = CC[, -ss], s = solver(step_max = 100000, h_init = 0.001),
            r1 = reg("none", lower = lower1, upper = upper1, step_max = 20),
            r2 = reg("none", upper = upper2, step_max = 20))
          rod <- rodeo(rattrod, opt(y = y[, -(ss+1)]), x0 = pmax(a$x0s[, 1], 1e-5),
            params = list(
              theta1 = as.vector(theta1[, reduced]),
              theta2 = as.vector(theta2[, reduced])))

          # If error code, use integral matching instead
          if (any(rod$codes > 2)) {
            ntw <- field(odeobj_, x = x0[-ss], 
              param = list(theta1 = a$params$theta1[, i],
              theta2 = a$params$theta2[, i]), differentials = TRUE)$f_dx != 0
          } else {
            ntw <- field(rattrod, x = x0[-ss],
              param = list(theta1 = rod$params$theta1[, 1],
                theta2 = rod$params$theta2[, 1]), differentials = TRUE)$f_dx != 0
            attr(ntw, "loss") <- rod$losses
          }
        } else {
          # empty model
          ntw <- Matrix::Matrix(0, 21, 21)
          colnames(ntw) <- rownames(ntw) <- names(x0[-ss])
        } #
        ntw
      })
    }
    return(ret)
  })
  
  
  # Save
  save(to.res, file = paste0("Results/", simname, "/res_", i, ".RData"))
  cat(paste0("Saved result"), "\n",
    file = log_dir, append = TRUE)
  
  to.res
}


stopCluster(cl)




