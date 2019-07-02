##---------------##
##   Estimator   ##
##---------------##

# makes y and x for glmnet
xy <- function(x_smo, inh, interact = FALSE) {
  delta_y <- x_smo[, -1]
  delta_y[, inh] <- 0
  delta_y <- apply(delta_y, 2, diff)
  
  d <- ncol(x_smo) - 1
  MA <- diag(d)
  ss1 <- numint(time = x_smo[, 1], x = x_smo, type = "power", A = MA, B = matrix(0, nrow = 1, ncol = d))
  ss1[, inh] <- 0
  # REcifprocal
  ss2 <- numint(time = x_smo[, 1], x = x_smo, type = "fracpower", A = matrix(0, ncol = d, nrow = d), B = MA)
  ss2[, inh] <- 0
  ss <- cbind(ss1, ss2)
  A = rbind(diag(d), diag(d) * 0) 
  B = rbind(0 * diag(d), diag(d))
  
  # fraction of first orders
  if (interact) {
    AA <- kronecker(matrix(1, nrow = d, ncol = 1), diag(d))
    BB <- kronecker(diag(d), matrix(1, nrow = d, ncol = 1))
    ss3 <- numint(time = x_smo[, 1], x = x_smo, type = "fracpower", A = AA, B = BB)
    ss3[, AA[, inh] != 0 | BB[, inh] != 0] <- 0
    ss <- cbind(ss, ss3)
    A <- rbind(A, AA)
    B <- rbind(B, BB)    
  }

  # MA <- diag(length(x0))
  # colnames(MA) <- rownames(MA) <- names(x0)
  # ss <- numint(time = x_smo[, 1], x = x_smo, type = "power", A = MA, B = MA) #[1, , drop = FALSE])
  coln <- sqrt(colSums(ss^2))
  coln[coln == 0] <- 1
  ss <- t(t(ss) / coln)
  xx <- ss
  
  return(list(delta_y = delta_y, xx = xx, ratobj = list(A = A, B = B)))
}

# Does the aim fit a list of data sets (environments)
aimer <- function(ys, inh, x_smooth = NULL, monotonise = FALSE) {
  if (is.null(x_smooth)) x_smooth <- ys

  # get delta_y and xx from yx
  yxs <- Map(xy, x_smo = x_smooth, inh = as.list(inh))
  
  # bind them
  delta_y <- do.call(rbind, lapply(yxs, function(x) x$delta_y))
  xx <- do.call(rbind, lapply(yxs, function(x) x$xx))
  ratobj <- yxs[[1]]$ratobj
  
  # glmnet part
  species_wise <- apply(delta_y, 2, function(yy) {
    # yy <- yy - mean(yy)
    foo <- glmnet(x = xx, y = yy, standardize = TRUE, intercept = TRUE, alpha = .5, nlambda = 100, dfmax = 22, 
      lambda.min.ratio = 0.05)
    ((t(ratobj$A) != 0 | t(ratobj$B) != 0) %*% (foo$beta != 0)) != 0
  })
  
  
  # combine species
  netw_gs <- sapply(seq_len(100), function(ilam) {
    netw_g <- Matrix(do.call(cbind, lapply(species_wise, function(x) x[, min(ilam, ncol(x))]))) != 0
    netw_g
  })
  
  # monotonise (if one edge proposed, then in it for all following lambdas)
  if (monotonise) {
    netw_gs <- lapply(Reduce('+', lapply(netw_gs, as.vector), accumulate = TRUE), function(x) {
      Matrix(x != 0, nrow = sqrt(length(x)))
    })
  }


  netw_gs
}

estimator <- function(y, inh, stab_remove) {
  # Split y by context
  tout <- y[, 1]
  envr <- cumsum(c(1, diff(tout) < 0))
  ys <- split(y, factor(envr))
  ys <- lapply(ys, matrix, ncol = ncol(y))
  if (length(ys) != length(inh)) stop("ys and inh does not match in length")
  
  # Classic, no stability removes
  aimer_classic <- aimer(ys, inh)
  
  # Stability selection, removes using stab_remove
  stab_rm <- stab_remove[,seq_along(ys),, drop = FALSE]
  yss <- apply(stab_rm, 3, function(i) {
    Map(function(y, remi) {
      y[-remi, , drop = FALSE]
    }, y = ys, remi = split(t(i), seq_len(ncol(i))))
  })
  s <- 0
  aimer_stab <- lapply(yss, function(y) {print(s <<- s + 1); aimer(y, inh)})
  
  # Narrow stability, 
  aimer_stab <- lapply(aimer_stab, function(stab) {
    # first sort out too large models (at most 75% filled, at least 5%)
    stab <- Filter(function(x) {sum(x) <= 0.75 * length(x) & sum(x) >= 0.05 * length(x)}, stab)
    
    # those who appear at some point
    stab <- Reduce('|', stab)
    
    stab
  })
  aimer_stab <- Reduce("+", aimer_stab) / length(aimer_stab)

  
  return(list(aimer_classic = aimer_classic, aimer_stab = aimer_stab))
}



no_neighbours <- function(n, l) {
  # n is number of indeces to choose from, l is number of chosen, no neighbours allowed
  while (TRUE) {
    proposal <- sort(sample(n, size = l))
    if (all(diff(proposal) > 1)) break
  }
  proposal
}


ratmak_generator <- function(prior) {
  # Prior is dxd logical matrix marking allowed parents (row for row)
  
}


