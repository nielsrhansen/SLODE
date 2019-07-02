##----------------##
##   Estimators   ##
##----------------##

## aicc for selecting the right model from maker or aim
aicc <- function(dfs_, rss, n) {
  gcv <- rss / (1 - dfs_ / n)^2
  ind <- which.min(gcv)
  sig <- pmax(rss[ind] / (n - dfs_[ind]), 1e-16)
  aic <- rss / sig + 2 * dfs_ + 2 * dfs_ * (dfs_ + 1) / (n - dfs_ - 1)
  return(aic)
}


## ranking models selected from different species
rank_models <- function(khats) {
  # khats is list (one element per species) each containing a (sparse) matrix of dim (p+1)xm
  # first row is rss (rss under normalisation), rest are parameter estimates
  # columns represents models
  # return is pxm' matrix
  
  ## extract taus
  khats <- lapply(khats, function(k) k[, !duplicated(k[1, ]), drop = FALSE])
  taus <- lapply(khats, function(k)  - k[1, ])
  betas <- lapply(khats, function(k) k[-1,, drop = FALSE] != 0)
  species <- seq_along(khats)
  
  ## order them, so the species which introduces largest increase in tau2 is next
  taus <- lapply(taus, function(tau) {
    tau <- diff(tau)
    tau[tau <= 0] <- 0
    tau
  })
  ord <- numeric()
  while (TRUE) {
    o <- which.max(unlist(lapply(taus, function(tau) tau[1])))
    if (taus[[o]][1] == 0) {
      break
    } else {
      ord <- c(ord, o)
      taus[[o]] <- c(taus[[o]][-1], 0)
    }
  }
  
  # for each species, find out what model to pick
  ord <- sapply(species, function(s) cumsum(ord == s), simplify = FALSE)
  ord <- do.call(cbind, ord)
  mod <- apply(ord, 1, function(o) {
    b <- sapply(species, function(s) {
      if (o[s] > 0) return(betas[[s]][, o[s]])
    }, simplify = FALSE)
    b <- do.call(cbind, b)
    apply(b, 1, any)
  })
  t(unique(t(mod)))
}


preaims <- function(y, x_smooth, makobj) {
  ## im
  ytout <- y[, 1]
  d <- ncol(makobj$A)
  t1 <- Sys.time()
  {
    # get design n' stuff
    # des <- latent::aim_design(makobj$A, makobj$B, x_smooth = x_smooth, tout = ytout)
    des <- episode::imd(makobj, opt(y), x = x_smooth)
    Y   <- des$Y
    X   <- des$X[[1]]
    
    # rescale observation at each time point 
    Ynrm  <- sqrt(colSums(matrix(des$Y, nrow = d)^2))
    if (!all(Ynrm == 0)) {
      Ynrm[Ynrm == 0] <- mean(Ynrm[Ynrm != 0]) 
    } else {
      Ynrm <- rep(1, length(Ynrm))
    }
    Y <- Y / rep(Ynrm, each = d)
    X <- X / rep(Ynrm, each = d)
    
    # find lasso fits
    lassos <- Map(function(x, y) {
      ynorm <- sqrt(sum(y^2))
      if (ynorm == 0) {
        ret <- rbind(tau = 1, Matrix::Matrix(0, nrow = ncol(x), ncol = 1))
      } else {
        y <- y / ifelse(ynorm == 0, 1, ynorm)
        xnorm <- sqrt(colSums(x^2))
        xnorm[xnorm == 0] <- 1
        x <- t(t(x) / xnorm)
        lasso <- glmnet::glmnet(x = x, y = y, intercept = TRUE, lower.limits = 0, pmax = 5, alpha = .75, nlambda = 30)
        ret <- rbind(tau = colSums((y - x %*% as.matrix(lasso$beta))^2), lasso$beta)
      }
      return(ret)
    }, x = lapply(split(X, rep_len(seq_len(d), length.out = nrow(X))), matrix, ncol = ncol(X)),
      y = split(Y, rep_len(seq_len(d), length.out = length(Y))))
    
    mods <- Matrix::Matrix(rank_models(lassos))
    
    # refit
    im1 <- Matrix::Matrix(apply(mods, 2, function(mod) {
      beta <- nnls::nnls(X[, mod, drop = FALSE], Y)$x
      ret <- numeric(ncol(X))
      ret[mod] <- beta
      ret
    }))
  }
  im1_time <- as.numeric(difftime(Sys.time(), t1, units = "sec"))
  
  
  # no rescaling  
  t1 <- Sys.time()
  {
    # get design n' stuff
    des <- episode::imd(makobj, opt(y), x = x_smooth)
    Y   <- des$Y
    X   <- des$X[[1]]
    
    # find lasso fits
    lassos <- Map(function(x, y) {
      ynorm <- sqrt(sum(y^2))
      if (ynorm == 0) {
        ret <- rbind(tau = 1, Matrix::Matrix(0, nrow = ncol(x), ncol = 1))
      } else {
        y <- y / ifelse(ynorm == 0, 1, ynorm)
        xnorm <- sqrt(colSums(x^2))
        xnorm[xnorm == 0] <- 1
        x <- t(t(x) / xnorm)
        lasso <- glmnet::glmnet(x = x, y = y, intercept = TRUE, lower.limits = 0, pmax = 5, alpha = .75, nlambda = 30)
        ret <- rbind(tau = colSums((y - x %*% as.matrix(lasso$beta))^2), lasso$beta)
      }
      return(ret)
    }, x = lapply(split(X, rep_len(seq_len(d), length.out = nrow(X))), matrix, ncol = ncol(X)),
      y = split(Y, rep_len(seq_len(d), length.out = length(Y))))
    
    mods <- Matrix::Matrix(rank_models(lassos))
    
    # refit
    im2 <- Matrix::Matrix(apply(mods, 2, function(mod) {
      beta <- nnls::nnls(X[, mod, drop = FALSE], Y)$x
      ret <- numeric(ncol(X))
      ret[mod] <- beta
      ret
    }))
  }
  im2_time <- as.numeric(difftime(Sys.time(), t1, units = "sec"))
  
  
  ## im with refitting first
  t1 <- Sys.time()
  {
    # get design n' stuff
    des <- episode::imd(makobj, opt(y), x = x_smooth)
    Y   <- des$Y
    X   <- des$X[[1]]
    
    # rescale observation at each time point
    Ynrm  <- sqrt(colSums(matrix(des$Y, nrow = d)^2))
    if (!all(Ynrm == 0)) {
      Ynrm[Ynrm == 0] <- mean(Ynrm[Ynrm != 0])
    } else {
      Ynrm <- rep(1, length(Ynrm))
    }
    Y <- Y / rep(Ynrm, each = d)
    X <- X / rep(Ynrm, each = d)
    
    # find lasso fits
    lassos <- glmnet::glmnet(x = X, y = Y, intercept = FALSE, 
      lower.limits = 0, pmax = 5 * d, alpha = .75, nlambda = 30)
    mods <- lassos$beta != 0
    
    # refit
    im3 <- Matrix::Matrix(apply(mods, 2, function(mod) {
      beta <- nnls::nnls(X[, mod, drop = FALSE], Y)$x
      ret <- numeric(ncol(X))
      ret[mod] <- beta
      ret
    }))
  }
  im3_time <- as.numeric(difftime(Sys.time(), t1, units = "sec"))
  
  
  # no scaling
  t1 <- Sys.time()
  {
    # get design n' stuff
    des <- episode::imd(makobj, opt(y), x = x_smooth)
    Y   <- des$Y
    X   <- des$X[[1]]
    
    # find lasso fits
    lassos <- glmnet::glmnet(x = X, y = Y, intercept = FALSE, 
      lower.limits = 0, pmax = 5 * d, alpha = .75, nlambda = 30)
    mods <- lassos$beta != 0
    
    # refit
    im4 <- Matrix::Matrix(apply(mods, 2, function(mod) {
      beta <- nnls::nnls(X[, mod, drop = FALSE], Y)$x
      ret <- numeric(ncol(X))
      ret[mod] <- beta
      ret
    }))
  }
  im4_time <- as.numeric(difftime(Sys.time(), t1, units = "sec"))
  
  
  return(list(
    k_hat = list(im1 = im1, im2 = im2, im3 = im3, im4 = im4),
    cpu_time = list(im1 = im1_time, im2 = im2_time, im3 = im3_time, im4 = im4_time)
  ))
}

# does some aims and makers
aimmaker <- function(y, x_smooth, makobj) {
  prea <- preaims(y, x_smooth, makobj)
  preb <- preaims(y, NULL, makobj)
  
  ## maker
  fit_maker <- function(k, y, x0 = NULL, ...) {
    model <- k != 0
    if (any(model)) {
      m <- mak(makobj$A[model, , drop = FALSE], makobj$B[model, , drop = FALSE])
      m$rs[[2]]$ctrl$step_max <- 5
      m$rs[[2]]$ctrl$step_screen <- 3
      m$rs[[2]]$ctrl$step_cycle <- 3
      m$rs[[2]]$screen <- TRUE
      m$rs[[2]]$reg_type <- "scad"
      m$rs[[2]]$a <- 3
      m$rs[[1]]$ctrl$step_max <- 5
      m$rs[[1]]$ctrl$step_screen <- 3
      m$rs[[1]]$ctrl$step_cycle <- 3
      e <- opt(y, lambda = 1e-8)
      if (is.null(x0)) x0 <- as.vector(t(y[c(0, which(diff(y[, 1]) < 0)) + 1, -1]))
      mm <- rodeo(m, e, params = list(rate = k[model]), x0 = x0)
      return(list(rss = mm$losses, ks = mm$params[[1]]))
    } else {
      return(list(rss = sum(y[, -1]^2), ks = numeric(0)))
    }
  }
  t1 <- Sys.time()
  
  prefits <- cbind(do.call(cbind, prea$k_hat), do.call(cbind, preb$k_hat))
  nparam  <- apply(prefits != 0, 2, sum)
  prefits <- prefits[, order(nparam), drop = FALSE]
  {
    prefits <- prefits[, !duplicated(t(as.matrix(prefits)) != 0)]
    postfits <- apply(prefits, 2, fit_maker, y = y)
    rss <- unlist(lapply(postfits, getElement, "rss"))
    ii <- apply(prefits != 0, 2, which)
    jj <- Map(rep.int, seq_along(ii), times = lapply(ii, length))
    maker1 <- Matrix::sparseMatrix(i = do.call(c, ii), j = do.call(c, jj), 
      x = unlist(lapply(postfits, function(mm) as.vector(mm$ks))),
      dims = dim(prefits))
  }
  maker1_time <- as.numeric(difftime(Sys.time(), t1, units = "sec"))
  
  
  ## refining maker
  t1 <- Sys.time()
  {
    dfs <- Matrix::colSums(maker1 != 0)
    sel <- sapply(unique(dfs), function(df) {
      which(dfs == df)[which.min(rss[dfs == df])]
    })
    maker2 <- maker1[, sel, drop = FALSE]
  }
  maker2_time <- as.numeric(difftime(Sys.time(), t1, units = "sec")) + maker1_time
  
  
  return(list(
    k_hat = c(prea$k_hat, preb$k_hat, list(maker1 = maker1, maker2 = maker2)),
    cpu_time = c(prea$cpu_time, preb$cpu_time, list(maker1 = maker1_time, maker2 = maker2_time))
  ))
}


# what time points to remove in stability: row = which to remove, col = series, slice = stab replicate
no_neighbours <- function(n, l) {
  # n is number of indeces to choose from, l is number of chosen, no neighbours allowed
  while (TRUE) {
    proposal <- sort(sample(n, size = l))
    if (all(diff(proposal) > 1)) break
  }
  proposal
}

## main function which for given makobj and y produces estimates of rates and initial state for each of the estimators considered, stab_removes holds what to remove in stability select
estimators <- function(y, makobj, stab_removes, stability = TRUE) {
  ## Expand time
  tout <- y[, 1]
  brks <- c(1, c(which(diff(tout) < 0), length(tout)) + 1)
  tout <- Map(seq, tout[brks[-length(brks)]], tout[(brks[-1] - 1)], length.out = (5 * (diff(brks) - 1) + 1))#by = tout %>% diff %>% extract({. > 0}) %>% min)
  if (is(tout, "list")) tout <- do.call(c, tout)
  
  
  ### Estimation on full data set ###
  
  # tsa
  t1 <- Sys.time()
  tsars_res <- tsars::tsars(makobj = makobj, y = y, max_reactions = 5, method = "exact")
  tsars_time <- as.numeric(difftime(Sys.time(), t1, units = "sec"))
  
  # GP
  t1 <- Sys.time()
  gp_res <- gpfit(y, tout = tout)
  gp_res <- gp_res$ms
  gp_res[, -1] <- pmax(gp_res[, -1], 1e-3)
  gp_time <- as.numeric(difftime(Sys.time(), t1, units = "sec"))
  
  # MCP
  m <- mak(A = makobj$A, B = makobj$B, r = reg(reg_type = "mcp"))
  t1 <- Sys.time()
  mm <- rodeo(m, opt(y, nlambda = 40, lambda_min_ratio = 1e-5), 
    params = list(rate = rep(0, nrow(m$A))), x0 = NULL)
  mcp_time <- as.numeric(difftime(Sys.time(), t1, units = "sec"))
  mcp_res <- mm$params[[1]]
  
  # SCAD
  m <- mak(A = makobj$A, B = makobj$B, r = reg(reg_type = "scad"))
  t1 <- Sys.time()
  mm <- rodeo(m, opt(y, nlambda = 40, lambda_min_ratio = 1e-5), 
    params = list(rate = rep(0, nrow(m$A))), x0 = NULL)
  scad_time <- as.numeric(difftime(Sys.time(), t1, units = "sec"))
  scad_res <- mm$params[[1]]
  
  # aimmaker
  aimmaker_clas <- aimmaker(y = y, x_smooth = gp_res, makobj = makobj)

  
  ### Stability ###
  if (stability) {
    ## Grouping of series    
    series  <- factor(cumsum(c(0, diff(y[, 1]) < 0)))
    ys  <- lapply(split(y, series), matrix, ncol = ncol(y)) # split by series
    ys  <- apply(stab_removes[,seq_len(length(ys)), , drop = FALSE], 3, function(i) {
      list(
        do.call(rbind, 
          Map(function(y, remi) {y[-remi, , drop = FALSE]}, y = ys, remi = split(t(i), seq_len(ncol(i))))))
    })
    ys <- lapply(ys, function(y) y[[1]]) # had to do/undo list from above, because of apply
    gp_ress <- lapply(ys, function(y) {
      gp_res <- gpfit(y, tout = tout)
      gp_res <- gp_res$ms
      gp_res[, -1] <- pmax(gp_res[, -1], 1e-3)
      
      if (any(!is.finite(gp_res))) {
        spl <- factor(cumsum(c(0, diff(gp_res[, 1]) < 0)))
        gps <- lapply(split(gp_res, spl), matrix, ncol = ncol(gp_res))
        gps <- lapply(gps, function(gp_) {
          cbind(time = gp_[, 1], apply(gp_[, -1, drop = FALSE], 2, function(z) {
            approx(x = gp_[, 1], y = z[], xout = gp_[, 1], rule = 2)$y
          }))
        })
        gp_res <- do.call(rbind, gps)
      }
      gp_res
    })
    
    
    # aimmaker
    aimmaker_makfixed <- function(y, x_smooth) {
      aimmaker(y = y, x_smooth = x_smooth, makobj = makobj)
    }
    aimmaker_stab <- Map(aimmaker_makfixed, y = ys, x_smooth = gp_ress)
  }
  
  
  clas <- list(k_hats = c(tsars = list(tsars_res), mcp = mcp_res, scad = scad_res, aimmaker_clas$k_hat),
    cpu_time = c(tsars = tsars_time, mcp = mcp_time, scad = scad_time, unlist(aimmaker_clas$cpu_time)))
  if (stability) {
    ret <- list(
      clas = clas,
      stab = list(k_hats = lapply(aimmaker_stab, getElement, "k_hat"),
        cpu_time = apply(do.call(rbind, lapply(aimmaker_stab, function(l) unlist(l$cpu))), 2, max)
      ))
  } else {
    ret <- list(clas = clas)
  }

  return(ret) 
}


