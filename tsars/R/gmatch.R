##---------------------------##
##  Gradient matching tools  ##
##---------------------------##


## Generates the full design matrix for gradient matching (ms is mean process)
X_full <- function(makobj, ms) {
  ff <- apply(ms[, -1, drop = FALSE], 1, field, o = makobj, param = rep(0, nrow(makobj$A)), differentials = TRUE)
  do.call(rbind, lapply(ff, function(f) f$f_dparam[[1]]))
}


## Fit a collection of linear submodels
# x is full design, y is obs, models logical matrix with p columns marking models
fit_models <- function(x, y, models, pos_coef = TRUE) {
  n = length(y); p = ncol(models); m = nrow(models);
  if (!(is.logical(models) | is(models, "lgCMatrix"))) stop("'models' is not logical.")
  if (n != nrow(x)) stop("Number of elements in y does not match number of rows in x.");
  if (p != ncol(x)) stop("Number of columns in models does not match number of columns in x.");

  if (pos_coef) {
    x <- apply(models, 1, function(model) {
      fit <- nnls::nnls(x[, model, drop = FALSE], y)
      return(c(ifelse(any(model), fit$deviance, sum(y^2)), fit$x))
    })
  } else {
    x <- apply(models, 1, function(model) {
      fit <- lm.fit(x[, model, drop = FALSE], y)
      return(c(sum(fit$residuals^2), fit$coefficients))
    })
  }
  x <- unlist(x)
  j <- apply(models, 1, function(model) {
    c(1, which(model) + 1)
    })
  j <- unlist(j)
  i <- apply(models, 1, sum) + 1
  i <- unlist(Map(rep.int, x = seq_len(m), times = i))
  x[!is.finite(x)] <- 0 # Replace NA's and such with 0
  nonzero <- x != 0;  i <- i[nonzero];  j <- j[nonzero];  x <- x[nonzero]

  return(Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(m, p + 1)))
}

## does best subset selection for models of size k=1,...,kmax, x design and y obs
bss <- function(x, y, kmax = 3, intercept = TRUE) {
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(n == length(y))

  ## Standardise (tau = 1 - |y - x beta|^2 / |y|^2 = 1 - |y - x beta|^2)
  if (intercept) y <- y - mean(y)
  ynorm <- sqrt(sum(y^2))
  y <- y / ifelse(ynorm == 0, 1, ynorm)
  xnorm <- sqrt(colSums(x^2))
  xnorm[xnorm == 0] <- 1
  x <- t(t(x) / xnorm)

  ## return object: kmax-x-p matrix with betahat, kmax vector of tau^2
  betahat  <- matrix(0, nrow = kmax, ncol = p)
  tau2  <- numeric(kmax)

  ## k = 1
  yx <- t(y) %*% x
  selected <- which.max(yx)
  tau2[1] <- max(0, yx[selected])^2
  betahat[1, selected] <- max(0, yx[selected])

  ## other k
  if (kmax > 1) {
    for (k in 2:kmax) {
      models    <- combn(x = 1:p, m = k, function(x) {1:p %in% x})
      fits      <- fit_models(x = x, y = y, models = t(models), pos_coef = TRUE)
      selected  <- which.min(fits[, 1])
      tau2[k]   <- 1 - fits[selected, 1]
      betahat[k, ] <- fits[selected, -1]
    }
  }

  ## original scale
  betahat <- t(t(betahat) * ynorm / xnorm)

  return(list(tau2 = tau2, betahat = betahat))
}



## does simulated annealing search for best models of size k=1,...,kmax, with x design and y obs
simanneal <- function(x, y, kmax = 5, niter = 1000, trace = FALSE, cooling = .015, temp = 1, coolfreq = 1, intercept = TRUE) {
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(n == length(y))

  ## Standardise
  if (intercept) y <- y - mean(y)
  ynorm <- sqrt(sum(y^2))
  y <- y / ifelse(ynorm == 0, 1, ynorm)
  xnorm <- sqrt(colSums(x^2))
  xnorm[xnorm == 0] <- 1
  x <- t(t(x) / xnorm)

  ## return object: kmax-x-p matrix with betahat, kmax vector of tau^2
  betahat  <- matrix(0, nrow = kmax, ncol = p)
  tau2  <- numeric(kmax)

  ## k = 1
  yx <- t(y) %*% x
  selected <- which.max(yx)
  tau2[1] <- max(0, yx[selected])^2
  betahat[1, selected] <- max(0, yx[selected]) * ynorm / xnorm[selected]

  ## other ks
  if (kmax > 1) {
    ## best result so far, used in sim anneal
    taumax <- tau2[1]
    selmax <- selected
    betmax <- max(0, yx[selected])

    for (k in 2:kmax) {
      # start from previous solution
      ctau          <- taumax
      selected      <- c(selmax, sample(setdiff(seq_len(p), selmax), 1))
      alternatives  <- setdiff(seq_len(p), selected)

      # draw randoms
      U <- runif(niter)
      cind <- rep_len(seq_along(selected), niter) # sample(seq_along(selected), niter, replace = TRUE)
      aind <- sample(seq_along(alternatives), niter, replace = TRUE)

      # sim anneal
      tm <- temp
      for (iter in seq_len(niter)) {
        # get random neighbour (alternative) fit
        afit <- nnls::nnls(x[, c(selected[-cind[iter]], alternatives[aind[iter]]), drop = FALSE], y)
        atau <- sum(afit$fitted^2)

        if (trace) cat(paste0("iter: ", sprintf("%5d", iter), "   atau: ", sprintf("%.6f",atau), "   ", "ctau: ", sprintf("%.6f",ctau), "   accept prob: ", sprintf("%.10f", min(1, exp((atau - ctau) / tm))), "\n"))

        if (atau > ctau | ifelse(tm == 0, FALSE, U[iter] <= exp((atau - ctau) / tm))) {
          # swap selected and alternative
          tmp <- alternatives[aind[iter]]
          alternatives[aind[iter]] <- selected[cind[iter]]
          selected[cind[iter]] <- tmp

          # change ctau
          ctau <- atau

          # update max stuff
          if (ctau > taumax) {
            taumax <- ctau
            selmax <- selected
            betmax <- afit$x
          }
        }

        # update temperature
        if (iter %% coolfreq == 0) tm <- tm * (1 - cooling)
      }

      # store it
      tau2[k] <- taumax
      betahat[k, selmax] <- betmax * ynorm / xnorm[selmax]
    }
  }

  return(list(tau2 = tau2, betahat = betahat))
}





#' Topological Sensitivity Analysis of MAK-Reaction Systems
#'
#' @description Evaluates algorithm of Babtie et. al (2014, PNAS).
#' @param makobj Object of class \code{\link{mak}}.
#' @param y Matrix holding data. First column is time (must be increasing), remainder are species.
#' @param max_reactions Integer, maximal number of reactions, default is 5. Be very cautious, the number of models searched easily explodes if this is set too high.
#' @param method Character giving method for subset selection. Options: "exact" and "simanneal".
#' @param tout The (augmented) time-points used for gradient matching. If \code{NULL}, the observed time points are used.
#' @details For each species all models based of the possible combinations of parents (at most \code{max_parents} per model) are found and fitted using gradient matching on a gaussian process fitted to data. Note the zero-model with no parents is trivial and thus not fitted.
#' @return A list with an element for each species. Each of these elements is a list with the following elements
#' \itemize{
#'  \item{betas}{Parameter estimates for each possible model for said species in sparse matrix (row = model, col = rate parameter).}
#'  \item{rss}{Vector of residual sum of squares for each model.}
#'  \item{models}{Logical sparse matrix marking what rate parameters where included in each model (row = model, col = rate parameter). }
#'  \item{parents}{Matrix holding indeces of the parents considered in each model (row = model, col = parents). Signs mark activating/inhibiting role of that parent. 0-index marks lag of parent.}
#' }
#' @export
tsars <- function(makobj, y, max_reactions = 3, method = "exact", tout = NULL, ...) {
  ## Checks
  d <- ncol(makobj$A)
  p <- nrow(makobj$A)
  species <- setNames(seq_len(d), colnames(y)[-1])
  if (ncol(y) != d + 1) stop("Number of columns of y must be one larger that number of species in mak.")

  ## Fit GP and get full design
  if (is.null(tout)) tout <- y[, 1]
  n <- length(tout)
  fitted_gp <- gpfit(t(t(y) - c(0, colMeans(y[, -1]))), tout)   # Subtract mean before gpfit, then add again
  fitted_gp$ms <- t(t(fitted_gp$ms) + c(0, colMeans(y[, -1])))
  fitted_gp$ms[, -1] <- pmax(fitted_gp$ms[,-1, drop = FALSE], 1e-8)
  X <- X_full(makobj, fitted_gp$ms)

  ## Sub-design on species-level
  mods <- (makobj$B != makobj$A) # & (makobj$A == 0)

  ## selection:
  if (method == "exact") {
    subsets <- Map(function(x, y, m) {
      ret <- bss(x = x[, m, drop = FALSE], y = y, kmax = max_reactions, ...)
      betahat <- matrix(0, nrow = nrow(ret$betahat), ncol = p)
      betahat[, m] <- ret$betahat
      ret$betahat <- Matrix::Matrix(betahat)
      ret
    }, x = lapply(split(X, rep(species, n)), matrix, ncol = ncol(X)),
      y = split(fitted_gp$dms[, -1, drop = FALSE], rep(species, each = n)),
      m = split(t(mods), species))
  } else if (method == "simanneal") {
    ## Run simulated annealing
    subsets <- Map(function(x, y, m) {
      ret <- simanneal(x = x[, m, drop = FALSE], y = y, kmax = max_reactions, temp = 1,
        cooling = 1-exp(log(0.05) / (50 * d^2)), niter = 100 * d^2, ...)
      betahat <- matrix(0, nrow = nrow(ret$betahat), ncol = p)
      betahat[, m] <- ret$betahat
      ret$betahat <- Matrix::Matrix(betahat)
      ret
    }, x = lapply(split(X, rep(species, n)), matrix, ncol = ncol(X)),
      y = split(fitted_gp$dms[, -1, drop = FALSE], rep(species, each = n)),
      m = split(t(mods), species))
  } else {
    stop("method not among options.")
  }

  ## rank them
  actives <- lapply(subsets, function(s) apply(s$betahat != 0, 1, which))
  taus <- lapply(subsets, function(s) s$tau2)

  ## order them, so the species which introduces largest increase in tau2 is next
  taus <- lapply(taus, function(tau) {
    tau <- diff(c(0, tau))
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
  ord <- sapply(species, function(s) cumsum(ord == s))
  ord <- apply(ord, 1, function(o) {
    sapply(species, function(s) {
      if (o[s] > 0) return(actives[[s]][[o[s]]])
    }, simplify = FALSE)
  })
  ord <- lapply(ord, do.call, what = c)
  ord <- lapply(ord, function(o) sort(unique(o)))

  ## refit estimates
  betahat <- lapply(ord, function(o) {
    nnls::nnls(X[, o, drop = FALSE], as.vector(t(fitted_gp$dms[, -1])))$x
  })
  betahat <- do.call(c, betahat)
  i <- as.integer(do.call(c, ord))
  j <- as.integer(do.call(c, Map(rep, seq_along(ord), lapply(ord, length))))

  j <- j[betahat > 0]
  i <- i[betahat > 0]
  betahat <- betahat[betahat > 0]

  subsets <- Matrix::sparseMatrix(i = i, j = j, x = betahat, dims = c(p, length(ord)))

  return(subsets)
}

