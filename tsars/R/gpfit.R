##---------------------------##
##   Gaussian Process Fits   ##
##---------------------------##

## changed to prevent underflow (using hyperbolics)
p2w   <- function(p) {return(log(sinh(p / 2)) + p / 2 + log(2))}
w2p   <- function(w) {return(log(cosh(w / 2)) + w / 2 + log(2))}
# w2pII <- function(w) {return(log(2) + log(cosh(w / 2)) + w / 2)}
# w2pII(-34)
# w2p(-34)


## profile -loglikelihood
pll <- function(w, y, K0, return_Ky = FALSE) {
  p <- w2p(w)
  K <- exp(-K0 / (2 * p[1]))
  diag(K) <- diag(K) + p[2]
  qrK <- try(qr(K), silent = TRUE)
  if (is(qrK, "qr")) {
    Ky  <- try(solve(qrK, y), silent = TRUE)
    if (is(Ky, "numeric")) {
      if (return_Ky) {
        return(Ky)
      } else {
        return(log(abs(sum(y * Ky))) + log(abs(prod(diag(qrK$qr))) + 1e-20) / nrow(K))
      }
    }
  }
  return(1e200)
}


## gaussian process fit
gp <- function(y, xin, xout) {
  K0  <- outer(xin, xin, "-")^2
  O   <- try(nlm(pll, c(0, 0), y = y, K0 = K0), silent = TRUE)
  if (!is(O, "try-error")) {
    p     <- w2p(O$estimate)
    Ky    <- pll(w = O$estimate, y = y, K0 = K0, return_Ky = TRUE)
    if (length(Ky) != 1) {
      Kout  <- outer(xout, xin, "-")
      m     <- as.vector(exp(- Kout^2 / (2 * p[1])) %*% Ky)
      dm    <- as.vector((-Kout * exp(- Kout^2 / (2 * p[1])) / p[1]) %*% Ky)
    } else {
      # return linear interpolation with derivate averaged between neighbours (simpson type of diff estimate)
      m   <- approx(x = xin, y = y, xout = xout)$y
      dm  <- diff(m) / diff(xout)
      dm  <- c(dm[1], (dm[-1] + dm[-length(dm)]) / 2, dm[length(dm)])
    }
  } else {
    # return linear interpolation with derivate averaged between neighbours (simpson type of diff estimate)
    m   <- approx(x = xin, y = y, xout = xout)$y
    dm  <- diff(m) / diff(xout)
    dm  <- c(dm[1], (dm[-1] + dm[-length(dm)]) / 2, dm[length(dm)])
    p   <- NA
  }
  return(list(m = m, dm = dm, p = p))
}


#' Gaussian Process Fits
#' @description Returns mean estimates of gaussian process and its derivative
#' @param ys Data stored column-wise, first column must be time and each consequtive column is abundance of a species.
#' @param tout Time vector of desired timepoints at which to evaluate mean of gaussian process and its derivative.
#' @return A list with two elements:
#' \item{ms}{Matrix of mean values, stored column wise.}
#' \item{dms}{Matrix of derivate mean values, stored columns wise.}
#' @details
#' Whenever time (first column of \code{ys}) decreases, the system restarts. Hence the data is assumed generated from s different contexts, where s - 1 is the number of decreases in the time vector. Each context has its own initial condition and parameter vector. The number of decreases in time in \code{ys} must match the number of decreases in \code{tout}.
#' @export
gpfit <- function(ys, tout) {
  tin <- ys[, 1]
  obs <- ys[, -1, drop = FALSE]

  brksi <- c(1, c(which(diff(tin) < 0), length(tin)) + 1)
  brkso <- c(1, c(which(diff(tout) < 0), length(tout)) + 1)

  if (length(brksi) != length(brkso)) stop("Time_in in y does not have same number of jumps as tout.")

  ms <- list(); dms <- list(); ps <- list()
  for (i in seq_len(length(brksi) - 1)) {
    indi <- seq(brksi[i], brksi[i + 1] - 1)
    indo <- seq(brkso[i], brkso[i + 1] - 1)
    res <- apply(obs[indi, , drop = FALSE], 2, gp, xin = tin[indi], xout = tout[indo])

    ms[[i]]  <- lapply(res, function(r) r$m)
    ms[[i]]  <- cbind(Time = tout[indo], do.call("cbind", ms[[i]]))
    dms[[i]] <- lapply(res, function(r) r$dm)
    dms[[i]] <- cbind(Time = tout[indo], do.call("cbind", dms[[i]]))
    ps[[i]]  <- do.call(cbind, lapply(res, function(r) r$p))
  }

  return(list(ms = do.call("rbind", ms), dms = do.call("rbind", dms), p = do.call("rbind", ps)))
}


#' Gaussian Process Evaluation
#' @description Returns mean of gaussian process and its derivative
#' @param ys Data stored column-wise, first column must be time and each consequtive column is abundance of a species.
#' @param tout Time vector of desired timepoints at which to evaluate mean of gaussian process and its derivative.
#' @param p A vector of parameters to feed to the fit. First coordinate is kernel variances, second observation variances.
#' @return A list with two elements:
#' \item{ms}{Matrix of mean values, stored column wise.}
#' \item{dms}{Matrix of derivate mean values, stored columns wise.}
#' @details
#' Whenever time (first column of \code{ys}) decreases, the system restarts. Hence the data is assumed generated from s different contexts, where s - 1 is the number of decreases in the time vector. Each context has its own initial condition and parameter vector. The number of decreases in time in \code{ys} must match the number of decreases in \code{tout}.
#' @export
gpe <- function(ys, tout, p) {
  tin <- ys[, 1]
  obs <- ys[, -1, drop = FALSE]

  brksi <- c(1, c(which(diff(tin) < 0), length(tin)) + 1)
  brkso <- c(1, c(which(diff(tout) < 0), length(tout)) + 1)

  if (length(brksi) != length(brkso)) stop("Time_in in y does not have same number of jumps as tout.")

  ms <- list(); dms <- list(); ll <- list()
  for (i in seq_len(length(brksi) - 1)) {
    indi <- seq(brksi[i], brksi[i + 1] - 1)
    indo <- seq(brkso[i], brkso[i + 1] - 1)
    res  <- apply(obs[indi, , drop = FALSE], 2, function(y) {
      K0      <- outer(tin[indi], tin[indi], "-")^2
      K       <- exp(-K0 / (2 * p[1]))
      diag(K) <- diag(K) + p[2]
      qrK <- try(qr(K), silent = TRUE)
      if (is(qrK, "qr")) {
        Ky  <- try(solve(qrK, y), silent = TRUE)
        if (length(Ky) != 1) {
          Kout  <- outer(tout[indo], tin[indi], "-")
          m     <- as.vector(exp(- Kout^2 / (2 * p[1])) %*% Ky)
          dm    <- as.vector((-Kout * exp(- Kout^2 / (2 * p[1])) / p[1]) %*% Ky)
        } else {
          # return linear interpolation with derivate averaged between neighbours (simpson type of diff estimate)
          m   <- approx(x = xin, y = y, xout = tout[indo])$y
          dm  <- diff(m) / diff(xout)
          dm  <- c(dm[1], (dm[-1] + dm[-length(dm)]) / 2, dm[length(dm)])
        }
        ll <- log(abs(sum(y * Ky))) + log(abs(prod(diag(qrK$qr))) + 1e-20) / nrow(K)
      } else {
        # return linear interpolation with derivate averaged between neighbours (simpson type of diff estimate)
        m   <- approx(x = xin, y = y, xout = tout[indo])$y
        dm  <- diff(m) / diff(xout)
        dm  <- c(dm[1], (dm[-1] + dm[-length(dm)]) / 2, dm[length(dm)])
        ll <- NA
      }
      return(list(m = m, dm = dm, ll = ll))
    })
    ms[[i]]  <- lapply(res, function(r) r$m)
    ms[[i]]  <- cbind(Time = tout[indo], do.call("cbind", ms[[i]]))
    dms[[i]] <- lapply(res, function(r) r$dm)
    dms[[i]] <- cbind(Time = tout[indo], do.call("cbind", dms[[i]]))
    ll[[i]]  <- do.call(c, lapply(res, function(r) r$ll))
  }

  return(list(ms = do.call("rbind", ms), dms = do.call("rbind", dms), ll = do.call("+", ll)))
}
