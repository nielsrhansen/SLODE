# test IM
rm(list = ls())

source("../ggplot_theme.R")

library(episode)
library(ggplot2)
library(reshape)
library(plyr)

# ground truth
d <- 4
A <- matrix(
  c(1, 0, 0, 1,
    0, 1, 0, 0,
    0, 1, 0, 0), ncol = 4, byrow = TRUE)
B <- matrix(
  c(0, 1, 0, 0,
    1, 0, 0, 1,
    1, 0, 1, 0), ncol = 4, byrow = TRUE)

x0 <- c(10, 2, 2, 10)
k <- c(1, 2, 0.75)

m <- mak(A, B, r = reg(reg_type = "scad"))
time <- seq(0, 1, by = .0001)
traj <- numsolve(o = m, time, x0 = x0, param = list(k = k))


sub_traj <- traj[1 + 0:100 * 100, ]

ggplot(melt(data.frame(traj), id.vars = "Time"), aes(x = Time, y = value, color = variable)) + geom_line()


bwds <- seq(0, .6, by = .1)
smoother <- function(y_traj) {
  sapply(bwds, function(bwd) {
    if (bwd <= 0) {
      # linear interpolation
      res <- apply(y_traj[, -1], 2, approx, x = y_traj[, 1], n = max(101, nrow(y_traj)))
    } else {
      res <- apply(y_traj[, -1], 2, ksmooth, x = y_traj[, 1], bandwidth = bwd, kernel = "normal", n.points = max(101, nrow(y_traj)))
    }
    # reorganise
    cbind(Time = res[[1]]$x, do.call(cbind, lapply(res, getElement, "y")))
  }, simplify = FALSE)
}


dattner <- function(x_smooth, y = x_smooth) {
  # Design
  des <- imd(m, op = opt(y), x = x_smooth)
  
  # Accumulate
  X <- apply(matrix(t(des$X[[1]]), nrow = length(A)), 1, cumsum)
  X <- matrix(t(X), ncol = nrow(A), byrow = TRUE)
  Y <- apply(matrix(des$Y, nrow = d), 1, cumsum)
  Y <- as.vector(t(Y))
  
  lm.fit(x = X, y = Y)$coefficients
}


MLE <- function(y_traj) {
  rod <- rodeo(x = m, op = opt(y_traj, nlambda = 5, lambda_min_ratio = 0.00001), 
    x0 = traj[1, -1], params = list(k = rep(0, 3)))
  setNames(c(rod$params$k[, 5], sum(rod$jerr[2, ])), c(paste0("x", 1:3), "jerr"))
}



# Simulation setup
BB <- 250
sigmas <- c(.1, .5, 1)
ns <- c(10, 25, 100)

# Noise
set.seed(270990)
eps <- array(rnorm(max(ns + 1) * d * BB), dim = c(max(ns + 1), d, BB))


## IM simulation
res <- list()
res_MLE <- list()
for (i in seq_along(sigmas)) {
  s <- 0
  res[[i]] <- apply(eps, 3, function(epsilon) {
    cat("IM ", i, " ", s <<- s + 1, "\n")
    # data
    y <- sub_traj
    y[, -1] <- y[, -1] + sigmas[i] * epsilon
    
    lmfit_n10 <- cbind(data.frame(do.call(rbind, lapply(smoother(y[1 + 0:10 * 10, ]), dattner))), bwd = bwds, n = 10)
    lmfit_n25 <- cbind(data.frame(do.call(rbind, lapply(smoother(y[1 + 0:25 * 4, ]), dattner))), bwd = bwds, n = 25)
    lmfit_n100 <- cbind(data.frame(do.call(rbind, lapply(smoother(y), dattner))), bwd = bwds, n = 100)
    
    rbind(lmfit_n10, lmfit_n25, lmfit_n100)
  })
  
  s <- 0
  res_MLE[[i]] <- apply(eps, 3, function(epsilon) {
    cat("MLE ", i, " ", s <<- s + 1, "\n")
    # data
    y <- sub_traj
    y[, -1] <- y[, -1] + sigmas[i] * epsilon

    mle_n10 <- cbind(as.data.frame(matrix(MLE(y[1 + 0:10 * 10, ]), nrow = 1)), n = 10)
    mle_n25 <- cbind(as.data.frame(matrix(MLE(y[1 + 0:25 * 4, ]), nrow = 1)), n = 25)
    mle_n100 <- cbind(as.data.frame(matrix(MLE(y), nrow = 1)), n = 100)
    
    rbind(mle_n10, mle_n25, mle_n100)
  })
  
}


## Reorganise
sig <- 0
IM <- do.call(rbind, 
  lapply(res, function(rr) {
    b <- 0
    rr <- do.call(rbind,
      lapply(rr, function(r) {
        b <<- b + 1
        cbind(r, b = b)
      }))
    sig <<- sig + 1
    cbind(rr, sigma = sigmas[sig])
  }))

sig <- 0
MLES <- do.call(rbind, 
  lapply(res_MLE, function(rr) {
    b <- 0
    rr <- do.call(rbind,
      lapply(rr, function(r) {
        b <<- b + 1
        cbind(r, b = b)
      }))
    sig <<- sig + 1
    cbind(rr, sigma = sigmas[sig])
  }))
names(MLES) <- c(paste0("x", 1:3), "jerr", "n", "b", "sigma")
save(IM, MLES, file = "IMdata.RData")
load("IMdata.RData")

# Summarise over replicates
IM <- melt(IM, id.vars = c("bwd", "n", "b", "sigma"))
IM <- ddply(IM, .(bwd, n, sigma, variable), summarise, 
  median = median(value),
  lo = quantile(value, probs = .05),
  hi = quantile(value, probs = .95))

MLES <- subset(MLES, jerr == 0)
MLES <- MLES[, setdiff(names(MLES), "jerr")]
MLES <- melt(MLES, id.vars = c("n", "b", "sigma"))
MLES <- ddply(MLES, .(n, sigma, variable), summarise, 
  median = median(value),
  lo = quantile(value, probs = .05),
  hi = quantile(value, probs = .95))

# Factor stuff
IM$n <- as.factor(IM$n)
MLES$n <- as.factor(MLES$n)

IM$sigma <- as.factor(IM$sigma)
levels(IM$sigma) <- paste0("sigma: ", levels(IM$sigma))
MLES$sigma <- as.factor(MLES$sigma)
levels(MLES$sigma) <- paste0("sigma: ", levels(MLES$sigma))

IM$variable <- as.factor(IM$variable)
levels(IM$variable) <- c("k[f]", "k[r]", "k[cat]")
MLES$variable <- as.factor(MLES$variable)
levels(MLES$variable) <- c("k[f]", "k[r]", "k[cat]")


MLES$bwd <- -0.1
pd <- position_dodge(0.05)

pdf(file = "../figures/im.pdf", width = 8, height = 8)
gg <- ggplot(IM, aes(x = bwd, y = median, color = n)) + geom_point(position = pd) + geom_line(position = pd) +
  facet_grid(variable ~ sigma,
    labeller = labeller(sigma = label_parsed, variable = label_parsed),
    scales = "free_y") + ylab("") +
  ggplot_theme + 
  xlab("Bandwidth") +
  geom_errorbar(mapping = aes(x = bwd, ymin = lo, ymax = hi, color = n), width = .1, position = pd) + 
  geom_hline(data = data.frame(value = k, variable = as.factor(levels(IM$variable))),
    mapping = aes(yintercept = value)) +
  geom_errorbar(data = MLES, mapping = aes(x = bwd, ymin = lo, ymax = hi, color = n), 
    width = .1, position = pd, linetype = "longdash") +
  geom_point(data = MLES, mapping = aes(x = bwd, y = median, color = n), position = pd)
print(gg)
dev.off()



