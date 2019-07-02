## Lotka-Volterra Example ##

rm(list = ls())

source("Sources/Smoother.R")
source("Sources/Main.R")
source("Sources/Estimator.R")
source("Sources/GroupCorD.R")
source("Sources/Recoverer.R")
source("Sources/Regressor.R")
source("Sources/Lineargen.R")
source("Sources/Generator.R")
source("Sources/Wrapper.R")
source("Sources/GroupGrp.R")
library(splines,MASS)
library(grplasso) # install.packages('grplasso');
library(locfit) # install.packages('locfit');
library(fda) # install.packages('fda');
library(episode)
library(ggplot2)
library(reshape)

source("../ggplot_theme.R")

RandomSeed = 1


## Simulation parameters
sim_parameters <- expand.grid(
  v = c(1, 3, 5, 7),
  n_series = c(2, 4, 8),
  sigma = c(.5, 1, 2),
  stringsAsFactors = FALSE)
n <- 51

## Generate Noiseless Data ##
k <- 5; d <- 2 * k
A <- rbind(diag(d), kronecker(diag(k), t(c(1, 1))))
B <- rbind(diag(rep(c(2, 0), k)), kronecker(diag(k), t(c(0, 2))))
m <- episode::mak(A = A, B = B, rx0 = reg(lower = -Inf))
netw_true <- Matrix(t(abs(B - A)) %*% abs(A)) != 0


## Larger MAK
AB <- sapply(1:d, function(i) {
  js <- setdiff(1:d, c(i, ifelse(i %% 2 == 0, i - 1, i + 1)))
  A <- sapply(js, function(j) {
    A <- numeric(d)
    A[c(i, j)] <- 1
    A
  }, simplify = FALSE)
  B <- sapply(js, function(j) {
    B <- numeric(d)
    B[i] <- 2     
    B
  }, simplify = FALSE)
  return(list(A = do.call(rbind, A), B = do.call(rbind, B)))
}, simplify = FALSE)
A_full <- rbind(A, diag(d), 
  do.call(rbind, lapply(AB, getElement, "A")))
B_full <- rbind(B, diag(rep(c(0, 2), k)),
  do.call(rbind, lapply(AB, getElement, "B")))
m_full <- mak(A_full, B_full)


# Draw noise
Bs <- 100
eps <- array(rnorm(d * n * 8 * Bs), dim = c(8 * n, d, Bs))


set.seed(RandomSeed)
time_seq <- seq(0, 5, length.out = n)
x0s <- matrix(runif(d * 8, max = 4), ncol = 8)

traj <- sapply(unique(sim_parameters$v), 
  function(v) numsolve(o = m, time = rep(time_seq, 8), param = c(rep(2, d), rep(v, k)), x0 = x0s), simplify = FALSE)

# Plot of data
traj_e <- traj[[4]]; traj_e[, 1] <- traj_e[, 1] + cumsum(c(0, - pmin(0, diff(traj_e[, 1]))))
ggplot(melt(data.frame(traj_e), id.vars = "Time"), aes(x = Time, y = value, color = variable)) +
  geom_line()


results <- list()
i <- 0
# eps <- eps[,,1:3]

for (sigma in unique(sim_parameters$sigma)) {
  for (n_series in unique(sim_parameters$n_series)) {
    for (tt in seq_along(traj)) {
      i <- i + 1
      
      y_raw <- traj[[tt]]
      
      ss <- 0
      results[[i]] <- apply(eps, 3, function(e) {
        y <- y_raw
        y[, -1] <- y[, -1] + sigma * e
        y <- y[cumsum(c(0, diff(y[, 1]) < 0)) < n_series, ]
        
        
        ### GRADE ###
        obs <- lapply(split(y, cumsum(c(0, diff(y[, 1]) < 0))), matrix, ncol = d + 1) 
        times_e <- seq(from = 0, to = max(time_seq), by = 0.0001)
        smthed <- smoothX(observations = obs, times_e = times_e, deg = deg, h = h, type_smooth = "smoothspline", type_data = "Perturbations")
        xhat <- smthed$Xhat
        
        lambda_N <- 30
        lambda_range <- c(-4,0)
        lambdas <- exp(seq(from = lambda_range[2], to = lambda_range[1], length.out = lambda_N)) 
        lambdas.int <- lambdas * 3
        
        # Fit GRADE
        fits <- GRADE(observations = obs, times_e = times_e, type_reg = "Y", type_data = "Perturbations", 
          type_basis = "monomial", type_smooth = "smoothspline", xhat = xhat, xprime = NULL, lambdas = lambdas, nbasis = 3)
        recint <- countGraph(graph_est = fits$neighbour.path, graph_true = as.matrix(netw_true), self = TRUE)  
        TFP_grade <- cbind(FPR = (recint$NP_sum  - recint$TP_sum) / sum(!netw_true), TPR = recint$TP_sum / sum(netw_true))
        
        
  
        ### AIM ###
        a <- aim(m_full, op = opt(y))
        rod <- rodeo(a)
        
        netw_guess  <- apply(a$params[[1]] != 0, 2, function(ind) {
          A_guess <- A_full[ind, , drop = FALSE] 
          B_guess <- B_full[ind, , drop = FALSE] 
          Matrix(t(abs(B_guess - A_guess)) %*% abs(A_guess)) != 0
        })
        TFP <- lapply(netw_guess, function(n) {
          c(FPR = sum(n & !netw_true), TPR = sum(n & netw_true))
        })
        TFP <- do.call(rbind, TFP)
        TFP <- t(t(TFP) / c(sum(!netw_true), sum(netw_true)))
        TFP_aim <- TFP
        
        
        
        print(rep(i, 100))
        print(ss <<- ss + 1)
        
        return(list(aim = TFP_aim, grade = TFP_grade))
      })
      save(results, file = paste0("results", i, ".RData"))
    }
  }
}
save(results, file = "results.RData")

setting <- list()
i <- 0
for (sigma in unique(sim_parameters$sigma)) {
  for (n_series in unique(sim_parameters$n_series)) {
    for (v in unique(sim_parameters$v)) {
      i <- i + 1
      setting[[i]] <- data.frame(sigma = sigma, E = n_series, v = v)
    }
  }
}


load("results.RData")
lip <- function(x, y, xout) {
  xo <- order(x)
  approx(x = c(0, x[xo], 1), y = c(0, y[xo], 1), xout = xout)$y
}

xs <- seq(0, 1, by = 0.01)
ave_res <- lapply(results, function(sim) {
  sim_res <- lapply(sim, function(x) {
    lapply(x, function(meth) {
      lip(meth[, "FPR"], meth[, "TPR"], xout = xs)
    })
  })
  sim_res <- do.call(function(...) Map(list, ...), sim_res)
  lapply(sim_res, function(s) {
    cbind(FPR = xs, TPR = Reduce("+", s) / length(s))
  })
})

med_res <- Map(function(sim, set) {
  sim_res <- lapply(sim, function(x) {
    lapply(x, function(meth) {
      lip(meth[, "FPR"], meth[, "TPR"], xout = xs)
    })
  })
  sim_res <- do.call(function(...) Map(list, ...), sim_res)
  lapply(sim_res, function(s) {
    bound <- do.call(cbind, s)
    m <- cbind(FPR = c(0, xs, 1), 
      TPR_lo = c(0, apply(bound, 1, quantile, probs = c(0.05)), 1),
      TPR_med = c(0, apply(bound, 1, quantile, probs = c(0.5)), 1),
      TPR_hi = c(0, apply(bound, 1, quantile, probs = c(0.95)), 1))
    cbind(data.frame(m), set)
  })
}, sim = results, set = setting)


med_res_all <- do.call(function(...) Map(list, ...), med_res)
med_res_all <- lapply(med_res_all, function(x) do.call(rbind, x))
med_res_all <- rbind(cbind(med_res_all$aim, Method = "AIM"), cbind(med_res_all$grade, Method = "GRADE"))


med_res_all$sigma <- as.factor(med_res_all$sigma)
levels(med_res_all$sigma) <- paste0("sigma: ", levels(med_res_all$sigma))

med_res_all$v <- as.factor(med_res_all$v)


pdf(paste0("../figures/grade_rocs.pdf"), paper = "a4")
for (vv in levels(med_res_all$v)) {
  gg <- ggplot(subset(med_res_all, v == vv), aes(x = FPR, y = TPR_med, color = Method)) + geom_line() +
    geom_line(aes(x = FPR, y = TPR_lo, color = Method), linetype = "longdash", alpha = .5) +
    geom_line(aes(x = FPR, y = TPR_hi, color = Method), linetype = "longdash", alpha = .5) +
    facet_grid(E ~ sigma, labeller = labeller(sigma = label_parsed, E = label_both)) +
    ylab(paste0("True Positive Rate (v = ", vv, ")")) + xlab("False Positive Rate") +
    scale_color_manual(values = c("#0072B2", "#D55E00")) +
    angle_theme +
    geom_abline(slope = 1, intercept = 0)
  print(gg)
  if (vv == "5") {
    gg_save <- gg
  }
}
dev.off()

pdf(paste0("../figures/grade_rocs_single.pdf"), height = 8, width = 8) 
print(gg_save)
dev.off()



