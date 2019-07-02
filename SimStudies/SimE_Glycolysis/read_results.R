##----------------##
##  Read results  ##
##----------------##

rm(list = ls())


## Libraries ##
library(tsars)
library(episode)
library(magrittr)
library(grid)
library(gridExtra)
library(plyr)
library(ggplot2)
library(reshape)


## Sources ##
simname <- "simprior_3"
# source("R/est.R")
source("R/assess.R")
source("R/Hynne_model.R")
load("trajs_inter.RData")

rat <- ratmak(A = A, C = CC, s = solver(step_max = 100000, h_init = 0.001))
netw <- field(o = rat, x = x0, param = list(theta1 = K1, theta2 = K2), differentials = TRUE)
netw <- netw$f_dx != 0
source("../ggplot_theme.R")


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


## Prepare result names ##
files <- dir(paste0("Results/", simname, "/"))
files <- files[grep("res_", files)]
nums  <- sub("res_", "", files)
nums  <- sub(".RData", "", nums)
nums  <- as.numeric(nums)
nums  <- nums[is.finite(nums)]
netws <- list()

for (i in nums) {
  ## Sim parameters
  sim_param <- sim_parameters[i, ]
  thrs <- seq(0, 1, length.out = 101)
  
  ## Load results ##
  load(paste0("Results/", simname, "/res_", i, ".RData"))
  
  if (sim_param$a == "no_prior") {
    # collect nonzero params
    nonz <- rep(0, nrow(to.res[[1]]$params))
    for (ss in seq_len(20)) {
      E <- Es[ss]
      nonz <- nonz + rowSums(to.res[[E]]$params != 0)
      if (ss %% 5 == 0) {
        netw_gs[[(ss %/% 5)]] <- list(nonz=nonz / (max(nonz) - min(nonz)), odeobj=to.res[[E]]$odeobj) 
      }
    }
    fp2 = list()
    ss <-0
    #threshold and get netw est
    for (ests in netw_gs) {
      ss <<- ss + 1
      fp2[[ss]] <- cbind(1, sapply(thrs, function(thr) {
        p <- ests$nonz
        p[p < thr] = 0
        nt <- field(ests$odeobj, param=p, x=x0, differentials = TRUE)$f_dx
        c(thr=thr, TPR = sum(nt & netw) / sum(netw), FPR = sum(nt & !netw) / sum(!netw))
      }), 0)
    }
  } else {
    netw_g <- Matrix(0, ncol = rat$d, nrow = rat$d)
    netw_gs <- list()
    for (ss in seq_len(20)) {
      E <- Es[ss] # context number, or which is inhibited
      netw_g[-E, -E] <- netw_g[-E, -E] + Reduce("+", to.res[[E]])
      if (ss %% 5 == 0) {
        netw_gs[[(ss %/% 5)]] <- apply(netw_g, 2, function(x) {
          if (sum(x^2) > 0) {
            (x - min(x)) / (max(x) - min(x))
          } else {
            x
          }
        })
      }
    }
    # Threshold them
    thrs <- seq(0, 1, length.out = 101)
    fp2 <- sapply(rev(thrs), function(thr) {
      lapply(netw_gs, function(nt) {
        nt <- nt > thr
        c(thr=thr, TPR = sum(nt & netw) / sum(netw), FPR = sum(nt & !netw) / sum(!netw))
      })
    })
    fp2 <- apply(fp2, 1, function(l) list(cbind(0, do.call(cbind, l), 1)))
    fp2 <- lapply(fp2, function(x) x[[1]])
  }
  
  # add params
  E <- 0
  fp2 <- lapply(fp2, function(fp) {
    E <<- E + 5
    data.frame(sim_param, 
               t(fp), E = E)
  })
  netws[[i]] <- do.call(rbind, fp2)
  
  aurocs <- lapply(fp2, function(fp) {
    sum(diff(fp$FPR) * fp$TPR[-length(fp$TPR)]) + 
      sum(diff(fp$TPR) * diff(fp$FPR)) / 2
  })
  
  cat("Setting:", i, "finished. AUROCS:", abs(signif(unlist(aurocs), 2)), " ", unlist(sim_param),  "\n")
  save(netws, sim_parameters, file = paste0("Results/", simname, "/cleaned", ".RData"))
}
# warnings()

# save(netws, sim_parameters, file = paste0("Results/", simname, "/cleaned", ".RData"))
# load(paste0("Results/", simname, "/cleaned", ".RData"))

netws <- do.call(rbind, netws)
netws$thr <- NULL

# Linear interpolate TP to match a fixed grid on FP
xs <- seq(0, 1, by = 0.01)
netws <- ddply(netws, .(sigma, a, eps, E), summarise,
  fpr = xs,
  tpr = approx(x = FPR, y = TPR, xout = xs)$y)
# Mean over replicates
netws <- ddply(netws, .(sigma, a, E, fpr), summarise,
  tpr = mean(tpr))

# aurocs
aurocs <- ddply(netws, .(sigma, a, E), summarise,
  auroc = sum(tpr[-length(tpr)] * diff(fpr)) + sum(diff(fpr) * diff(tpr)) / 2)

netws$a <- as.factor(netws$a)
netws$E <- as.factor(netws$E)
if (!any(grepl("E:", levels(netws$E)))) {
  levels(netws$E) <- paste0("E: ", levels(netws$E))
}
netws$sigma <- as.factor(netws$sigma)
if (!any(grepl("sigma:", levels(netws$sigma)))) {
  levels(netws$sigma) <- paste0("sigma: ", levels(netws$sigma))
}


pdf(paste0("../figures/rocs.pdf"), height = 8, width = 8.25 * 8 / 10)
ggplot(netws, 
  aes(x = fpr, y = tpr, linetype = a)) + geom_line(color = col_scheme["AIM"]) +
  geom_abline(slope = 1, intercept = 0) + 
  xlab("False Positive Rate") + ylab("True Positive Rate") +
  facet_grid(E ~ sigma, labeller = label_parsed) + angle_theme +
  scale_color_discrete(guide = FALSE) +
  scale_linetype_manual(values = c("solid", "twodash", "dotted"),
    labels = c("46 complexes  ", "92 complexes  ", "approx. model"),
    name = "")
dev.off()




## Network Figure ##
library(igraph)
diag(netw) <- FALSE
gg <- graph_from_adjacency_matrix(netw)
nums <- which(sim_parameters$sigma == 0.25 & sim_parameters$a == "extra")
ntws <- list()
s <- 0
for (i in nums) {
  ## Sim parameters
  sim_param <- sim_parameters[i, ]
  
  ## Load results ##
  load(paste0("Results/", simname, "/res_", i, ".RData"))
  
  ## Summarise networks
  netw_g <- Matrix(0, ncol = rat$d, nrow = rat$d)
  for (ss in seq_len(20)) {
    E <- Es[ss] # context number, or which is inhibited
    netw_g[-E, -E] <- netw_g[-E, -E] + Reduce("+", to.res[[E]])
  }
  
  s <<- s + 1
  ntws[[s]] <- netw_g
  
  cat("Setting:", i, "finished. \n")
}


# save(ntws, file = "ntws.RData")
# load("ntws.RData")

nt2 <- lapply(ntws, function(x) {
  apply(x, 2, function(z) {
    sz <- sum(z)
    if (sz > 0) z <- z / sz
    z
  })
})
nt2 <- Reduce("+", nt2)


## Get the most popular edges for each component
tt <- 2
i <- 0
m <- Matrix(apply(nt2, 2, function(x) {
  i <<- i + 1
  z <- logical(length(x))
  x[i] <- -1 # self edges
  zo <- rev(order(x))
  z[zo[seq_len(tt)]] <- TRUE
  z
}))
rownames(m) <- colnames(m) <- rownames(netw)



igraph_options(edge.arrow.size = .5, vertex.size = 10, vertex.label.cex = .5)
igraph_options(edge.arrow.size = .5, vertex.size = 15, vertex.label.cex = .7)

pdf(paste0("../figures/graph_true.pdf"), height = 7, width = 7)
plot(gg, layout = layout_in_circle(t(gg)),
  vertex.shape = rep("rectangle", 22),
  vertex.color = "white")
dev.off()

gg_guess <- graph_from_adjacency_matrix(t(m))
E(gg_guess)$weight <- attr(E(gg_guess), "vnames") %in% attr(E(gg), "vnames")
E(gg_guess)$color[E(gg_guess)$weight == 1] <- 'green'
E(gg_guess)$color[E(gg_guess)$weight == 0] <- 'gray'
pdf(paste0("../figures/graph_2.pdf"), height = 7, width = 7)
plot(gg_guess, layout = layout_in_circle(gg),
  vertex.shape = rep("rectangle", 22),
  vertex.color = "white")
dev.off()





