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
simname <- "sim_01"
source("R/simulation_parameters.R")
source("R/estimators.R")
source("R/assessment.R")
source("R/plotting.R")


## Prepare result names ##
files <- dir(paste0("Results/", simname, "/"))
files <- files[grep("res_", files)]
nums  <- sub("res_", "", files)
nums  <- sub(".RData", "", nums)
nums  <- as.numeric(nums)
nums  <- nums[is.finite(nums)]


rates <- list()
netws <- list()
AUCs  <- list()
mspes <- list()
cpus  <- list()
s <- 1

#i<-nums[9]; i <- 9
for (i in nums) {
  ## Sim parameters
  sim_param <- sim_parameters[i, ]
  
  ## Load results ##
  load(paste0("Results/", simname, "/res_", i, ".RData"))

  ## Remove the errors
  to.res <- to.res[unlist(lapply(to.res, function(x) !is(x, "try-error")))]
  if (length(to.res) == 0) next
  
  ## Rename those methods that repeat
  ren <- function(nn) {
    for (i in 2:length(nn)) {
      if (nn[i] %in% nn[seq_len(i - 1)]) nn[i] <- paste0(nn[i], rep("_", sum(nn[i] == nn[seq_len(i - 1)])))
    }
    nn
  }
  to.res <- lapply(to.res, function(x) {
    names(x$rate$clas) <- ren(names(x$rate$clas))
    names(x$network$clas) <- ren(names(x$network$clas))
    names(x$MSPE) <- ren(names(x$MSPE))
    names(x$CPU$clas) <- ren(names(x$CPU$clas))
    x
  })
  
  ## Cut down stability, so does not report more variables than classic ##
  cutdown <- function(meth, rmax) {
    if (length(meth$FP) == 0) {
      meth$FP <- numeric(0)
      meth$TP <- numeric(0)
      meth
    } else {
      report_meth <- meth$FP + meth$TP
      report_ind  <- (report_meth <= rmax)
      meth <- lapply(meth, function(x) x[report_ind])
      meth
    }
  }
  to.res <- lapply(to.res, function(ret) {
    report_max <- lapply(ret$rate$clas, function(meth) {
      max(c(meth$FP + meth$TP, 2))
    })

    ret
  })

  ##----------------------------------------------------------##
  ##  Smooth and save first ten curves for network and rates  ##
  ##----------------------------------------------------------##

  ## Rates ## 
  precis <- lapply(to.res, function(x) {
    lapply(x$rate[-c(1,2)], 
      function(est) lapply(est, 
        function(meth) {
          if (length(meth$TP) == 0 | length(meth$FP) == 0) {
            c(x$rate$CP / (x$rate$CP + x$rate$CN))
          } else {
            ss <- meth$TP + meth$FP
            c(meth$TP[ss > 0] / ss[ss > 0], x$rate$CP / (x$rate$CP + x$rate$CN)) 
          }
        })) # All also reports full model (we cannot define precision for empty model)
  })
  # now its ordered in three levels: replicate, method, estimator. 
  precis <- do.call(function(...) Map(list, ...), precis)                               # now: method, replicate, estimator
  precis <- lapply(precis, function(meth) do.call(function(...) Map(list, ...), meth))  # now: method, estimator, replicate
  
  recall <- lapply(to.res, function(x) {
    lapply(x$rate[-c(1,2)], 
      function(est) lapply(est, 
        function(meth) {
          if (length(meth$TP) == 0 | length(meth$FP) == 0) {
            c(1)
          } else {
            ss <- meth$TP + meth$FP
            c(meth$TP[ss > 0] / x$rate$CP, 1) 
          }
        })) # All also reports full model 
  })
  # now its ordered in three levels: replicate, method, estimator. We desire: method, estimator, replicate
  recall <- do.call(function(...) Map(list, ...), recall) 
  recall <- lapply(recall, function(meth) do.call(function(...) Map(list, ...), meth))
  
  rates[[i]] <- Map(get_precis_recall, recall, precis)
  rates[[i]] <- lapply(rates[[i]], function(rat) cbind(rat, sim_param))

  
  ## Network ##
  tprs <- lapply(to.res, function(x) {
    lapply(x$network[-c(1,2)], 
      function(est) lapply(est, 
        function(meth) {
          if (length(meth$TP) == 0 | length(meth$FP) == 0) {
            c(0, 1)
          } else {
            c(0, meth$TP / x$network$CP, 1) 
          }
        }))   # all estimators must report empty and full model
  })
  tprs <- do.call(function(...) Map(list, ...), tprs) 
  tprs <- lapply(tprs, function(meth) do.call(function(...) Map(list, ...), meth))
  
  fprs <- lapply(to.res, function(x) {
    lapply(x$network[-c(1,2)], 
      function(est) lapply(est, 
        function(meth) {
          if (x$network$CN > 0 & length(meth$FP) > 0) {
            return(c(0, meth$FP / x$network$CN, 1))                  # all estimators must report empty and full model
          } else {
            return(c(0, rep(1 / sim_param$d^2, length(meth$FP)), 1)) # how to define FPR if no true negatives?
          }}))   
  })
  fprs <- do.call(function(...) Map(list, ...), fprs) 
  fprs <- lapply(fprs, function(meth) do.call(function(...) Map(list, ...), meth))

  netws[[i]] <- Map(get_rocs, fprs, tprs)
  netws[[i]] <- lapply(netws[[i]], function(net) cbind(net, sim_param))

  ## AUC
  AUCs[[i]]  <- Map(function(ff, tt) {
    meth <- Map(function(fff, ttt) {
      unlist(Map(auroc, fpr = fff, tpr = ttt))
    }, fff = ff, ttt = tt)
    meth <- melt(cbind(data.frame(do.call(cbind, meth), sim_param)), id.vars = c("alpha", "d", "k_type", "n", "n_series", "sigma", "w_type"))
    meth
  }, ff = fprs, tt = tprs)
  
  
  
  ## MSPE and CPU Time
  mspe <- lapply(to.res, function(x) x$MSPE)
  mspe <- do.call(cbind, mspe)
  # mspe <- t(t(mspe) / mspe["im1",])
  mspes[[i]] <- data.frame(method = rep(rownames(mspe), ncol(mspe)),
    value = as.vector(mspe))
  mspes[[i]] <- cbind(mspes[[i]], sim_param)
  
  
  cpu <- lapply(to.res, function(x) unlist(x$CPU))
  cpu <- as.data.frame(do.call(rbind, cpu))
  cpu <- melt(cpu)
  cpus[[i]] <- cbind(cpu, sim_param)
  
  
  cat(sprintf("Finished %3d of %d\n", s, length(nums)))
  s <- s + 1
}
#warnings()

rates <- do.call(function(...) Map(list, ...), rates[nums]) # Before: setup, method. Now: method, setup. 
rates <- lapply(rates, function(meth) do.call(rbind, meth)) # Now: method (each entry with combined data.frame over setups and results)
netws <- do.call(function(...) Map(list, ...), netws[nums])
netws <- lapply(netws, function(meth) do.call(rbind, meth))
mspes <- do.call(rbind, mspes)
cpus  <- do.call(rbind, cpus)
AUCs <- do.call(function(...) Map(list, ...), AUCs[nums]) # Before: setup, method. Now: method, setup. 
AUCs <- lapply(AUCs, function(meth) do.call(rbind, meth))



## AUC GRAPH
areas <- subset(AUCs$clas, variable %in% c("im1", "maker2", "scad", "tsars"))
areas$variable <- droplevels(areas$variable)
areas <- ddply(areas, .(alpha, d, k_type, n, n_series, sigma, w_type, variable), summarise,
  med = median(value),
  lo = quantile(value, probs = 0.05), #mean(value) - 1.96 * sd(value) / sqrt(length(value)),
  hi = quantile(value, probs = 0.95)) #mean(value) + 1.96 * sd(value) / sqrt(length(value)))
areas$n_series <- as.factor(areas$n_series)
levels(areas$n_series) <- paste0("E: ", levels(areas$n_series))

areas$alpha <- as.factor(areas$alpha)
levels(areas$alpha) <- paste0("alpha: ", levels(areas$alpha))

areas$sigma <- as.factor(areas$sigma)
levels(areas$sigma) <- paste0("sigma: ", levels(areas$sigma))

areas$d <- as.factor(areas$d)

# rename and reorder levels to lexiographic
levels(areas$variable)
levels(areas$variable) <- c("EGM", "SCAD", "IM", "AIM")
areas$variable <- factor(areas$variable, sort(levels(areas$variable)))

pd <- position_dodge(0.45)
pdf(paste0("../figures/auroc.pdf"), width = 8, height = 10)
gg <- ggplot(areas, aes(x = d, y = med, color = variable)) + 
  geom_point(position = pd, size = 2) +
  # geom_errorbar(position = pd, width = .2) +
  facet_grid(n_series + alpha ~ sigma, labeller = label_parsed) +
  xlab("Number of species (d)") +
  ylab("AUROC") + #coord_cartesian(ylim = c(0.45, .95)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", alpha = .5) +
  ggplot_theme +
  scale_color_manual(name = "Method", values = col_scheme[levels(areas$variable)])
print(gg)
dev.off()




## rates precision-recalls
rates <- lapply(rates, function(rate) {
  rate$n_series  <- factor(rate$n_series, levels = levels(as.factor(sim_parameters$n_series)))
  levels(rate$n_series) <- paste0("E: ", levels(rate$n_series))
  
  rate$sigma     <- factor(rate$sigma, levels = levels(as.factor(sim_parameters$sigma)))
  levels(rate$sigma) <- paste0("sigma: ", levels(rate$sigma))
  
  rate$method    <- factor(rate$method)
  rate$sample    <- factor(rate$sample)
  rate
})


# Extract methods, rename and re-order
rate <- subset(rates$clas, method %in% paste0("", c("im1", "maker2", "scad", "tsars")))
rate$method <- droplevels(rate$method)
levels(rate$method)
levels(rate$method) <- c("IM", "AIM", "SCAD", "EGM")
rate$method <- factor(rate$method, sort(levels(rate$method)))


pdf(paste0("../figures/prec_reca.pdf"), paper = "a4")
# w_type_val <- "a"; k_type_val <- "a"; alpha_val <- 1; d_val <- 9
for (w_type_val in unique(rate$w_type)) {
  for (k_type_val in unique(rate$k_type)) {
    for (alpha_val in unique(rate$alpha)) {
      for (d_val in unique(rate$d)) {
        spec <- subset(rate, 
          d == d_val & alpha == alpha_val & k_type == k_type_val & w_type == w_type_val)
        
        gg <- ggplot(spec, aes(x = recall, y = precis, color = method))
        
        gg <- gg +
          geom_path() +
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
          facet_grid(n_series ~ sigma, drop = FALSE, labeller = label_parsed) +
          ggtitle(substitute(paste("d", "=", d_val, " and ", alpha, "=", alpha_val), list(d_val = d_val, alpha_val = alpha_val))) +
          xlab("Recall") + ylab("Precision") + ggplot_theme +
          scale_x_continuous(breaks = c(0, 0.5, 1)) +
          scale_color_manual(name = "Method", values = col_scheme[levels(spec$method)])
        
        print(gg)
        
        if (alpha_val == 1 & d_val == 9) {
          rates_gg <- gg
        }
      }
    }
  }
}
dev.off()

pdf(paste0("../figures/prec_reca_single.pdf"), height = 8, width = 8)
print(rates_gg + ggtitle("")) 
dev.off()




## network roc
netws <- lapply(netws, function(rate) {
  rate$n_series  <- factor(rate$n_series, levels = levels(as.factor(sim_parameters$n_series)))
  levels(rate$n_series) <- paste0("E: ", levels(rate$n_series))
  
  rate$sigma     <- factor(rate$sigma, levels = levels(as.factor(sim_parameters$sigma)))
  levels(rate$sigma) <- paste0("sigma: ", levels(rate$sigma))
  
  rate$method    <- factor(rate$method)
  rate$sample    <- factor(rate$sample)
  rate
})


# Extract methods, rename and re-order
netw <- subset(netws$clas, method %in% paste0("", c("im1", "maker2", "scad", "tsars")))
netw$method <- droplevels(netw$method)
levels(netw$method)
levels(netw$method) <- c("IM", "AIM", "SCAD", "EGM")
netw$method <- factor(netw$method, sort(levels(netw$method)))

pdf(paste0("../figures/netw_rocs.pdf"), paper = "a4")
# w_type_val <- "a"; k_type_val <- "a"; alpha_val <- 1; d_val <- 9
for (w_type_val in unique(netw$w_type)) {
  for (k_type_val in unique(netw$k_type)) {
    for (alpha_val in unique(netw$alpha)) {
      for (d_val in unique(netw$d)) {
        spec <- subset(netw, 
          d == d_val & alpha == alpha_val & k_type == k_type_val & w_type == w_type_val)
        
        gg <- ggplot(spec, aes(x = FPR, y = TPR, color = method))
        
        gg <- gg +
          geom_path() +
          geom_abline(intercept = 0, slope = 1) +
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
          facet_grid(n_series ~ sigma, drop = FALSE, labeller = label_parsed) +
          ggtitle(substitute(paste("d", "=", d_val, " and ", alpha, "=", alpha_val), list(d_val = d_val, alpha_val = alpha_val))) + 
          ggplot_theme +
          xlab("False Positive Rate") + ylab("True Positive Rate") + 
          scale_x_continuous(breaks = c(0, 0.5, 1)) +
          scale_color_manual(name = "Method", values = col_scheme[levels(spec$method)])
        
        print(gg)
        
        if (alpha_val == 1 & d_val == 9) {
          netw_gg <- gg
        }
      }
    }
  }
}
dev.off()

pdf(paste0("../figures/netw_rocs_single.pdf"), height = 8, width = 8)
print(netw_gg + ggtitle("")) 
dev.off()





## mspe violin
# mspes <- mspes_backup
mspes_backup <- mspes
mspes$n_series  <- factor(mspes$n_series, levels = levels(as.factor(sim_parameters$n_series)))
mspes$sigma     <- factor(mspes$sigma, levels = levels(as.factor(sim_parameters$sigma)))
mspes$method    <- factor(mspes$method)
mspes_old <- mspes
mspes <- subset(mspes, 
  method %in% c("im1", "maker2", "scad", "tsars"))
mspes$method <- droplevels(mspes$method)

# Levels
levels(mspes$method) 
levels(mspes$method) <- c("IM", "AIM", "SCAD", "EGM")
mspes$method <- factor(mspes$method, levels = sort(levels(mspes$method)))

mspes <- ddply(mspes, .(alpha, d, k_type, n, n_series, sigma, w_type, method), summarise,
  med = median(value),
  lo = quantile(value, probs = 0.05),
  hi = quantile(value, probs = 0.95))

mspes$alpha <- as.factor(mspes$alpha)
levels(mspes$alpha) <- paste0("alpha: ", levels(mspes$alpha))

levels(mspes$sigma) <- paste0("sigma: ", levels(mspes$sigma))

mspes$d <- as.factor(mspes$d)
# levels(mspes$d) <- paste0("d: ", levels(mspes$d))


pd <- position_dodge(0.45)
pdf(paste0("../figures/mspes.pdf"), width = 8, height = 5)
gg <- ggplot(subset(mspes, n_series == 4), aes(x = d, y = med, color = method)) + 
  geom_point(position = pd, size = 2) +
  # geom_errorbar(mapping = aes(ymin = lo, ymax = hi), position = pd, width = .2) +
  facet_grid(alpha ~ sigma, labeller = label_parsed) +
  xlab("Number of species (d)") +
  ylab("MSE") + coord_cartesian(ylim = c(0, 55)) + 
  # geom_hline(yintercept = 1, linetype = "dotted", alpha = .5) +
  scale_color_manual(name = "Method", values = col_scheme[levels(mspes$method)]) + ggplot_theme
print(gg)
dev.off()



## CPU 
# cpus <- cpus_backup
cpus_backup <- cpus
cpus$n_series  <- factor(cpus$n_series, levels = levels(as.factor(sim_parameters$n_series)))
cpus$sigma     <- factor(cpus$sigma, levels = levels(as.factor(sim_parameters$sigma)))
cpus$method    <- factor(cpus$variable)

cpus <- subset(cpus, 
  method %in% c(paste0("clas.", c("im1", "maker2", "scad", "tsars"))))
cpus$method <- droplevels(cpus$method)

# Method levels
levels(cpus$method) 
levels(cpus$method) <- c("EGM", "SCAD", "IM", "AIM")
cpus$method <- factor(cpus$method, levels = sort(levels(cpus$method)))

cpus <- ddply(cpus, .(alpha, d, k_type, n, n_series, sigma, w_type, method), summarise,
  med = median(value),
  lo = quantile(value, probs = 0.05),
  hi = quantile(value, probs = 0.95))

cpus$alpha <- as.factor(cpus$alpha)
levels(cpus$alpha) <- paste0("alpha: ", levels(cpus$alpha))

levels(cpus$sigma) <- paste0("sigma: ", levels(cpus$sigma))

levels(cpus$n_series) <- paste0("E: ", levels(cpus$n_series))

pd <- position_dodge(0.45)
pdf(paste0("../figures/cpus.pdf"), width = 8, height = 10)
gg <- ggplot(cpus, aes(x = d, y = med, color = method, ymin = lo, ymax = hi)) + 
  geom_point(position = pd) + geom_line(position = pd) +
  geom_errorbar(position = pd, width = .2) +
  facet_grid(alpha + sigma ~ n_series, labeller = label_parsed) +
  xlab("d") +
  ylab("CPU Time [s]") + 
  scale_color_manual(name = "Method", values = col_scheme[levels(cpus$method)]) + 
  ggplot_theme + scale_y_log10() +
  scale_x_continuous(breaks = c(7, 9, 11))
print(gg)
dev.off()
