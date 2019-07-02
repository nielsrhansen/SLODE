##----------------##
##   Assessment   ##
##----------------##

get_network <- function(makobj, ind) {
  if (all(!ind)) {
    matrix(FALSE, ncol(makobj$A), ncol(makobj$A)) 
  } else {
    (abs(t((makobj$B - makobj$A)[ind, , drop = FALSE])) %*% abs(makobj$A[ind, , drop = FALSE])) != 0
  }
}

# fpr and tpr are equal length vectors, 
auroc <- function(fpr, tpr) {
  o <- order(fpr)
  
  # Order
  fpr <- fpr[o]
  tpr <- tpr[o]
  
  # Sum of rectangular area and triangular area
  ret <- sum((tpr[-length(tpr)] + diff(tpr) / 2) * diff(fpr))
  
  ret 
}


# takes k_hat (matrix) and returns FP and TP of rates
k_hat_2_rate <- function(khat, truth) {
  list(FP = apply(khat != 0 , 2, function(x) {sum(x & (truth$k == 0))}),
    TP = apply(khat != 0 , 2, function(x) {sum(x & (truth$k != 0))}))
}
# khats is list of khat matrix (for one estimator) stab level
stab_2_rate <- function(khats, truth) {
  # Stability part
  stab <- list()
  sup <- lapply(khats, function(khat) {
    apply(khat != 0, 1, any)
  })
  sup <- Reduce("+", sup) / length(sup)
  prob <- seq(min(sup), max(sup), length.out = 100)
  prob <- prob[prob > 0]
  
  stab$FP <- sapply(rev(prob), function(p) {
    sum((sup >= p) & (truth$k == 0))
  })
  stab$TP <- sapply(rev(prob), function(p) {
    sum((sup >= p) & (truth$k != 0))
  })
  
  return(stab)
}
# takes k_hat matrix and reduces chooses a subset in the middle (wrt to number of param)
narrow <- function(k_hat, a, b) {
  stopifnot(a < b)
  nparam <- apply(k_hat != 0, 2, sum)
  k_hat[, which(nparam >= a & nparam <= b), drop = FALSE]
}

k_hat_2_netw <- function(khat, truth) {
  true_network <- get_network(truth$makobj, truth$k != 0)
  
  list(FP = apply(khat, 2, function(k) {
    sum(get_network(truth$makobj, k != 0) & !true_network)
  }),
    TP = apply(khat, 2, function(k) {
      sum(get_network(truth$makobj, k != 0) & true_network)
    }))
}
stab_2_netw <- function(khats, truth) {
  true_network <- get_network(truth$makobj, truth$k != 0)

  # Stability part
  stab <- list()
  sup <- lapply(khats, function(khat) {
    apply(khat != 0, 1, any)
  })
  sup <- Reduce("+", sup) / length(sup)
  prob <- seq(min(sup), max(sup), length.out = 100)
  prob <- prob[prob > 0]
  
  stab$FP <- sapply(rev(prob), function(p) {
    sum(get_network(truth$makobj, sup >= p) & !true_network)
  })
  stab$TP <- sapply(rev(prob), function(p) {
    sum(get_network(truth$makobj, sup >= p) & true_network)
  })

  return(stab)
}

## main function which for given est_return (return from estimators) and truth produces the assessment parameters
assess <- function(est_return, truth, validation_y, stability = TRUE) {
  # structure of est_return:
  # clas and stab
  # each has k_hats and cpu_time
  # k_hat and cpu_time
  # k_hat has three elements (tsars, aim, maker)
  # each of which has a number of matrices (in a list) with hat(k) values (the first one is using full data and thus used in mspe)
  
  
  ## Rate classification
  if (stability) {
    khats_stab <- do.call(function(...) Map(list, ...), est_return$stab$k_hats)
    ret <- list(
      CN = sum(truth$k == 0),
      CP = sum(truth$k != 0),
      clas = lapply(est_return$clas$k_hats, k_hat_2_rate, truth = truth),
      stab = lapply(khats_stab, stab_2_rate, truth = truth),
      stab_narrow = lapply(
        lapply(khats_stab, function(e) lapply(e, narrow, a = floor(truth$d * .5), b = 3 * truth$d)),
        stab_2_rate, truth = truth))
  } else {
    ret <- list(
      CN = sum(truth$k == 0),
      CP = sum(truth$k != 0),
      clas = lapply(est_return$clas$k_hats, k_hat_2_rate, truth = truth))
  }

    
  ## MSPE
  mspe <- lapply(est_return$clas$k_hats, # get first element in each (here full data is used) 
    function(ks) {
      mpe <- apply(ks, 2, function(k) {
        mean((validation_y - numsolve(o = truth$makobj, time = validation_y[, 1], param = list(k = k), x0 = truth$x0))^2)
      })
      sel <- which.min(mpe)
      mean((truth$trajectory[, -1] - numsolve(o = truth$makobj, time = truth$trajectory[, 1], 
        param = list(k = ks[, sel]), x0 = truth$x0)[, -1])^2)
  })
  
  
  ## Selected k for MSPE
  selk <- lapply(est_return$clas$k_hats, # get first element in each (here full data is used) 
    function(ks) {
      mpe <- apply(ks, 2, function(k) {
        mean((validation_y - numsolve(o = truth$makobj, time = validation_y[, 1], param = list(k = k), x0 = truth$x0))^2)
      })
      sel <- which.min(mpe)
      ks[, sel]
    })
  
  
  ## network classification
  true_network <- get_network(truth$makobj, truth$k != 0)
  if (stability) {
    network <- list(
      CN = sum(!true_network),
      CP = sum(true_network),
      clas = lapply(est_return$clas$k_hat, k_hat_2_netw, truth = truth),
      stab = lapply(khats_stab, stab_2_netw, truth = truth),
      stab_narrow = lapply(
        lapply(khats_stab, function(e) lapply(e, narrow, a = floor(truth$d * .5), b = 3 * truth$d)),
        stab_2_netw, truth = truth))    
  } else {
    network <- list(
      CN = sum(!true_network),
      CP = sum(true_network),
      clas = lapply(est_return$clas$k_hat, k_hat_2_netw, truth = truth))
  }

  return(list(rate = ret,
    network = network,
    MSPE = unlist(mspe),  
    CPU = lapply(est_return, getElement, "cpu_time"),
    selk = selk))
}



### Interpret results
# fprs and tprs are lists (methods) of lists (samples) for false pos rates / true pos rates respectively 
get_rocs <- function(fprs, tprs, nsample = 0) {
  stopifnot(length(fprs) == length(tprs))
  stopifnot(all.equal(lapply(fprs, length), lapply(tprs, length)))
  stopifnot(all(unlist(lapply(lapply(fprs, length), ">=", nsample))))
  
  # Smooth
  smo <- Map(function(fp, tp) {
    sx <- seq(0, 1, by = 0.05)
    sys <- Map(function(x, y) {
      approx(x = x, y = y, xout = sx)$y
    }, x = fp, y = tp)
    sy <- do.call(cbind, sys)
    sy <- apply(sy, 1, median, na.rm = TRUE)
    
    data.frame(FPR = c(0, sx, 1), TPR = c(0, sy, 1))
  }, fprs, tprs)
  smo <- Map(function(data, name) {
    data$method <- name
    data
  }, smo, names(smo))
  smo <- do.call(rbind, smo)
  smo$sample <- "smooth"
  
  # Samples
  if (nsample > 0) {
    sam <- Map(function(fpr, tpr) {
      Map(function(fp, tp) {
        data <- data.frame(FPR = fp, TPR = tp)
        data <- unique(data)
        data[order(data$FPR), ]
      }, fpr, tpr)
    }, fprs, tprs) 
    sam <- Map(function(l, name) {
      len <- unlist(lapply(l, nrow))
      l <- l[order(len, decreasing = TRUE)[seq_len(nsample)]]
      l <- Map(function(x, s) {
        x$sample <- s
        x
      }, l, seq_along(l))
      l <- do.call(rbind, l)
      l$method <- name
      l
    }, sam, names(sam))
    sam <- do.call(rbind, sam)
    smo <- rbind(sam, smo)
  }

  smo
}


# recall and precis are lists (methods) of lists (samples) for recalls / precisions respectively 
get_precis_recall <- function(recall, precis, nsample = 0) {
  stopifnot(length(recall) == length(precis))
  stopifnot(all.equal(lapply(recall, length), lapply(precis, length)))
  stopifnot(all(unlist(lapply(lapply(recall, length), ">=", nsample))))
  
  # Smooth
  smo <- Map(function(reca, prec) {
    sx <- seq(0, 1, by = 0.05)
    sys <- Map(function(x, y) {
      approx(x = c(0, x), y = c(0, y), xout = sx)$y
    }, x = reca, y = prec)
    sy <- do.call(cbind, sys)
    sy <- apply(sy, 1, median, na.rm = TRUE)
    
    data.frame(recall = sx, precis = sy)
  }, recall, precis)
  smo <- Map(function(data, name) {
    data$method <- name
    data
  }, smo, names(smo))
  smo <- do.call(rbind, smo)
  smo$sample <- "smooth"
  
  # Samples
  if (nsample > 0) {
    sam <- Map(function(reca, prec) {
      Map(function(rec, pre) {
        data <- data.frame(recall = rec, precis = pre)
        data <- unique(data)
        data[order(data$recall), ]
      }, reca, prec)
    }, recall, precis) 
    sam <- Map(function(l, name) {
      len <- unlist(lapply(l, nrow))
      l <- l[order(len, decreasing = TRUE)[seq_len(nsample)]]
      l <- Map(function(x, s) {
        x$sample <- s
        x
      }, l, seq_along(l))
      l <- do.call(rbind, l)
      l$method <- name
      l
    }, sam, names(sam))
    sam <- do.call(rbind, sam)
    smo <- rbind(sam, smo)
  }

  smo
}
