## estimator for prior ##

aimer <- function(ys, A, CC, sc1, sc2) {
  y <- do.call(rbind, ys)
  d <- ncol(A)
  a_subs <- list()
  
  for (s in seq_len(d)) {
    A_sub <- A[ A[, s] != 0, , drop = FALSE]
    ratt_sub <- ratmak(A = A_sub, C = CC, s = solver(step_max = 100000, h_init = 0.001), 
      r1 = reg(
        contexts = apply(sc1, 2, function(sc) matrix(sc, nrow = nrow(CC))[, A[, s] != 0 , drop = FALSE])),
      r2 = reg(
        contexts = apply(sc2, 2, function(sc) matrix(sc, nrow = nrow(CC))[, A[, s] != 0 , drop = FALSE])))
    a_subs[[s]] <- episode::aim(ratt_sub, opt(y = y), adapts = NULL)
    cat("s:", s, "\n")
  }
  
  netw_gs <- lapply(a_subs, function(a_) {
    netw_g <- sapply(seq_len(ncol(a_$params$theta1)), function(i) {
      field(a_$o, x = abs(a_$x0s[seq_len(d), 1]), 
        param = list(theta1 = a_$params$theta1[, i],
          theta2 = a_$params$theta2[, i]), differentials = TRUE)$f_dx != 0 
    })
    netw_g
  })
  
  netw_gs
}

# estimates network by parents
est <- function(ys, A, CC, sc1, sc2, stab_remove) {
  # classic
  aimer_clas <- aimer(ys, A = A, CC = CC, sc1 = sc1, sc2 = sc2)
  
  # with stability select
  yss <- apply(stab_remove[,seq_along(ys),, drop = FALSE], 3, function(i) {
    Map(function(y, remi) {
      y[-remi, , drop = FALSE]
    }, y = ys, remi = split(t(i), seq_len(ncol(i))))
  })
  s <- 0
  aimer_stab <- lapply(yss, function(y) {print(s <<- s + 1); aimer(y, A = A, CC = CC, sc1 = sc1, sc2 = sc2)})
}

no_neighbours <- function(n, l) {
  # n is number of indeces to choose from, l is number of chosen, no neighbours allowed
  while (TRUE) {
    proposal <- sort(sample(n, size = l))
    if (all(diff(proposal) > 1)) break
  }
  proposal
}
