# Libraries
library(episode)
library(reshape)
library(ggplot2)
library(glmnet)

# Source
source("R/Hynne_model.R")


set.seed(2709)

# rational mass action kinetics
rat <- ratmak(A = A, C = CC, s = solver(step_max = 100000, h_init = 0.001))
d <- length(x0)

sc1 <- matrix(1, nrow = length(K1), ncol = d)
sc2 <- matrix(1, nrow = length(K2), ncol = d)
for (target in seq_len(d)) {
  inh_react <- CC[, target] != 0
  sc1[, target] <- as.vector(1 - outer(inh_react, rep(1, ncol(K1))))
  sc2[, target] <- as.vector(1 - outer(inh_react, rep(1, ncol(K2))))
}

rat <- ratmak(A = A, C = CC, s = solver(step_max = 100000, h_init = 0.001))
netw <- field(o = rat, x = x0, param = list(theta1 = K1, theta2 = K2), differentials = TRUE)
netw <- netw$f_dx != 0

x0s <- replicate(ncol(sc1), setNames(x0[sample(d, d)], names(x0)))

ti <- seq(0, 5, by = .01)
eq_traj <- numsolve(o = rat, time = rep(ti, d), x0 = x0s, 
                    param = list(theta1 = as.vector(K1) * sc1, theta2 = as.vector(K2) * sc2))
eq_traj

# Prepare correct complexes and extra complexes
A_correct <- A
sc1_correct <- sc1
sc2_correct <- sc2
A_extra <- rbind(A, A[-1, sample(ncol(A))])
sc1_extra <- matrix(1, nrow = nrow(A_extra) * nrow(CC), ncol = d)
sc2_extra <- matrix(1, nrow = nrow(A_extra) * nrow(CC), ncol = d)
for (target in seq_len(d)) {
  inh_react <- CC[, target] != 0
  sc1_extra[, target] <- as.vector(1 - outer(inh_react, rep(1, nrow(A_extra))))
  sc2_extra[, target] <- as.vector(1 - outer(inh_react, rep(1, nrow(A_extra))))
}

# Order to choose E (choose whatever species entering most complexes (true and extra) not yet inhibited)
AA <- Matrix(A_extra != 0)
apply(AA, 2, sum)
Es <- numeric()
#Es[1] <- which.max(apply(AA, 2, sum))
for (i in 1:20) {
  Es[i] <- which.max(apply(AA[!apply(AA[, Es, drop = FALSE], 1, any), , drop = FALSE], 2, sum))
}
Es <- rev(Es)


trajs <- split(eq_traj, c(0, cumsum(diff(eq_traj[, 1]) < 0)))
trajs <- lapply(trajs, matrix, ncol = d + 1)
trajs <- lapply(trajs, function(tt) {colnames(tt) <- c("Time", names(x0)); tt})
ggplot2::ggplot(reshape::melt(data.frame(trajs[[2]]), id.vars = "Time"), 
                aes(x = Time, y = value, color = variable)) +
  geom_line() + facet_wrap(~ variable)
save(trajs, sc1_correct, sc2_correct, A_correct, A_extra, sc1_extra, sc2_extra, Es, file = "trajs_inter.RData")


