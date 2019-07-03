##--------------------------------##
##  Small ODE example: EnvZ/OmpR  ##
##--------------------------------##

rm(list = ls())

library(episode)
library(ggplot2); library(reshape)
library(igraph)
source("../ggplot_theme.R")

## Define system ##
x0 <- c('(EnvZ-P)OmpR' = 4, 
        'EnvZ(OmpR-P)' = 5, 
        'EnvZ-P' = 6,
        'EnvZ' = 7,
        'OmpR' = 8,
        'OmpR-P' = 9)
set.seed(27 + 09)
x0 <- c('(EnvZ-P)OmpR' = 0, 
        'EnvZ(OmpR-P)' = 0, 
        'EnvZ-P' = 0,
        'EnvZ' = 0,
        'OmpR' = 0,
        'OmpR-P' = 0)
x0[] <- runif(6, min = 5, max = 10)
A <- matrix(c(0, 0, 1, 0, 1, 0,
              1, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 1,
              0, 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0,
              0, 0, 1, 0, 0, 0), 
            byrow = TRUE, ncol = 6)
B <- matrix(c(1, 0, 0, 0, 0, 0,
              0, 0, 1, 0, 1, 0,
              0, 0, 0, 1, 0, 1,
              0, 1, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 1,
              0, 0, 0, 1, 1, 0,
              0, 0, 1, 0, 0, 0,
              0, 0, 0, 1, 0, 0), 
            byrow = TRUE, ncol = 6)
k <- c('k1' = 0, 
       'k-1' = 0,
       'kt' = 0,
       'k2' = 0,
       'k-2' = 0,
       'kp' = 0,
       'kk' = 0,
       'k-k' = 0)
k[] <- rnorm(8, 3)

colnames(B) <- colnames(A) <- names(x0)
rownames(B) <- rownames(A) <- names(k)

m <- mak(A = A, B = B)
ti <- seq(0, 1, by = 0.01)
sc <- cbind(1, 1, c(rep(0, 3), rep(1, 5)), c(rep(1, 3), rep(0, 3), rep(1, 2)))
trajs <- numsolve(m, time = rep(ti, 4), 
                  x0 = cbind(x0, 1.5 * sample(x0), x0, x0),
                  param = k * sc)
trajs <- lapply(split(trajs, c(0, cumsum(diff(trajs[, 1]) < 0))), 
                matrix, ncol = ncol(trajs))
trajs <- lapply(trajs, function(tt) {colnames(tt) <- c("Time", names(x0)); tt})
# original system, perturbed system, intervened 1, intervened 2


# ## trajectory
# ggplot(melt(data.frame(trajs[[2]]), id.vars = "Time"), aes(x = Time, y = value, color = variable)) +
#   geom_line()
library(Matrix)
## true network
netw <- field(m, x = x0, param = k, differentials = TRUE)$f_dx != 0
diag(netw) <- FALSE
netw
m
gg <- graph_from_adjacency_matrix(t(netw))
plot(gg, layout = layout_(gg, in_circle()), vertex.shape = rep("circle", 6),
     vertex.color = "skyblue", edge.arrow.size=.5)



## larger frame model: x -> y and x + y -> z and z -> x + y
d <- ncol(A)
A <- B <- matrix(0, ncol = d, nrow = 0)
for (x in seq_len(d)) {
  for (y in setdiff(seq_len(d), x)) {
    A_ <- matrix(0, nrow = 1 + d - 2, ncol = d)
    B_ <- matrix(0, nrow = 1 + d - 2, ncol = d)
    A_[1, x] <- 1
    B_[1, y] <- 1
    s <- 1
    for (z in setdiff(seq_len(d), c(x, y))) {
      s <<- s + 1
      if (x > y) {
        A_[s, c(x, y)] <- 1
        B_[s, z] <- 1
      } else {
        B_[s, c(x, y)] <- 1
        A_[s, z] <- 1
      }
    }
    A <- rbind(A, A_)
    B <- rbind(B, B_)
  }
}



## Sample data and plot data
data <- trajs
data <- lapply(data, function(dat) {
  dat[, -1] <- dat[, -1] + matrix(rnorm(length(dat[, -1]), sd = 1), nrow = nrow(dat))
  dat[0:30 * 3 + 1, ]
})

## plot data
ii <- 0
data_ <- lapply(data, function(dat) {
  ii <<- ii + 1
  cbind(dat, Type = ii)
})
dat <- data.frame(do.call(rbind, data_))
names(dat) <- c("Time", names(x0), "Type")
dat$Type <- factor(dat$Type)
levels(dat$Type) <- c("Original System", "Perturbed System", "(EnvZ-P)OmpR Inhibited", "EnvZ(OmpR-P) Inhibited")
gp <- ggplot(melt(dat, id.vars = c("Time", "Type")), 
             aes(x = Time, y = value, color = variable)) +
  geom_point() + geom_line() + ylab("Abundance") +
  scale_color_discrete("") + facet_wrap(~Type) + xlab("Time [min]") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 17, face = "bold"),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5))
pdf(paste0("../figures/EnvZ_data.pdf"), height = 6, width = 9)
print(gp)
dev.off()


## Recovery
sc_full <- cbind(1, 1, 
                 !apply(cbind(A[, c(1)], (B - A)[, c(1)]) != 0, 1, any),
                 !apply(cbind(A[, c(2)], (B - A)[, c(2)]) != 0, 1, any))


## Recovery via AIM ##
m_full <- mak(A, B)

m_full

x_smooth <- lapply(data, function(dat) {
  cbind(Time = dat[, 1], apply(dat[, -1], 2, function(y) {
    ll <- loess(y ~ x, data.frame(y = y, x = dat[, 1]), 
                family = "symmetric")
    predict(ll, newdata = data.frame(x = dat[, 1]))
  }))
})

## Perturbed
a_per <- aim(m_full, 
             op = opt(do.call(rbind, data[c(1, 2)])),
             x = do.call(rbind, x_smooth[c(1, 2)]))#, adapts = NULL)

## Intervened
a_int <- aim(mak(A, B, 
                 r = reg(contexts = sc_full[, c(3, 4)])), 
             op = opt(do.call(rbind, data[c(3, 4)])), 
             x = do.call(rbind, x_smooth[c(3, 4)]))#, adapts = NULL)

## Both
a_bot <- aim(mak(A, B, 
                 r = reg(contexts = sc_full)), 
             op = opt(do.call(rbind, data)), 
             x = do.call(rbind, x_smooth))

aas <- list(per = a_per, int = a_int, bot = a_bot)
rods <- lapply(aas, rodeo)

# No refitting, since only one smoother
netw_guess <- lapply(aas, function(a) {
  sel <- which.max(Matrix::colSums(a$params$rate!=0) >= 8)
  netw_guess <- field(m_full, x = x0, param = a$params[[1]][, sel],
                      differentials = TRUE)$f_dx != 0
  diag(netw_guess) <- FALSE
  netw_guess
})





## igraph it
gg_s <- lapply(netw_guess, function(ntw) {
  gg_guess <- graph_from_adjacency_matrix(t(ntw))  
  gg_u <- union(gg, gg_guess)
  E(gg_guess)
  E(gg_u)$weight <- as.numeric((attr(E(gg_u), "vnames") %in% attr(E(gg_guess), "vnames")) & 
                                 (attr(E(gg_u), "vnames") %in% attr(E(gg), "vnames"))) - 
    as.numeric((attr(E(gg_u), "vnames") %in% attr(E(gg_guess), "vnames")) &
                 !(attr(E(gg_u), "vnames") %in% attr(E(gg), "vnames")))
  E(gg_u)$color <- rep('gray', length(E(gg_u)))
  E(gg_u)$color[E(gg_u)$weight == 1] <- 'green'
  E(gg_u)$color[E(gg_u)$weight == 0] <- 'gray'
  E(gg_u)$color[E(gg_u)$weight == -1] <- 'red'
  gg_u
})




lay <- layout_on_grid(gg)
igraph_options(vertex.size = 65, vertex.label.cex = .95, edge.arrow.size = 1)



lay[, 2] <- c(.25, 0, .25, .75, 1, .75)
pdf(paste0("../figures/EnvZ_est_per.pdf"), height = 7, width = 10)
plot(gg_s$per, layout = lay,
     vertex.shape = rep("crectangle", 6),
     vertex.color = "white")
dev.off()

pdf(paste0("../figures/EnvZ_est_int.pdf"), height = 7, width = 10)
plot(gg_s$int, layout = lay,
     vertex.shape = rep("crectangle", 6),
     vertex.color = "white")
dev.off()

pdf(paste0("../figures/EnvZ_est_both.pdf"), height = 7, width = 10)
plot(gg_s$bot, layout = lay,
     vertex.shape = rep("crectangle", 6),
     vertex.color = "white")
dev.off()

pdf(paste0("../figures/EnvZ_truth.pdf"), height = 7, width = 10)
plot(gg, layout = lay,
     vertex.shape = rep("crectangle", 6),
     vertex.color = "white")
dev.off()



# figure of ranking #
xs <- sapply(seq(.25, 1.25, length.out = 5), function(sp) {
  lapply(data, function(dat) {
    cbind(Time = dat[, 1], apply(dat[, -1], 2, function(y) {
      ll <- loess(y ~ x, data.frame(y = y, x = dat[, 1]), span = sp, 
                  family = "symmetric")
      predict(ll, newdata = data.frame(x = dat[, 1]))
    }))
  }) 
}, simplify = FALSE)

as_both <- lapply(xs, function(x_smooth) {
  a <- aim(mak(A, B, 
               r = reg(contexts = sc_full)), 
           op = opt(do.call(rbind, data)), 
           x = do.call(rbind, x_smooth))
  rod <- rodeo(a)
  rod
})

ntws <- lapply(as_both, function(rod) {
  apply(rod$params$rate, 2, function(k) {
    netw_guess <- field(m_full, x = x0, param = k,
                        differentials = TRUE)$f_dx != 0
    diag(netw_guess) <- FALSE
    netw_guess
  })
})

## all lambda sequences are (up to a scale equivalent)
common_lambda <- as_both[[3]]$op$lambda


ss <- 0
rank <- lapply(as_both[-1], function(rod) {
  ss <<- ss + 1
  cbind(loss = rod$loss, df = apply(rod$params$rate != 0, 2, sum), coor = seq_along(rod$loss), smooth = ss, lambda = common_lambda)[1:18, ]
})

rank <- do.call(rbind, rank)

choice <- sapply(sort(unique(rank[, "df"])), function(df) {
  red <- rank[rank[, "df"] == df, , drop = FALSE]
  red[which.min(red[, "loss"]), c("coor", "smooth")]
})

dfr <- data.frame(rank)

dfr$opt <- apply(rank[, c("coor", "smooth")], 1, function(x) {
  any(apply(choice == x, 2, all))
})

head(dfr)
dfr$label <- rep("", nrow(dfr))
dfr$label[dfr$opt] <- dfr$df[dfr$opt]

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
dfr$smooth <- factor(dfr$smooth)

dfr1 <- dfr
dfr1$smooth <- factor(dfr1$smooth, levels = rev(levels(dfr1$smooth)))
g1 <- ggplot(dfr1, aes(x = log(lambda), y = smooth, fill = smooth, alpha = df, color = opt)) +
  geom_tile(width = 0.365, height = 0.9, size = 1) + 
  ggplot_theme + ylab("Smoother") + xlab(expression(log(lambda))) +
  geom_text(aes(label = label), alpha = 1) +
  scale_fill_manual(values = cbPalette, guide = FALSE) +
  theme(legend.position = "right") +
  scale_alpha_continuous(name = "Number of\nparameters") +
  scale_color_manual(values = c("white", "black"), guide = FALSE)
g1

dfr1$smooth <- factor(dfr1$smooth, levels = rev(levels(dfr1$smooth)))
g2 <- ggplot(dfr1, aes(x = df, y = loss, color = smooth)) + 
  geom_point(aes(shape = opt), size = 2, stroke = 2) + ggplot_theme +
  theme(legend.position = "right") +
  xlab("Number of parameters") + ylab("Loss") +
  scale_shape_manual(values = c(4, 1), guide = FALSE) +
  scale_color_manual(values = rev(cbPalette), name = "Smoother   ")

g2
library(gridExtra)


gg <- grid.arrange(g1, g2)

pdf(paste0("../figures/Alg_visual.pdf"), height = 7, width = 12)
grid.arrange(g1, g2)
dev.off()
