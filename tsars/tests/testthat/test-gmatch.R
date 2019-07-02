library(tsars)
library(aim)
library(magrittr)
library(reshape)
library(ggplot2)

context("gmatch")

set.seed(03042017)

test_that("X full", {
  m <- mak_enzyme(3, empty = FALSE)
  shuff <- sample(seq_len(nrow(m$A)))
  m$A <- m$A[shuff, ]; m$B <- m$B[shuff, ];
  k <- c(5:1, rep(0, nrow(m$A) - 5))
  x0 <- c(E = 8, S = 2, ES = 0.5)#, H = 5, J = 1.4)
  Time <- seq(0, .1, by = .001)

  ## Simulate data
  tys <- numsolve(m, Time, param = list(k = k), x0)
  ys  <- tys %>% extract(Time %>% length %>% '-'(1) %>% '/'(10) %>% floor %>% seq(0, .) %>% '*'(10) %>% '+'(1), , drop = FALSE) # every tenth obs
  ys[, -1] %<>% '+'(length(.) %>% rnorm(sd = 1))

  ## Get GP fit
  gp_res <- gpfit(ys, ys[, 1])

  ## Plot
  tdys <- cbind(tys %>% extract(, 1),
    tys %>% extract(, -1, drop = FALSE) %>% apply(1, field, param = k, o = m) %>% t) %>%
    set_colnames(c("Time", names(x0)))

  colnames(tys) <- colnames(tdys)
  tydf <- rbind(tys %>% data.frame %>% cbind(fun = "id"),
    tdys %>% data.frame %>% cbind(fun = "diff"))
  tydf$fun %<>% factor
  tydf %<>% melt(id.vars = c("Time", "fun"))

  ydf  <- ys %>% data.frame %>% cbind(fun = "id")
  ydf$fun %<>% factor
  ydf %<>% melt(id.vars = c("Time", "fun"))

  gpdf <- gp_res %$%
    rbind(ms %>% data.frame %>% cbind(fun = "id"),
      dms %>% data.frame %>% cbind(fun = "diff"))
  gpdf %<>% melt(id.vars = c("Time", "fun"))

  makp <- ggplot(tydf, aes(x = Time, y = value)) +
    geom_line(linetype = "dashed") +
    geom_point(data = ydf, aes(x = Time, y = value)) +
    geom_line(data = gpdf, aes(x = Time, y = value)) +
    facet_grid(fun ~ variable, scales = "free_y") + theme_bw() + ylab("")
  #print(makp)

  ## Design test
  X <- X_full(m, gp_res$ms)
  expect_identical(dim(X), c(nrow(ys) * ncol(m$A), nrow(m$A)))
  expect_identical(X %*% k %>% as.vector,
    gp_res$ms %>% extract(, -1, drop = FALSE) %>%
      apply(1, field, param = k, o = m) %>%
      as.vector)

})


