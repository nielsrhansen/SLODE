library(tsars)
library(aim)
library(magrittr)
library(ggplot2)
library(reshape)

context("gpfit")

set.seed(30032017)


test_that("Trigonometric test", {
  ## Truth and obs
  tys   <- seq(-1, 1, by = .1) %>%
    { c(., cos(.), exp(.)) } %>%
    matrix(ncol = 3) %>%
    set_colnames(c("Time", "cos", "exp"))
  tdys  <- tys[, 1] %>%
    { c(., -sin(.), exp(.)) } %>%
    matrix(ncol = 3) %>%
    set_colnames(c("Time", "cos", "exp"))
  ys  <- tys

  ## Add noise and get gps
  ys[, -1] %<>% '+'(length(.) %>% rnorm(sd = .5))
  gp_res <- ys %>% gpfit(tout = seq(-1.5, 1.5, by = .01))

  ## Check names are transfered
  expect_identical(
    lapply(gp_res[-3], colnames),
    setNames(replicate(2, c("Time", colnames(ys)[-1]), simplify = FALSE), c("ms", "dms"))
    )

  ## Prepare for plot
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
  expect_equal(mean(gpdf$value), 0.75806063353071007072)

  trip <- ggplot(tydf, aes(x = Time, y = value)) +
    geom_line(linetype = "dashed") +
    geom_point(data = ydf, aes(x = Time, y = value)) +
    geom_line(data = gpdf, aes(x = Time, y = value)) +
    facet_grid(fun ~ variable) + theme_bw() + ylab("")
  # print(trip)
  # ggsave("tests/testthat/plots/trip_ref.pdf", plot = trip, width = 7, height = 7)
  # ggsave("tests/testthat/plots/trip.pdf", plot = trip, width = 7, height = 7)


})


test_that("MAK test", {
  A <- matrix(c(1, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 1, 0,
    0, 1, 0, 0,
    0, 0, 0, 1), ncol = 4, byrow = TRUE)
  B <- matrix(c(0, 0, 1, 0,
    1, 1, 0, 0,
    1, 0, 0, 1,
    0, 0, 0, 1,
    0, 0, 1, 0), ncol = 4, byrow = TRUE)
  m <- mak(A, B)
  k <- c(1, 2, 0.5, 0, 0); x0 <- c(E = 2, S = 8, ES = 0.5, P = 0.5)
  Time <- seq(0, 10, by = .1)

  ## Simulate data
  tys <- numsolve(m, Time, x0, k)
  ys  <- tys %>% extract(c(1, 2, 4, 8, 16, 32, 64, 101), , drop = FALSE)
  ys[, -1] %<>% '+'(length(.) %>% rnorm(sd = .5))

  ## Get GP fit
  gp_res <- gpfit(t(t(ys) - c(0, colMeans(ys[, -1]))), Time)
  gp_res$ms <- t(t(gp_res$ms) + c(0, colMeans(ys[, -1])))


  ## Evaluate true differential
  tdys <- cbind(tys %>% extract(, 1),
                tys %>% extract(, -1, drop = FALSE) %>% apply(1, field, param = k, o = m) %>% t) %>%
    set_colnames(c("Time", names(x0)))

  ## Prepare for plot
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
  expect_equal(mean(gpdf$value), 1.23815568863593883364)

  makp <- ggplot(tydf, aes(x = Time, y = value)) +
    geom_line(linetype = "dashed") +
    geom_point(data = ydf, aes(x = Time, y = value)) +
    geom_line(data = gpdf, aes(x = Time, y = value)) +
    facet_grid(fun ~ variable, scales = "free_y") + theme_bw() + ylab("")
  # print(makp)
  # ggsave("tests/testthat/plots/makp_ref.pdf", plot = makp, width = 7, height = 7)
  # ggsave("tests/testthat/plots/makp.pdf", plot = makp, width = 7, height = 7)








  # multiple series and gpe test
  gp_res <- rbind(ys) %>% gpe(tout = rep(seq(-1.5, 15, by = .01), 1), p = c(.0100, 1))

  gp_res$ll

  ## Check names are transfered
  expect_identical(
    lapply(gp_res[-3], colnames),
    setNames(replicate(2, c("Time", colnames(ys)[-1]), simplify = FALSE), c("ms", "dms"))
  )

  ## Prepare for plot
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

  trip <- ggplot(tydf, aes(x = Time, y = value)) +
    geom_line(linetype = "dashed") +
    geom_point(data = ydf, aes(x = Time, y = value)) +
    geom_line(data = gpdf, aes(x = Time, y = value)) +
    facet_grid(fun ~ variable, scales = "free_y") + theme_bw() + ylab("")
  print(trip)

})
