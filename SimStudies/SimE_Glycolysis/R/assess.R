
# linear interpolation
lip <- function(x, y, xout, pad = TRUE) {
  ox <- order(x)
  x <- x[ox]
  y <- y[ox]
  if (pad) {
    padl <- c(0, x[1]); padr <- c(x[length(x)], 1)
    x <- c(padl, x, padr) 
    y <- c(padl, y, padr)
  }
  approx(x = x, y = y, xout = xout)$y
}

# fpr and tpr are equal length vectors, 
auroc <- function(fpr, tpr) {
  xs <- seq(0, 1, by = 0.05)
  ys <- lip(x = fpr, y = tpr, xout = xs, pad = FALSE)
  sum((ys[-length(ys)] + diff(ys) / 2) * diff(xs))
}
