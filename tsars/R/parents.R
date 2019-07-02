##------------------------------------##
##  Generate all signed parents set   ##
##------------------------------------##

#' Generate all signed parent sets of certain size.
#'
#' @param d Positive integer giving total number of species.
#' @param i Positive integer, index of target.
#' @param np Positive integer, number of parents.
#' @return A matrix, each column representing a unique signed parent set.
#' @examples
#' parents(5, 2, 2)
#'
#' library(magrittr)
#' p_max   <- 2
#' target  <- 1
#' m       <- mak_enzyme(5)
#' pas     <- sapply(seq_len(p_max), parents, d = ncol(m$A), i = target)
#' models  <- pas %>%
#'   lapply(. %>% apply(2, find_reactions, i = target, makobj = m)) %>%
#'   do.call(cbind, .)
#' # each column is one of the possible models where 'target' has between 1 and p_max parents
#' # each row represent a row (reaction) in A and B
#'
#' # Non-magrittr version
#' models  <- do.call(cbind, lapply(pas, function(pa) apply(pa, 2, find_reactions, i = target, makobj = m)))
#' @export
parents <- function(d, i, np) {
  if (!(i %in% seq_len(d))) stop(paste0("i not in 1:d"))
  if (np < 1) stop(paste0("np is smaller than 1"))
  if (np > d - 1) stop(paste0("np is larger than d-1"))

  p <- combn(setdiff(seq_len(d), i), np)
  s <- t(do.call(expand.grid, replicate(np, c(1, -1), simplify = FALSE)))
  return(matrix(apply(s, 2, "*", x = p), nrow = np))
}


