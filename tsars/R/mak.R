##-------------------------------------##
##  Stoichiometric Matrices (A and B)  ##
##-------------------------------------##


#' Stoichiometric Matrices of Simple Enzymatic Reactions
#'
#' @description Generate stoichiometric matrices of mass action kinetics reaction network containing all reactions on the form
#' \deqn{X + Y -> X + Z}
#' where \eqn{Y \neq Z \neq X} and \eqn{X, Y, Z \in 1,...,d} (\eqn{0,...,d} if \code{empty = TRUE})
#' @param d Positive integer giving total number of species.
#' @param empty Logical indicating if empty elements on LHS and RHS are accepted.
#' @param ... Other arguments passed to mak.
#' @return An object of class \code{mak} holding the two stoichiometric matrices A and B.
#' @export
mak_enzyme <- function(d, empty = FALSE, ...) {
  if (d < 2) stop("d must be at least 2.")
  yzs <- combn(ifelse(empty, 0, 1):d, 2)
  yzs <- apply(cbind(yzs, yzs[2:1,]), 2, function(yz) {
    sapply(setdiff(seq(ifelse(empty, 0, 1), d), yz[2]), c, yz)
  })
  yzs <- matrix(yzs, nrow = 3)

  A <- apply(yzs[1:2, ], 2, function(xyz) {
    table(factor(xyz, levels = ifelse(empty, 0, 1):d))
  })
  B <- apply(yzs[2:3, ], 2, function(xyz) {
    table(factor(xyz, levels = ifelse(empty, 0, 1):d))
  })
  if (empty) {
    A <- A[-1, ]
    B <- B[-1, ]
  }
  return(mak(t(A), t(B), ...))
}

#' Join species of stoichiometric matrices
#'
#' @description cbinds the stoichiometric matrices of two \code{mak}-objects. Remaining variables are inherited from the first object.
#' @param mak1 \code{mak}-object.
#' @param mak2 \code{mak}-object.
#' @return An object of class \code{mak} holding the two stoichiometric matrices A and B.
#' @export
join_species <- function(mak1, mak2) {
  if (!is(mak1, "mak")) stop("mak1 is not of class 'mak'.")
  if (!is(mak2, "mak")) stop("mak2 is not of class 'mak'.")
  stopifnot(dim(mak1$A) == dim(mak1$B),
    dim(mak2$A) == dim(mak2$B))
  if (nrow(mak1$A) != nrow(mak2$A)) stop("The two mak-objects do not have same number of rows.")
  mak1$A <- cbind(mak1$A, mak2$A)
  mak1$B <- cbind(mak1$B, mak2$B)
  return(mak1)
}


#' Join reactions of stoichiometric matrices
#'
#' @description rbinds the stoichiometric matrices of two \code{mak}-objects and removes duplets. Remaining variables are inherited from the first object.
#' @param mak1 \code{mak}-object.
#' @param mak2 \code{mak}-object.
#' @return An object of class \code{mak} holding the two stoichiometric matrices A and B.
#' @export
join_reactions <- function(mak1, mak2) {
  if (!is(mak1, "mak")) stop("mak1 is not of class 'mak'.")
  if (!is(mak2, "mak")) stop("mak2 is not of class 'mak'.")
  stopifnot(dim(mak1$A) == dim(mak1$B),
    dim(mak2$A) == dim(mak2$B))
  if (ncol(mak1$A) != ncol(mak2$A)) stop("The two mak-objects do not have same number of columns.")
  AB <- cbind(rbind(mak1$A, mak2$A), rbind(mak1$B, mak2$B))
  AB <- unique(AB)
  mak1$A <- AB[, seq_len(ncol(mak1$A))]
  mak1$B <- AB[, -seq_len(ncol(mak1$A))]
  return(mak1)
}


#' Find Reactions in Stochiometric Matrices
#'
#' @description Find reactions for which a given target species has parents among a certain set.
#' @param makobj Object of class \code{\link{mak}}.
#' @param i Index of target, single valued.
#' @param sj Set of signed parents.
#' @details The target species is indicated via \code{i}. The signed parents \code{sj} must be indeces of the species, but a change of sign is allowed. Positive sign means positive interaction, while negative sign means negative interaction. The target  \code{i} has positive interaction with parent \code{j} in reaction \code{r} iff.  \eqn{B_{ri} - A_{ri} > 0} and \eqn{A_{rj} \neq 0}. Similarly, \code{i} has positive interaction with \code{j} in reaction \code{r} iff. \eqn{B_{ri} - A_{ri} < 0} and \eqn{A_{rj} \neq 0}. No interaction takes place if \eqn{A_{rj}=0}.
#'
#' A reaction is accepted iff it either activates or inhibits \code{i}. A reaction \code{r} activates \code{i} iff. \code{i} interacts only with \code{i} and the positive elements of \code{sj} and at least one positive element of \code{sj} interacts with \code{i}. Similarly a reaction inhibits  \code{i} iff. \code{i} interacts only with \code{i} and the negative elements of \code{sj} and at least one negative element of \code{sj} interacts with \code{i}.
#'
#' @return A logical vector indicating which rows of the stoichiometric matrices represent reactions satisfying target conditions.
#' @export
find_reactions <- function(makobj, i, sj) {
  if (!is(makobj, "mak")) stop("makobj is not of class 'mak'.")
  s <- sign(sj)
  j <- abs(sj)
  if (length(i) != 1) stop("i must have length 1.")
  d <- ncol(makobj$A)
  if (!all(c(i, j) %in% seq_len(d))) stop(paste0("i or some abs(sj) not in 1:", d))

  act <- (makobj$B - makobj$A)[, i] > 0                                         # i must increase
  if (length(setdiff(seq_len(d), c(i, j[s == 1]))) > 0) {
    act <- act & apply(makobj$A[, -c(i, j[s == 1]), drop = FALSE] == 0, 1, all) # no non-i + non-act in lhs of react
  }
  if (length(j[s == 1]) > 0) {
    act <- act & apply(makobj$A[, j[s == 1], drop = FALSE] != 0, 1, any)        # at least one activator must be in lhs of react
  } else {
    act <- act & FALSE                                                          # if no activators none will work
  }

  inh <- (makobj$B - makobj$A)[, i] < 0
  if (length(setdiff(seq_len(d), c(i, j[s == -1]))) > 0) {
    inh <- inh & apply(makobj$A[, -c(i, j[s == -1]), drop = FALSE] == 0, 1, all)# no non-i + non-inh in lhs of react
  }
  if (length(j[s == -1]) > 0) {
    inh <- inh & apply(makobj$A[, j[s == -1], drop = FALSE] != 0, 1, any)       # at least one inhibitor must be in lhs of react
  } else {
    inh <- inh & FALSE                                                          # if no inhibitors none will work
  }

  return(act | inh)
}
