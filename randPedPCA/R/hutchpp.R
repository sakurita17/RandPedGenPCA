




makeRandMat  = function(m, n) 2*matrix(sample(1:2, size = m*n, replace = TRUE, prob = c(0.5, 0.5)), nrow = m, ncol = n) - 3


#' Hutch++ trace estimation
#'
#' @param B An object related to the matrix A for which the trace is to be estimated
#' @param num_queries Number of random vectors to draw
#' @param sketch_frac Hutch++ detail
#' @param center Whether or not to implicitly centre
#' @param oraculum The oracle function to be used
#'
#' @details
#' The Hutch++ algorithm (Meyer et al. 2021, https://doi.org/10.48550/arXiv.2010.09649)
#' estimates the trace of a matrix A by evaluating matrix
#' vector products of A and (sub-gaussian) random vectors. This is used on a
#' matrix B which is related to A through some function. The oracle function has
#' to be chosen so that oracle(B, G) returns the product A %*% G. By default,
#' the oracle function is set to work on a pedigree's L inverse matrix. But this
#' implementation is general and should work - given a custom oracle function -
#' on other input too.
#'
#' In the context of pedigree PCA, this is used to estimate the trace of an
#' (implicitly) centred additive relationship matrix.
#'
#' There logical parameter center allows for a pedigree's L matrix to be
#' (implicitly) centred. This is important because centring changes the total
#' variance of the data and thus the trace of A.
#'
#'
#' @return An estimate of A's trace - numeric
#' @export
#'
#' @examples
#' hutchpp(pedLInv)
#' hutchpp(pedLInv, center=TRUE)
hutchpp <- function(B,
                    num_queries=10,
                    sketch_frac = 2/3,
                    center=FALSE,
                    oraculum=oraculumLi){
  oraculum <- match.fun(oraculum)
  dimension <- dim(B)[1]
  S_num_queries <- round(num_queries * sketch_frac / 2)

  Hutch_num_queries <- num_queries - S_num_queries
  S <- makeRandMat(dimension, S_num_queries)
  Q <- qr.Q(Matrix::qr(oraculum(B, S, center), 0))
  G <- makeRandMat(dimension, Hutch_num_queries)
  G <- G - Q %*% (t(Q) %*% G)

  trace_est <- sum(diag(t(Q) %*% oraculum(B, Q, center))) + sum(diag(t(G) %*% oraculum(B, G, center))) / Hutch_num_queries


  attr(trace_est, "center") <- center
  return(trace_est)
}

