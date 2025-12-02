library(roxygen2)
library(pedigreeTools)
library(Matrix)

#' Compute genomic relationship matrix G
#'
#' @param M Genotype matrix (individuals x markers).
#' @param eps Small numeric value added to the diagonal.
#'
#' @return Genomic relationship matrix G.
#' @export
getG <- function(M, eps = eps) {
  # Alles frequencies
  p <- apply(M, 2, mean)/2
  q <- 1 - p
  
  # Creating P matrix
  #P <- matrix(rep(p*2, nrow(M)), ncol = ncol(M), nrow = nrow(M), byrow = TRUE)
  
  #Z <- M - P
  #Z <- as.matrix(M - P)
  Z <- sweep(x = as.matrix(M), MARGIN = 2, STATS = 2*p, FUN = "-")
  
  # scaled factor k = 2*sum*(p*q)
  k <- 2*sum(p*q)
  
  # G matrix
  G <- (Z%*%t(Z))/k 
  diag(G) <- diag(G) + eps
  
  return(G)
} 


#' Build H matrix directly
#'
#' @param A Pedigree relationship matrix.
#' @param G Genomic relationship matrix.
#' @param geno Vector of genotyped IDs.
#'
#' @return Hybrid relationship matrix H (direct method).
#' @export
directh <- function(A, G, geno.idx) {
  
  geno.idx <- colnames(A)[colnames(A) %in% geno.idx]
  
  all_ids  <- rownames(A)
  if (!all(geno.idx %in% all_ids))
    stop("Some geno IDs not found in A.")
  
  non.geno.idx <- setdiff(all_ids, geno.idx)
  
  A11 <- A[non.geno.idx, non.geno.idx, drop = FALSE]
  A22 <- A[geno.idx,     geno.idx,     drop = FALSE]
  A12 <- A[non.geno.idx, geno.idx,     drop = FALSE]
  A21 <- t(A12)
  
  ## You need the *correct* algebra here; assuming your formula is what you want:
  B <- solve(A22, A21)         # A22^{-1} A21
  D <- solve(A22, G - A22)     # A22^{-1} (G - A22)
  
  H11 <- A11 + A12 %*% D %*% B
  H12 <- A12 %*% solve(A22, G)
  H21 <- t(H12)
  H22 <- G
  
  Htop    <- cbind(H11, H12)      # rows: non_geno,  cols: non_geno, geno
  Hbottom <- cbind(H21, H22)      # rows: geno,      cols: non_geno, geno
  H       <- rbind(Htop, Hbottom) # rows: c(non_geno, geno)
  
  # Now reorder H back to the original order of A:
  H.direct <- H[rownames(A),rownames(A)]
  
  return(H.direct)
  
}


#' Apply A^{-1} using Cholesky factor
#'
#' @param Li Cholesky factor of A^{-1}.
#' @param X Matrix or vector to multiply.
#'
#' @return A^{-1}X
#' @export
oraculumLi <- function(Li, X) { 
  r <- solve(t(Li), X) 
  l <- solve(Li, r) 
  l
}


#' Randomized range finder for H
#'
#' @param Li Cholesky factor of A^{-1}.
#' @param M Genotype matrix.
#' @param geno Vector of genotyped IDs.
#' @param rank Target rank.
#' @param depth Number of power iterations.
#' @param numVectors Number of random vectors.
#' @param eps Numeric value added to the diagonal.
#' @param cent Logical; center Q
#'
#' @return Low-rank basis Q approximating range(H).
#' @export
randRangeFinder_H <- function(Li, M, geno.idx, rank, depth, numVectors, eps = eps, cent = FALSE) {
  
  geno.idx <- rownames(Li)[rownames(Li) %in% geno.idx]
  non.geno.idx <- setdiff(rownames(Li), geno.idx)
  dim <- nrow(Li)
  
  # 1. Random test matrix Omega
  Omega <- matrix(rnorm(dim * numVectors), nrow = dim, ncol = numVectors)
  
  # 2. First multiply: Y = H * Omega
  Y <- oraculumH(Li, M, non.geno.idx, geno.idx, Omega, eps = eps)
  
  # 3. Different iterations
  if (depth > 0) {
    for (i in seq_len(depth)) {
      qrObject <- qr(Y)
      Q <- qr.Q(qrObject)
      Y <- oraculumH(Li, M, non.geno.idx, geno.idx, Q, eps = eps)  # Y = H Q
    }
  }
  
  # 4. Final orthonormal basis Q
  qrObject <- qr(Y)
  Q <- qr.Q(qrObject)
  
  Q[, 1:rank, drop = FALSE]
}


#' Apply H to a matrix X (indirect Colleau method)
#'
#' @param Li Cholesky factor of A^{-1}.
#' @param M Genotype matrix (only genotyped rows used).
#' @param non.geno IDs not genotyped.
#' @param geno Genotyped IDs.
#' @param X Matrix to multiply.
#' @param eps Diagonal inflation.
#'
#' @return H X (indirect computation)
#' @export
oraculumH <- function(Li, M, non.geno.idx, geno.idx, X, eps = eps) {
  n <- nrow(Li)
  if (is.vector(X)) {
    X <- matrix(X, nrow = n, ncol = 1)
  }
  
  numVectors <- ncol(X)
  
  # -------------------------
  # Step 1: Z = A^{-1} X  (via Li)
  # -------------------------
  Z  <- oraculumLi(Li, X)
  Z1 <- Z[non.geno.idx, , drop = FALSE]
  Z2 <- Z[geno.idx, , drop = FALSE]
  
  # -------------------------
  # Step 2: Y2 = GA_22^{-1} z_2 
  #-------------------------
  
  ## 2.1: R = A^{12} Z_2
  subMatrix <- matrix(0, nrow = nrow(Li), ncol = numVectors,
                      dimnames = list(rownames(Li), colnames(X)))
  subMatrix[geno.idx, ] <- as.matrix(Z2)
  tmp <- Li %*% subMatrix
  R2.indirect <- crossprod(Li, tmp)
  R2.indirect <- R2.indirect[non.geno.idx, , drop = FALSE]
  
  ## 2.2: M = A^{22} Z_2
  subMatrix[,] <- 0
  subMatrix[geno.idx, ] <- as.matrix(Z2)
  tmp <- Li %*% subMatrix
  M2.indirect <- crossprod(Li, tmp)
  M2.indirect <- M2.indirect[geno.idx, , drop = FALSE]
  
  ## 2.3: Q = (A^{11})^{-1} R
  Li11   <- Li[, non.geno.idx]
  Ainv11 <- crossprod(Li11)
  Q2.indirect <- solve(Ainv11, R2.indirect)  # CM method
  
  ## 2.4: P = -A^{21} Q
  subMatrix[,] <- 0
  subMatrix[non.geno.idx, ] <- as.matrix(Q2.indirect)
  tmp <- Li %*% subMatrix
  P2.indirect <- crossprod(Li, tmp)
  P2.indirect <- -(P2.indirect[geno.idx, , drop = FALSE])
  
  ## 2.5: f(Z_2) = P + M
  f.Z2 <- P2.indirect + M2.indirect
  
  ## Build G and apply (G + eps I) via Colleau
  if (!identical(rownames(M), rownames(f.Z2))) {
    stop("M and f.Z2 do not have the same row order. Please correct before continuing.")
  }
  
  M.geno <- M[geno.idx, , drop = FALSE]
  p <- apply(M.geno, 2, mean)/2
  q <- 1 - p
  Zg <- sweep(x = as.matrix(M.geno),
              MARGIN = 2, STATS = 2*p, FUN = "-")
  k <- 2 * sum(p * q)
  
  C  <- crossprod(Zg, f.Z2)
  Y2_raw <- Zg %*% C
  Y2 <- (Y2_raw / k) + eps * f.Z2
  
  Y2.indirect <- Y2 
  
  # -------------------------
  # Step 3: D2 = Y2 - Z2
  # -------------------------
  D2 <- Y2.indirect - Z2
  
  # -------------------------
  # Step 4: build D1 
  # -------------------------
  
  ## 4.1: R = A^{12} D_2
  subMatrix[,] <- 0
  subMatrix[geno.idx, ] <- as.matrix(D2)
  tmp <- Li %*% subMatrix
  R1.indirect <- crossprod(Li, tmp)
  R1.indirect <- R1.indirect[non.geno.idx, , drop = FALSE]
  
  ## 4.2: M = A^{22} D_2
  subMatrix[,] <- 0
  subMatrix[geno.idx, ] <- as.matrix(D2)
  tmp <- Li %*% subMatrix
  M1.indirect <- crossprod(Li, tmp)
  M1.indirect <- M1.indirect[geno.idx, , drop = FALSE]
  
  ## 4.3: Q = (A^{11})^{-1} R
  Li11   <- Li[, non.geno.idx]  
  Ainv11 <- crossprod(Li11)
  Q1.indirect <- solve(Ainv11, R1.indirect)
  
  ## 4.4: P = -A^{21} Q
  subMatrix[,] <- 0
  subMatrix[non.geno.idx, ] <- as.matrix(Q1.indirect)
  tmp <- Li %*% subMatrix
  P1.indirect <- crossprod(Li, tmp)
  P1.indirect <- -(P1.indirect[geno.idx, , drop = FALSE])
  
  ## 4.5: f(D_1) = P + M
  f.D1 <- P1.indirect + M1.indirect
  
  ## 4.6: D1 = A^{-1} A_{12} f(D_2)
  subMatrix[,] <- 0
  subMatrix[geno.idx, ] <- as.matrix(f.D1)
  D1.indirect  <- oraculumLi(Li, subMatrix)
  D1.indirect  <- D1.indirect[non.geno.idx, , drop = FALSE]
  
  # -------------------------
  # Step 5: HX = [Y1; Y2]
  # -------------------------
  Y1.indirect <- Z1 + D1.indirect
  
  H.indirect <- matrix(data = 0, 
                       nrow = nrow(Li), 
                       ncol = numVectors,
                       dimnames = list(rownames(Li), 
                                       colnames(X)))
  H.indirect[non.geno.idx, ] <- as.matrix(Y1.indirect)
  H.indirect[geno.idx, ]     <- as.matrix(Y2.indirect)
  
  return(H.indirect)
}


#' PCA on H using randomized SVD
#'
#' @param Li Cholesky factor of A^{-1}.
#' @param M Genotype matrix.
#' @param geno Vector of genotyped IDs.
#' @param rank Target PCA rank.
#' @param depth Number of power iterations.
#' @param numVectors Number of random vectors.
#' @param eps Diagonal inflation.
#' @param cent Logical; center Q?
#'
#' @return List with u, d, v (truncated SVD of H).
#' @export
rhpca <- function(Li, M, geno.idx, rank, depth, numVectors, eps = eps, cent = FALSE) {
  
  geno.idx <- rownames(Li)[rownames(Li) %in% geno.idx]
  non.geno.idx <- setdiff(rownames(Li), geno.idx)
  dim <- nrow(Li)
  
  # 1. Get Q that approximates the range of H
  Q <- randRangeFinder_H(Li, M, geno.idx,
                         rank, depth, numVectors,
                         eps = eps, cent = cent)
  
  # 2. Optionally center Q if requested
  if (cent) {
    Q <- scale(Q, center = TRUE, scale = FALSE)
  }
  
  # 3. Compute HX = H Q
  HX <- oraculumH(Li, M, non.geno.idx, geno.idx, Q, eps = eps)
  
  # 4. Small SVD on B = Q' H Q ≈ Q' HX
  B <- crossprod(Q, HX)   # B is (rank x rank) approx to H in Q basis
  svdObject <- svd(B)
  
  U <- Q %*% svdObject$u
  D <- svdObject$d
  V <- svdObject$v
  
  list(
    u = U[, 1:rank, drop = FALSE],
    d = D[1:rank],
    v = V[, 1:rank, drop = FALSE]
  )
}






