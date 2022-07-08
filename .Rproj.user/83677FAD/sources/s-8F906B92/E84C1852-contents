#' Sparse Correspondence Analysis based on the "projected Penalized Matrix Decomposition"
#'
#' @param DATA the I times J contingency table
#' @param dimensions, integer, the number of dimensions to return (default to 2)
#' @param doublecentering, logical: should the data be double-centered (default to TRUE)
#' @param s1, vector of size 'dimensions" containing the left side regularization parameters, the coefficients in this vector should belong to the interval [1, sqrt(I)], (default to rep(1, dimensions))
#' @param s2, vector of size 'dimensions" containing the left side regularization parameters, the coefficients in this vector should belong to the interval [1, sqrt(J)], (default to rep(1, dimensions))
#'
#' @return An object containing all the necessary outputs for sparse CA
#' @export
#'
#' @examples
#' sCAwithPMD(HairEyeColor[,,1],dimensions = 3L, s1 = rep(0.5 * sqrt(dim(HairEyeColor)[1]), 3), rdsRight = rep(0.5 * sqrt(dim(HairEyeColor)[2]), 3))
#'
sCAwithPMD <- function(
  DATA,
  dimensions = 2L,
  doublecentering = TRUE,
  s1 = rep(1, dimensions),
  s2 = rep(1, dimensions)) {

  # mRP <- ExPosition::makeRowProfiles(DATA, weights = NULL, masses = NULL, hellinger = FALSE)
  N <- sum(DATA)
  X <- 1/N * DATA
  Lv <- rowSums(X)
  Rv <- colSums(X)
  if (doublecentering) X <- X - Lv %*% t(Rv)
  LW <- 1/Lv
  RW <- 1/Rv

  res.gpmd <- GPMD(
    X,
    LW = LW,
    RW = RW,
    k = dimensions,
    rdsLeft = s1,
    rdsRight = s2
  )

  res <- gpmd.out(res.gpmd, X = X, LW = LW, RW = RW)
  res$X.preproc <- X

  return(res)

}

#' @export
#' @keywords internal
gpmd.out <- function(
    res,
    X,
    LW,
    RW) {

  if (!is.vector(LW)) stop("Weights should be vectors!")
  if (!is.vector(RW)) stop("Weights should be vectors!")

  out <- list()
  # SVD
  out$svd$d <- res$d
  out$svd$u <- res$u[, , drop = FALSE]
  out$svd$v <- res$v[, , drop = FALSE]
  # Eigenvalues
  out$eig <- res$d^2
  # GSVD
  out$gsvd$d <- res$d
  out$gsvd$p <- res$p
  out$gsvd$q <- res$q
  out$gsvd$LW <- LW
  out$gsvd$RW <- RW
  # \Dr^{-1}(\Z -\ler\lec\transpose)\Dc^{-1}\Q
  out$fi <- diag(LW) %*% X %*% diag(RW) %*% res$q
  out$fj <- diag(RW) %*% t(X) %*% diag(LW) %*% res$p
  # Dr times F squared divided by eigenvalue
  # out$ci <- diag(1/LW) %*% (out$fi) ^ 2 %*% diag(1/out$eig)
  # out$cj <- diag(1/RW) %*% (out$fj) ^ 2 %*% diag(1/out$eig)
  out$ci <- diag(1/LW) %*% (diag(LW) %*% res$p %*% diag(res$d)) ^ 2 %*% diag(1/out$eig)
  out$cj <- diag(1/RW) %*% (diag(RW) %*% res$q %*% diag(res$d)) ^ 2 %*% diag(1/out$eig)
  # out$ci <- (res$p) ^ 2 -> error because of Hadamard product
  # out$cj <- (res$q) ^ 2 -> error because of Hadamard product
  ## Sparsity
  out$sparsity$rdsLeft <- res$rdsLeft
  out$sparsity$rdsRight <- res$rdsRight
  out$sparsity$SI <- res$SI

  rownames(out$fi) <- rownames(out$ci) <- rownames(out$gsvd$p) <- rownames(out$svd$u) <- rownames(X)
  rownames(out$fj) <- rownames(out$cj) <- rownames(out$gsvd$q) <- rownames(out$svd$v) <- colnames(X)

  colnames(out$fi) <- colnames(out$fj) <- colnames(out$ci) <- colnames(out$cj) <- colnames(out$gsvd$p) <- colnames(out$gsvd$q) <- colnames(out$svd$u) <- colnames(out$svd$v) <- paste0("Component ", seq_along(out$svd$d))
  return(out)
}
