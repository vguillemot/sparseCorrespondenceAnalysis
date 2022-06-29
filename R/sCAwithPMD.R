#' Group sparse Correspondence Analysis
#'
#' @param DATA the contingency table
#' @param k the number of dimensions
#' @param tol
#' @param doublecentering
#' @param init
#' @param initLeft
#' @param initRight
#' @param seed
#' @param rdsLeft
#' @param rdsRight
#' @param orthogonality
#' @param OrthSpaceLeft
#' @param OrthSpaceRight
#' @param projPriority
#' @param projPriorityLeft
#' @param projPriorityRight
#' @param itermaxALS
#' @param itermaxPOCS
#' @param epsALS
#' @param epsPOCS
#'
#' @return
#' @export
#'
#' @examples
#' sCAwithPMD(HairEyeColor[,,1],components = 3L, rdsLeft = rep(0.5 * sqrt(dim(HairEyeColor)[1]), 3), rdsRight = rep(0.5 * sqrt(dim(HairEyeColor)[2]), 3))
#'    
sCAwithPMD <- function(
  DATA, components = 2L,
  doublecentering = TRUE,
  rdsLeft = rep(1, components),
  rdsRight = rep(1, components)) {
  
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
    k = components,
    rdsLeft = rdsLeft,
    rdsRight = rdsRight
  )
  
  res <- gpmd.out(res.gpmd, X = X, LW = LW, RW = RW)
  res$X.preproc <- X

  return(res)

}

gpmd.out <- function(
    res,
    X,
    LW,
    RW) {
  
  if (!is.vector(LW)) stop("Weights should be vectors, damn it!")
  if (!is.vector(RW)) stop("Weights should be vectors, damn it!")
  
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
