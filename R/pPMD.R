#' Projected PMD of a matrix.
#'
#' @param Data the rectangular matrix to decompose ;
#' @param k the desired rank of the singular decomposition (default to 2) ;
#' @param rdsLeft a vector of radiuses
#' (>0) of the $L_1$ or $L_G$ balls for each of the k left vectors ;
#' @param rdsRight a vector of radiuses
#' (>0) of the $L_1$ or $L_G$ balls for each of the k right vectors ;
#' @param tol.si tolerance for the sparsity
#'
#' @return Pseudo-singular vectors and values
#' @examples
#' X <- matrix(rnorm(20), 5, 4)
#' pPMD(X)
#' pPMD(
#'   X,
#'   k = 3L,
#'   rdsLeft = rep(0.5 * sqrt(5), 3),
#'   rdsRight = rep(0.5 * sqrt(4), 3))
#' @author Vincent Guillemot, Ju-Chi Yu
#' @export

pPMD <- function(Data, k = 2L,
                 rdsLeft = rep(1, k), rdsRight = rep(1, k),
                 tol.si = .Machine$double.eps) {

  # # Test that the arguments are valid
  # garb <- runTestsPMD(Data, k, init, initLeft, initRight, seed,
  #                  rdsLeft, rdsRight,
  #                  grpLeft, grpRight,
  #                  orthogonality, OrthSpaceLeft, OrthSpaceRight,
  #                  projPriority,
  #                  projPriorityLeft,
  #                  projPriorityRight)

  I <- nrow(Data)
  J <- ncol(Data)

  U <- matrix(NA, I, k)
  V <- matrix(NA, J, k)

  d <- rep(NA, k)

  for (r in 1:k) {
    ## PMD
    suppressWarnings(res.pmd <- PMD(
      x = Data, # X = Data,
      sumabsu = rdsLeft[r], # rdsLeft[r]
      sumabsv = rdsRight[r], # rdsRight = rdsRight[r]
      K = 1, center = FALSE,
      type = "standard",
      trace = FALSE))

    if (all(res.pmd$u == 0) | all(res.pmd$v == 0)) {
      stop ("Too many components are estimated. Try extracting fewer components.")
    }
    ## Get weights and pseudo singular values
    U[, r] <- res.pmd$u
    V[, r] <- res.pmd$v
    d[r] <- res.pmd$d

    Ur <- U[ , r, drop = FALSE]
    Vr <- V[ , r, drop = FALSE]

    ## Post-hoc orthogonalization
    Data <- (diag(I) - Ur %*% t(Ur)) %*% Data %*% (diag(J) - Vr %*% t(Vr))
  }

  oD <- order(d, decreasing = TRUE)
  # oD <- 1:R
  res <- list(d = d[oD], u = U[, oD], v = V[, oD])

  res$rdsLeft <- rdsLeft
  res$rdsRight <- rdsRight

  return(res)
}

#' @export
#' @keywords internal
runTestsPMD <- function(X, k, init, initLeft, initRight,
                         rdsLeft, rdsRight,
                         grpLeft, grpRight,
                         projPriority,
                         projPriorityLeft,
                         projPriorityRight,
                         itermaxALS, itermaxPOCS,
                         epsALS, epsPOCS) {

  ##### Test X ####
  if (nrow(X)==1 & ncol(X)==1)
    stop("You are attempting a gsGSVD of a scalar.")

  if (any(is.na(X)))
    stop("X should not contain missing values")

  ##### Test R ####
  if (!is.integer(k)) stop("R should be an integer.")
  if (k <= 1) stop("K should be > 1.")

  ##### Test initialization ####
  if (is.null(init)) {
    if (is.null(initLeft) | ! is.matrix(initLeft))
      stop("initLeft should be a matrix.")
    if (is.null(initRight)  | ! is.matrix(initRight))
      stop("initRight should be a matrix.")
  }
  if (! init %in% c("svd", "rand"))
    stop("init should be either svd or rand.")

  return(NULL)
}

#' @export
#' @keywords internal
initializePMD <- function(X, I, J, k, init, initLeft, initRight, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (any(c(init, initLeft, initRight) == "svd")) {
    svdx <- svd(X, nu=k, nv=k)
  }

  if (is.null(init)) {
    if (initLeft == "svd") {
      U0 <- svdx$u
    } else if (initLeft == "rand") {
      U0 <- 1/(I-1) * MASS::mvrnorm(n = I, mu = rep(0,k),
                              Sigma = diag(k), empirical = TRUE)
    } else {
      U0 <- initLeft
    }

    if (initRight == "svd") {
      V0 <- svdx$u
    } else if (initRight == "rand") {
      V0 <- 1/(I-1) * MASS::mvrnorm(n = I, mu = rep(0,k),
                              Sigma = diag(k), empirical = TRUE)
    } else {
      V0 <- initRight
    }
  } else if (init == "svd") {
    U0 <- svdx$u
    V0 <- svdx$v
  } else if ( init=="rand") {
    U0 <- 1/(I-1) * MASS::mvrnorm(n = I, mu = rep(0,k),
                            Sigma = diag(k), empirical = TRUE)
    V0 <- 1/(J-1) * MASS::mvrnorm(n = J, mu = rep(0,k),
                            Sigma = diag(k), empirical = TRUE)
  } else {
    stop("Unkown error, contact support!")
  }

  return(list(U0 = U0, V0 = V0))
}

#' @export
#' @keywords internal
sparseIndex <- function(res.ppmd, singularValues, tol = 1e-10) {
  R <- length(res.ppmd$d)
  singularValues <- singularValues[1:R]
  U <- res.ppmd$u
  V <- res.ppmd$v
  U.sq <- U^2
  V.sq <- V^2
  ctrLeft <- U.sq
  ctrRight <- V.sq

  I <- NROW(ctrLeft)
  J <- NROW(ctrRight)
  rdsLeft <- res.ppmd$rdsLeft
  rdsRight <- res.ppmd$rdsRight

  r1 <- compute.fit(singularValues, res.ppmd$d)

  # Compute the sparsity part of the index
  n0inU <- cumsum(colSums(ctrLeft <= tol))
  n0inV <- cumsum(colSums(ctrRight <= tol))
  radiusIndexLeftG <- cumgmean(rdsLeft / sqrt(I))
  radiusIndexRightG <- cumgmean(rdsRight / sqrt(J))
  radiusIndexLeftA <- cummean(rdsLeft / sqrt(I))
  radiusIndexRightA <- cummean(rdsRight / sqrt(J))

  r2 <- n0inU / (I * (1:R))
  r3 <- n0inV / (J * (1:R))
  r4 <- (n0inU + n0inV) / ((I + J) * (1:R))
  # Combine
  SI <- r1 * r4
  SIleft <- r1 * r2
  SIright <- r1 * r3

  return(list(
    SI = SI, SIleft = SIleft, SIright = SIright,
    r1 = r1, r2 = r2, r3 = r3, r4 = r4,
    n0inU = n0inU, n0inV = n0inV,
    rdsLeft = rdsLeft, rdsRight = rdsRight,
    radiusIndexLeftG = radiusIndexLeftG,
    radiusIndexRightG = radiusIndexRightG,
    radiusIndexLeftA = radiusIndexLeftA,
    radiusIndexRightA = radiusIndexRightA))
}

#' @export
#' @keywords internal
compute.fit <- function(d, pseudo.d) {
  r1 <- cumsum(pseudo.d ^ 2) / cumsum(d ^ 2)
  return(r1)
}

#' @export
#' @keywords internal
gmean <- function(x, na.rm = TRUE) {
  if (any(abs(x) < 2*.Machine$double.eps)) return(0)
  if (any(x < 0)) return(0)
  exp(mean(log(x), na.rm = na.rm))
}

#' @export
#' @keywords internal
cummean <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  return(cumsum(x) / seq_along(x))
}

#' @export
#' @keywords internal
cumgmean <- function(x, na.rm = TRUE) {
  if (any(abs(x) < 2*.Machine$double.eps)) return(0)
  if (any(x < 0)) return(0)
  exp(cummean(log(x), na.rm = na.rm))
}

