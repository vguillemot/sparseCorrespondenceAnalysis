#' @title Sparse Generalized Singular Value Decomposition
#' @description Constrained SVD of a matrix (wrapper of c++ functions).
#'
#' @param X a (data) matrix;
#' @param LW PARAM_DESCRIPTION
#' @param RW PARAM_DESCRIPTION
#' @param k the desired rank of the singular decomposition, Default: 0
#' @param tol PARAM_DESCRIPTION, Default: .Machine$double.eps
#' @param rdsLeft The radius (>0) of the
#' $L_1$ ball for each left vector, Default: rep(1, k)
#' @param rdsRight The radius (>0) of the $L_1$ balls for each right vector, Default: rep(1, k)
#' @param tol.si Tolerance for the computation of the Sparse Index, set by default to .Machine$double.eps
#'
#' @return Pseudo-singular vectors and values
#' @details DETAILS
#' @examples
#' X <- matrix(rnorm(20), 5, 4)
#' GPMD(X, LW = diag(5), RW = diag(4))
#' @rdname sparseGSVD
#' @author Vincent Guillemot, Ju-Chi Yu
#' @export

GPMD <- function(
    X, LW, RW,
    k = 0, # tol = .Machine$double.eps,
    rdsLeft = rep(1, k), rdsRight = rep(1, k),
    tol.si = .Machine$double.eps){

  X_dimensions <- dim(X)

  ## stolen from MASS::ginv()
  if (length(X_dimensions) > 2 || !(is.numeric(X) || is.complex(X))){
    stop("gsvd: 'X' must be a numeric or complex matrix")
  }
  if ( !is.matrix(X) ){
    X <- as.matrix(X)
  }
  
  # a few things about LW for stopping conditions
  LW_is_missing <- missing(LW)

  if ( !LW_is_missing ){

    LW_is_vector <- is.vector(LW)

    if ( !LW_is_vector ){
      if ( nrow(LW) != ncol(LW) | nrow(LW) != X_dimensions[1] ){
        stop("gsvd: nrow(LW) does not equal ncol(LW) or nrow(X)")
      }
      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(LW)){
        stop("gsvd: LW is empty (i.e., all 0s")
      }
    }

    if(LW_is_vector){
      if(length(LW)!=X_dimensions[1]){
        stop("gsvd: length(LW) does not equal nrow(X)")
      }
      # if you gave me all zeros, I'm stopping.
      # if(all(abs(LW)<=tol)){
      if(!are_all_values_positive(LW)){
        stop("gsvd: LW is not strictly positive values")
      }
    }
  }

  # a few things about RW for stopping conditions
  RW_is_missing <- missing(RW)
  if ( !RW_is_missing ){

    RW_is_vector <- is.vector(RW)

    if ( !RW_is_vector ){
      if( nrow(RW) != ncol(RW) | nrow(RW) != X_dimensions[2] ){
        stop("gsvd: nrow(RW) does not equal ncol(RW) or ncol(X)")
      }
      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(RW)){
        stop("gsvd: RW is empty (i.e., all 0s")
      }
    }

    if(RW_is_vector){
      if(length(RW)!=X_dimensions[2]){
        stop("gsvd: length(RW) does not equal ncol(X)")
      }
      # if you gave me all zeros, I'm stopping.
      # if(all(abs(RW)<=tol)){
      if(!are_all_values_positive(RW)){
        stop("gsvd: RW is not strictly positive values")
      }
    }
  }


  ## convenience checks *could* be removed* if problematic
  # convenience checks & conversions; these are meant to minimize LW's memory footprint
  if(!LW_is_missing){
    if( !LW_is_vector){

      if( is_diagonal_matrix(LW) ){
        LW <- diag(LW)
        LW_is_vector <- T  # now it's a vector
      }
    }

    if( LW_is_vector & all(LW==1) ){
      LW_is_missing <- TRUE
      LW <- substitute() # neat! this makes it go missing
    }
  }

  # convenience checks & conversions; these are meant to minimize RW's memory footprint
  if(!RW_is_missing){
    if( !RW_is_vector ){

      if( !RW_is_vector & is_diagonal_matrix(RW) ) {
        RW <- diag(RW)
        RW_is_vector <- TRUE  # now it's a vector
      }
    }

    if( RW_is_vector & all(RW == 1) ) {
      RW_is_missing <- TRUE
      RW <- substitute() # neat! this makes it go missing
    }
  }

  # this manipulates X as needed based on XLW
  if ( !LW_is_missing ) { ## plain SVD
    if ( LW_is_vector ) {
      sqrt_LW <- sqrt(LW)
      X <- X * sqrt_LW
    } else {              ## GSVD
      LW <- as.matrix(LW)
      X <- sqrt_psd_matrix(LW) %*% X
    }
  }


  # this manipulates X (or Y) as needed based on XRW
  if ( !RW_is_missing ) { ## plain SVD
    if( RW_is_vector ) {
      sqrt_RW <- sqrt(RW)
      # X <- sweep(X,2, sqrt_RW,"*") ## replace the sweep with * & t()
      X <- t(t(X) * sqrt_RW)
    } else {              ## GSVD
      RW <- as.matrix(RW)
      X <- X %*% sqrt_psd_matrix(RW)
    }
  }


  if ( !is.integer(k) ){
    k <- as.integer(k) # round down to an integer
  }

  # all the decomposition things
  if (k <= 0) {
    k <- min(X_dimensions)
  }

  res <- pPMD(Data = X, k = k,
              rdsLeft = rdsLeft, 
              rdsRight = rdsRight,
              tol.si = tol.si)

  res$d_full <- res$d
  res$l_full <- res$d_full^2
  # res$tau <- (res$l_full/sum(res$l_full)) * 100
  components.to.return <- min(length(res$d_full), k) # a safety check
  res$d <- res$d_full[1:components.to.return]
  res$l <- res$d^2
  res$u <- res$u[,1:components.to.return, drop = FALSE]
  res$v <- res$v[,1:components.to.return, drop = FALSE]

  # make scores according to weights
  if(!LW_is_missing){
    if(LW_is_vector){

      # res$p <- sweep(res$u,1,1/sqrt_LW,"*")
      res$p <- res$u / sqrt_LW
      # res$fi <- sweep(sweep(res$p,1,LW,"*"),2,res$d,"*")
      res$fi <- t(t(res$p * LW) * res$d)

    } else {

      # res$p <- (LW %^% (-1/2)) %*% res$u
      res$p <- invsqrt_psd_matrix(LW) %*% res$u
      # res$fi <- sweep((LW %*% res$p),2,res$d,"*")
      res$fi <- t(t(LW %*% res$p) * res$d)

    }
  } else {

    res$p <- res$u
    # res$fi <- sweep(res$p,2,res$d,"*")
    res$fi <- t(t(res$p) * res$d)

  }

  if(!RW_is_missing){
    if(RW_is_vector){

      # res$q <- sweep(res$v,1,1/sqrt_RW,"*")
      res$q <- res$v / sqrt_RW
      # res$fj <- sweep(sweep(res$q,1,RW,"*"),2,res$d,"*")
      res$fj <- t(t(res$q * RW) * res$d)

    }else{

      res$q <- invsqrt_psd_matrix(RW) %*% res$v
      # res$fj <- sweep((RW %*% res$q),2,res$d,"*")
      res$fj <- t(t(RW %*% res$q) * res$d)

    }
  } else {

    res$q <- res$v
    # res$fj <- sweep(res$q,2,res$d,"*")
    res$fj <- t(t(res$q)  * res$d)

  }

  rownames(res$fi) <- rownames(res$u) <- rownames(res$p) <- rownames(X)
  rownames(res$fj) <- rownames(res$v) <- rownames(res$v) <- colnames(X)
  # Compute Sparsity Index
  res$rdsLeft <- rdsLeft
  res$rdsRight <- rdsRight
  res.SI <- sparseIndex(res.ppmd = res, singularValues = svd(X, 0, 0)$d, tol = tol.si)
  res$SI <- res.SI
  
  # class(res) <- c("sGSVD", "sSVD", "list")
  return(res)

}



#' @export
#'
#' @title \code{is_diagonal_matrix}: test if a matrix is a diagonal matrix.
#'
#' @description \code{is_diagonal_matrix} takes a matrix and tests if it is a diagonal matrix.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#'
#' @return A boolean. TRUE if the matrix is a diagonal matrix, FALSE if the matrix is not.


## enhanced version of this: https://stackoverflow.com/questions/11057639/identifying-if-only-the-diagonal-is-non-zero
is_diagonal_matrix <- function(x,tol=.Machine$double.eps){
  
  if( length(dim(x)) != 2 ){
    stop("is_diagonal_matrix: x is not a matrix.")
  }
  if( !is.numeric(x) ){
    stop("is_diagonal_matrix: x is not numeric.")
  }
  if( dim(x)[1] != dim(x)[2] ){
    stop("is_diagonal_matrix: x is not a square matrix.")
  }
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }
  
  x[ x^2 < tol ] <- 0
  return(all(x[lower.tri(x)] == 0, x[upper.tri(x)] == 0))
}


#' @export
#'
#' @title \code{is_empty_matrix}: test if a matrix contains all 0s.
#'
#' @description \code{is_empty_matrix} takes a matrix and tests if it contains all 0s.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#'
#' @return A boolean. TRUE if the matrix contains all 0s, FALSE if the matrix does not.


is_empty_matrix <- function(x,tol=.Machine$double.eps){
  
  if( length(dim(x)) != 2 ){
    stop("is_empty_matrix: x is not a matrix.")
  }
  if( !is.numeric(x) ){
    stop("is_empty_matrix: x is not numeric.")
  }
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }
  
  x[abs(x) < tol] <- 0
  
  if(sum(abs(x))==0){
    return(TRUE)
  }else{
    return(FALSE)
  }
  
}

#' @export
#'
#' @title \code{are_all_values_positive}: check for strictly positive values in object
#'
#' @description \code{are_all_values_positive} test presumably a vector, but any numeric object, for strictly positive values
#'
#' @param x A numeric object
#'
#' @return A boolean. TRUE if all values are positive, FALSE otherwise

are_all_values_positive <- function(x){
  
  !(any(is.null(x)) |
      any(is.infinite(x)) |
      any(is.na(x)) |
      any(is.nan(x)) |
      any(x < 0))
  
}

#' @export
#'
#' @title \code{sqrt_psd_matrix}: square root of a positive semi-definite matrix
#'
#' @description \code{sqrt_psd_matrix} takes a square, positive semi-definite matrix and returns the square root of that matrix (via the eigenvalue decomposition by way of \code{tolerance_eigen}).
#'
#' @details Note that \code{sqrt_psd_matrix} uses \code{tolerance_eigen(...,tol=1e-13)} with the tolerance fixed at 1e-13. This enforces an "acceptably zero" value for eigenvalues, and also stops the computation if any eigenvalues are not zero or positive (within the tolerance parameter). If the computation is halted, that means the matrix was not positive semi-definite.
#'
#' @param x A square matrix (presumably positive semi-definite)
#'
#' @return A matrix. The square root of the \code{x} matrix
#' @seealso \code{\link{tolerance_eigen}}

sqrt_psd_matrix <- function(x){
  
  ## checks: just that they are a matrix & square & numeric
  if( length(dim(x)) != 2 ){
    stop("sqrt_psd_matrix: x is not a matrix.")
  }
  if( !is.numeric(x) ){
    stop("sqrt_psd_matrix: x is not numeric.")
  }
  if( dim(x)[1] != dim(x)[2] ){
    stop("sqrt_psd_matrix: x is not a square matrix.")
  }
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }
  
  ## tolerance_eigen
  res <- tolerance_eigen(x, tol = 1e-13)
  
  ## rebuild
  return(t(t(res$vectors) * sqrt(res$values) ) %*% t(res$vectors))
  
}


#' @export
#'
#' @title \code{invsqrt_psd_matrix}: inverse square root of a positive semi-definite matrix
#'
#' @description \code{invsqrt_psd_matrix} takes a square, positive semi-definite matrix and returns the inverse square root of that matrix (via the eigenvalue decomposition by way of \code{tolerance_eigen}).
#'
#' @details Note that \code{invsqrt_psd_matrix} uses \code{tolerance_eigen(...,tol=1e-13)} with the tolerance fixed at 1e-13. This enforces an "acceptably zero" value for eigenvalues, and also stops the computation if any eigenvalues are not zero or positive (within the tolerance parameter). If the computation is halted, that means the matrix was not positive semi-definite.
#'
#' @param x A square matrix (presumably positive semi-definite)
#'
#' @return A matrix. The inverse square root of the \code{x} matrix
#' @seealso \code{\link{tolerance_eigen}}

invsqrt_psd_matrix <- function(x){
  
  ## checks: just that they are a matrix & square & numeric
  if( length(dim(x)) != 2 ){
    stop("invsqrt_psd_matrix: x is not a matrix.")
  }
  if( !is.numeric(x) ){
    stop("invsqrt_psd_matrix: x is not numeric.")
  }
  if( dim(x)[1] != dim(x)[2] ){
    stop("invsqrt_psd_matrix: x is not a square matrix.")
  }
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }
  
  ## tolerance_eigen
  res <- tolerance_eigen(x, tol = 1e-13)
  
  ## rebuild
  return(t(t(res$vectors) * (1/sqrt(res$values)) ) %*% t(res$vectors))
  
}

