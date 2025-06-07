#' FASt CDF computed via Algorithm 1 and Cholesky decomposition
#'
#' @param x a vector
#' @param mu mean vector
#' @param Sigma vcov matrix
#' @param log.p bool
#' @param eps rescale eigen
#' @param tolerance tolerance
#'
#' @return prob or log prob CDF
#'
#' @importFrom Rfast cholesky
#'
#'
#' @noRd
#' @examples
#' FAStCDF_algo1(c(0,0), log = FALSE )
FAStCDF_algo1 = function(x,
                       mu = rep(0, length(x)),
                       Sigma = diag(length(x)),
                       log.p = FALSE,
                       eps = 2,
                       tolerance = 1e-5) {
  
  
  
  # an attempt to circumvent degenerate scenarios
  ind <- which(is.infinite(x))
  if(length(ind)>0) {x   <- x[-ind]}
  n   <- length(x)
  if( n == 0 ){
    return(1.0)
  }
  
  x <- x - mu
  
  # define scaling
  scalingOfEig = eps # this could be 2, 10, etc...
  
  # find smallest eigenvector
  eigs = eigen(Sigma, symmetric = T, only.values = T)
  smallestEig = eigs$values[n]
  
  # find Choleski of tildeSigma
  tildeSigma = Sigma - smallestEig / scalingOfEig * diag(1, n, n)
  
  L = t(chol(tildeSigma))
  
  paramsEP = getParamsEP_priorMean_algo1_cpp(
    X = L,
    y = rep(1, n),
    priorMean = sqrt(scalingOfEig / smallestEig) *
      solve(L) %*% x,
    nu2 = scalingOfEig / smallestEig,
    tolerance = tolerance
  )
  
  log_p <- paramsEP$logML
  x <- ifelse(log.p, log_p, exp(log_p))
  
  return(x)
}


