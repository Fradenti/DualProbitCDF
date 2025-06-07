
# Loading libraries --------------------------------------------------------

# Hermite quadrature for true probabilities
library("statmod")
library("TruncatedNormal")
library("tlrmvnmvt")
library("mvtnorm")
library("tidyverse")
library("EPmvnCDF")
Rcpp::sourceCpp("cpp_source/ridgeway_cpp.cpp")

# Useful functions --------------------------------------------------------

## Constant correlation -----------------------------------------------------

nnode <- 200
nodeWeight <- gauss.quad(nnode, "hermite")
intfct <- function(x, b, rho) {
  y <- rep(0, length(x))
  for (i in 1:length(x))
    y[i] <- -.5 * log(pi) + sum(pnorm((b +
                                         sqrt(2 * rho) * x[i]) /
                                         sqrt(1 - rho), log.p = TRUE))
  return(exp(y))
}

constRhoProb <- function(b, rho) {
  sum(exp(log(nodeWeight$weights) + log(intfct(nodeWeight$nodes, b, rho))))
}

## Function to generate vcov matrices ------------------------------------------
create_covM = function(n,
                       rho = NULL,
                       seed = 123,
                       type = c("const", "dense_JSS", "dense_XtX", "fungible")) {
  type = match.arg(type)
  set.seed(23 * seed)
  
  if (type == "fungible") {
    lambda <- runif(n)
    lambda <- lambda * n / sum(lambda)
    covM   <- fungible::rGivens(lambda, Seed = i)$R
    
  } else if (type == "const") {
    covM       <- matrix(rho, n, n)
    diag(covM) <- 1
    
  } else if (type == "dense_JSS") {
    nx <- 20
    ny <- 20
    n <- nx * ny
    vecx <- c(1:nx) - 1
    vecy <- c(1:ny) - 1
    geom <- cbind(kronecker(vecx, rep(1, ny)), kronecker(rep(1, nx), vecy))
    geom <- geom + matrix(runif(n * 2), n, 2)
    geom <- geom / max(nx, ny)
    distM <- as.matrix(dist(geom)) # unfeasible if n >> 0
    covM <- matern(distM, 0.1, 1.0)
    
  } else if (type == "dense_XtX") {
    A = matrix(runif(n ^ 2) * 2 - 1, ncol = n)
    Sigma = t(A) %*% A
    covM <- cov2cor(Sigma)
  }
  
  return(covM)
  
}
