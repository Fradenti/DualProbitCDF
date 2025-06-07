#' FASt cdf algorithm
#'
#'@description This function computes the value of a generic cumulative distribution function of a multivariate Gaussian
#' exploiting its connection with the marginal likelihood of a tailored dual Bayesian probit model.
#'
#'@details
#' Let us consider, without loss of generality, an \eqn{m}-dimensional multivariate Gaussian random variable with mean 0 and covariance matrix \eqn{\Sigma}. Then, for \eqn{u = (u_1,\ldots,u_m)^\top\in \mathbb{R}^m}, its cumulative distribution function is defined as
#'
#' \deqn{\Phi_m(u;\Sigma) = \int_{{-\infty}}^{u} \int_{C_u}\dfrac{1}{\sqrt{(2\pi)^m|\Sigma|}} \exp \left\{-\dfrac{1}{2} w^\top \Sigma^{-1} w \right\} d w,}
#' where \eqn{C_u=\left\{w=(w_1,\ldots,w_m)^\top\in \mathbb{R}^m\colon w_i\le u_i,\ \forall i=1,\ldots,m \right\}}.
#'
#'
#' Given any positive definite \eqn{m\times m} covariance matrix \eqn{\Sigma}, call \eqn{\lambda_{m}} its smallest eigenvalue and define \eqn{\tilde{\Sigma}=\Sigma-\epsilon\lambda_{m}I_m}, where \eqn{\epsilon\in (0,1)}.
#' Consider now any factorization of \eqn{\tilde{\Sigma}} of the form \eqn{\tilde{\Sigma}=P\tilde\Lambda P^\top}, where \eqn{P} and \eqn{\tilde\Lambda} are \eqn{m\times m} matrices
#' with the former invertible and the latter diagonal and positive definite. Then, for any \eqn{u\in\mathbb{R}^m}, \eqn{\Phi_m(u;\Sigma)} \strong{equals the marginal likelihood of a dual Bayesian probit model}.
#'
#' Specifically, consider a generic observation vector \eqn{y\in\{0,1\}^n} and design matrix \eqn{X\in\mathbb{R}^{n\times p}}, with \eqn{i}-th row denoted with \eqn{x_i^\top}, indicating the covariate vector for observation \eqn{i}.
#' Then, this model can be defined as
#' \deqn{y_i\mid \beta\overset{ind}{\sim} \text{Bern}\left\{\Phi_1((2y_i-1)x_i^\top\beta;1)\right\}, \quad i=1,\ldots,n,\quad   \beta\sim\mathcal{N}_p\left({\xi},\Omega\right),}
#' where, in particular, \eqn{n=p=m}, \eqn{y_i=1} for all \eqn{i=1,\ldots,m,} and
#' \deqn{{\xi} =\left(\epsilon \lambda_{m}\right)^{-1/2}P^{-1}u,\quad \Omega= \left(\epsilon \lambda_{m}\right)^{-1}\tilde\Lambda,\quad X = P.}
#'
#' If the Gaussian of interest is not centered, one can assume \eqn{u = x-\mu}.
#'
#' The function implements several algorithms to estimate this dual Bayesian probit model via Expectation Propagation, and then computing its marginal likelihood.
#'
#' @param x the vector of quantiles of interest.
#' @param mu mean vector of the multivariate Gaussian of interest.
#' @param Sigma variance-Covariance matrix of the multivariate Gaussian of interest.
#' @param log.p logical. If \code{TRUE}, the log-cumulative probability is returned. Default is \code{FALSE}.
#' @param eps scalar, rescaling coefficient for the smallest eigenvalue (see Fasano and Denti, 2025+ for details).
#' @param tolerance scalar, a small positive value specifying the stopping criterion.
#' @param algorithm Can be either c("A","B"),
#' @param method Can be either c("chol","eigen")
#'
#' @return The estimated cumulative (log-) probability.
#' @export
#'
#'
#' @references Fasano, A., and Denti, F. (2025+). Multivariate Gaussian cumulative distribution functions as the marginal likelihood of their dual Bayesian probit models.
#' \emph{ArXiv}
#'
#'
#' @examples
#' FASt_cdf(c(0,0,0), log = FALSE)
FASt_cdf = function(x,
                    mu = rep(0, length(x)),
                    Sigma = diag(length(x)),
                    log.p = FALSE,
                    eps = 2,
                    tolerance = 1e-5,
                    algorithm = c("B", "A"),
                    method = c("chol", "eigen")) {
  method <- match.arg(method)
  algorithm <- match.arg(algorithm)
  
  # Checks ----------------------------------------------------------------
  
  if (length(mu) != length(x)) {
    stop(
      paste(
        "Error: the x provided has dimension",
        length(x),
        "while the mean vector has dimension",
        length(mu)
      ),
      "."
    )
  }
  if ((nrow(Sigma) != length(x)) | (ncol(Sigma) != length(x))) {
    stop(paste(
      "Error: Sigma has to be a squared matrix, while it has dimensions",
      dim(Sigma)
    ))
  }
  
  
  # Functions calls ----------------------------------------------------------------
 if (algorithm == "A") {
    if (method == "chol") {
      res <- FAStCDF_algo1(x = x, mu = mu, Sigma = Sigma,log.p = log.p, eps = eps, tolerance = tolerance)
    } else if (method == "eigen") {
      res <- FAStCDF_algo1_eig(x = x, mu = mu, Sigma = Sigma, log.p = log.p, eps = eps, tolerance = tolerance)
    }
    
  } else if (algorithm == "B") {
    
    if (method == "chol") {
      res <- FAStCDF_chol(x = x, mu = mu, Sigma = Sigma, log.p = log.p, eps = eps, tolerance = tolerance)
    } else if (method == "eigen") {
      res <- FAStCDF_eig(x = x, mu = mu, Sigma = Sigma, log.p = log.p, eps = eps, tolerance = tolerance)
    }
  }
  return(res)
}