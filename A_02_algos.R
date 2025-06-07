
# ONE-TIME RUNS -----------------------------------------------------

# Truncated Normal -------------------------------------------------------------
TRUNCN <- function(covM, b_vec, type){
  n          <- nrow(covM)
  startTime2 <- Sys.time()
  trnc       <- TruncatedNormal::pmvnorm(mu = rep(0, n), 
                                         sigma = covM, 
                                         ub = b_vec)
  time   <- difftime(Sys.time(), startTime2, units=("secs"))[[1]]
  prob   <- trnc[[1]]
  relerr <- attr(trnc, "relerr")
  
  Results   <-  data.frame(log2prob = log2(prob), 
                           time = time, 
                           type = type,
                           algo = "Botev",
                           b = b_vec[1],
                           dim = n,
                           info = relerr) # relative error
  return(Results)
}

## TileLowRank or Genz  --------------------------------------------------------
TLRGenz <- function(covM, b_vec, type){
  n <- nrow(covM)
  
  if(type  != "const"){
    startTime1 <- Sys.time()
    prob       <- tlrmvnmvt::pmvn(lower = rep(-Inf,  n), 
                             upper = b_vec, 
                             sigma = covM)[1]
    time       <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
    
    R1 <- data.frame(log2prob = log2(prob), 
                     time = time, 
                     type = type,
                     algo = "Genz",
                     b = b_vec[1],
                     dim = n,
                     info = "none")
    
  }else{
    startTime1 <- Sys.time()
    prob    <- tlrmvnmvt::pmvn(lower = rep(-Inf,  n), 
                               upper = b_vec, 
                               sigma = covM, 
                               algorithm = tlrmvnmvt::TLRQMC(m = round(sqrt(n))))
    time <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
    R1   <- data.frame(log2prob = log2(prob), 
                       time = time, 
                       type = type,
                       algo = "TLRank",
                       b = b_vec[1],
                       dim = n,
                       info = "none")
  }
  
  return(R1)
}

## EP Chol Algo2  --------------------------------------------------------
EP_CHOL = function(covM, b_vec, eps, tol, type){
  n <- nrow(covM)
  startTime3 <- Sys.time()
  ep_prob  <- EPmvnCDF::FASt_cdf(x = b_vec,
                                 Sigma = covM,
                                 log.p = TRUE,
                                 eps = eps,
                                 tolerance = tol, 
                                 algorithm = "B",
                                 method = "chol")
  tim3 <- difftime(Sys.time(), startTime3, units=("secs"))[[1]]
  
  R3 = data.frame(
    log2prob =ep_prob/log(2), 
    time = tim3, 
    type = type,
    algo = "EP_Chol_algoB",
    b = b_vec[1],
    dim = n,
    info = paste0("eps=",eps," tol=",tol))
  return(R3) 
}

## EP Eig Algo2  --------------------------------------------------------
EP_EIG = function(covM, b_vec, eps, tol, type){
  n <- nrow(covM)
  startTime3 <- Sys.time()
  ep_prob_eig <- EPmvnCDF::FASt_cdf(x = b_vec,
                                    Sigma = covM,
                                    log.p = TRUE,
                                    eps = eps,
                                    tolerance = tol, 
                                    algorithm = "B",
                                    method = "eig")
  tim3 <- difftime(Sys.time(), startTime3, units=("secs"))[[1]]
  
  R3 = data.frame(
    log2prob =ep_prob_eig/log(2), 
    time = tim3, 
    type = type,
    algo = "EP_Eigen_algoB",
    b = b_vec[1],
    dim = n,
    info = paste0("eps=",eps," tol=",tol))
  return(R3)  
}

## EP Chol Algo1  --------------------------------------------------------
EP_CHOL_algo1 = function(covM, b_vec, eps, tol, type){
  n <- nrow(covM)
  startTime3 <- Sys.time()
  ep_prob_eig <- EPmvnCDF::FASt_cdf(x = b_vec,
                                    Sigma = covM,
                                    log.p = TRUE,
                                    eps = eps,
                                    tolerance = tol, 
                                    algorithm = "A",
                                    method = "chol")
  tim3 <- difftime(Sys.time(), startTime3, units=("secs"))[[1]]
  
  R3 = data.frame(
    log2prob =ep_prob_eig/log(2), 
    time = tim3, 
    type = type,
    algo = "EP_Chol_algoA",
    b = b_vec[1],
    dim = n,
    info = paste0("eps=",eps," tol=",tol))
  return(R3)  
}

## EP Eig Algo1  --------------------------------------------------------
EP_EIG_algo1 = function(covM, b_vec, eps, tol, type){
  n <- nrow(covM)
  startTime3 <- Sys.time()
  ep_prob_eig <- EPmvnCDF::FASt_cdf(x = b_vec,
                                             Sigma = covM,
                                             log.p = TRUE,
                                             eps = eps,
                                             tolerance = tol,
                                             algorithm = "A",
                                             method = "eig")
  tim3 <- difftime(Sys.time(), startTime3, units=("secs"))[[1]]
  
  R3 = data.frame(
    log2prob =ep_prob_eig/log(2), 
    time = tim3, 
    type = type,
    algo = "EP_Eigen_algoA",
    b = b_vec[1],
    dim = n,
    info = paste0("eps=",eps," tol=",tol))
  return(R3)  
}

## EP Eig Algo1  --------------------------------------------------------
Orthants_Ridgeway = function(covM, b_vec, M = 1e4, verb = F, ESS_fraction = .5, type){
  n <- nrow(covM)
  startTime0 <- Sys.time()
  orthant_p <- ridgway_smc_cpp(x = b_vec,
                               Sigma = covM,
                               ESS_fraction = ESS_fraction,
                               M = M, verbose = verb)
  tim3 <- difftime(Sys.time(), startTime0, units=("secs"))[[1]]
  
  R3 = data.frame(
    log2prob = orthant_p$logCDF/log(2), 
    time = tim3, 
    type = type, #paste0("M=", M),
    algo = "RIDGEWAY",
    b = b_vec[1],
    dim = n,
    info =paste0("M = ", M,"-ESSf = ",ESS_fraction))
  return(R3)  
}
