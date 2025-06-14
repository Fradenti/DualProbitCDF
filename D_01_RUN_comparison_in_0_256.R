# Dense ----------------------------------------------------------------

Rcpp::sourceCpp("cpp_source/ridgeway_cpp.cpp")
cov_dense <- readRDS("RDS/00_all_covs_dense.RDS")
cov256 <- cov_dense[[4]]
n <- nrow(cov256)
NSIM <- 10
u    <- 0

NSAMPLE <- c(5000,10000,20000,50000)
PROBS_GENZ <- TIMES_GENZ <- matrix(NA,NSIM,length(NSAMPLE))

for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime1 <- Sys.time()
    prob       <- log2(tlrmvnmvt::pmvn(lower = rep(-Inf,  n), 
                                upper = rep(u,n), 
                                sigma = cov256,
                                algorithm = tlrmvnmvt::GenzBretz(N = NSAMPLE[g]/20) )[1])
    time       <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
    TIMES_GENZ[k,g] <- time
    PROBS_GENZ[k,g]  <- prob
    cat(paste(k,"--",g,"\n"))
}
}
saveRDS(PROBS_GENZ, "RDS/Comparison_in_0/Dense_GENZ_log2P_256.RDS")
saveRDS(TIMES_GENZ, "RDS/Comparison_in_0/Dense_GENZ_times_256.RDS")

PROBS_BOTEV <- TIMES_BOTEV <- matrix(NA,NSIM,length(NSAMPLE))
for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
        set.seed(123*k)
        startTime1 <- Sys.time()
        prob       <-   log2(
                            TruncatedNormal::pmvnorm(mu = rep(0, n), 
                                                     sigma = cov256, 
                                                     ub = rep(u,n),
                                                     B = NSAMPLE[g])
                          )
        time       <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
        TIMES_BOTEV[k,g] <- time
        PROBS_BOTEV[k,g]  <- prob
        cat(paste(k,"--",g,"\n"))
      }
}
saveRDS(PROBS_BOTEV, "RDS/Comparison_in_0/Dense_BOTEV_log2P_256.RDS")
saveRDS(TIMES_BOTEV, "RDS/Comparison_in_0/Dense_BOTEV_times_256.RDS")


PROBS_RIDGE <- TIMES_RIDGE <- matrix(NA,NSIM,length(NSAMPLE))
for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime0 <- Sys.time()
    orthant_p <- ridgway_smc_cpp(x = rep(u,n),
                                 Sigma = cov256,
                                 ESS_fraction = .5,
                                 M = NSAMPLE[g], 
                                 verbose = FALSE)
    time <- difftime(Sys.time(), startTime0, units=("secs"))[[1]]
    TIMES_RIDGE[k,g] <- time
    PROBS_RIDGE[k,g]  <- (orthant_p$logCDF)/log(2)
    cat(paste(k,"--",g,"\n"))
  }
}
saveRDS(PROBS_RIDGE, "RDS/Comparison_in_0/Dense_RIDGE_log2P_256.RDS")
saveRDS(TIMES_RIDGE, "RDS/Comparison_in_0/Dense_RIDGE_times_256.RDS")

# Fungible ----------------------------------------------------------------

Rcpp::sourceCpp("cpp_source/ridgeway_cpp.cpp")
cov_fung <- readRDS("RDS/00_all_covs_fungi.RDS")
cov256 <- cov_fung[[4]]
n <- nrow(cov256)
NSIM <- 10

NSAMPLE <- c(5000,10000,20000,50000)
PROBS_GENZ <- TIMES_GENZ <- matrix(NA,NSIM,length(NSAMPLE))

for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime1 <- Sys.time()
    prob       <- log2(tlrmvnmvt::pmvn(lower = rep(-Inf,  n), 
                                       upper = rep(u,n), 
                                       sigma = cov256,
                                       algorithm = tlrmvnmvt::GenzBretz(N = NSAMPLE[g]/20))[1])
    time       <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
    TIMES_GENZ[k,g] <- time
    PROBS_GENZ[k,g]  <- prob
    cat(paste(k,"--",g,"\n"))
  }
}
saveRDS(PROBS_GENZ, "RDS/Comparison_in_0/fung_GENZ_log2P_256.RDS")
saveRDS(TIMES_GENZ, "RDS/Comparison_in_0/fung_GENZ_times_256.RDS")

PROBS_BOTEV <- TIMES_BOTEV <- matrix(NA,NSIM,length(NSAMPLE))
for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime1 <- Sys.time()
    prob       <-   log2(
      TruncatedNormal::pmvnorm(mu = rep(0, n), 
                               sigma = cov256, 
                               ub = rep(u,n),
                               B = NSAMPLE[g])
    )
    time       <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
    TIMES_BOTEV[k,g] <- time
    PROBS_BOTEV[k,g]  <- prob
    cat(paste(k,"--",g,"\n"))
  }
}
saveRDS(PROBS_BOTEV, "RDS/Comparison_in_0/fung_BOTEV_log2P_256.RDS")
saveRDS(TIMES_BOTEV, "RDS/Comparison_in_0/fung_BOTEV_times_256.RDS")


PROBS_RIDGE <- TIMES_RIDGE <- matrix(NA,NSIM,length(NSAMPLE))
for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime0 <- Sys.time()
    orthant_p <- ridgway_smc_cpp(x = rep(u,256),
                                 Sigma = cov256,
                                 ESS_fraction = .5,
                                 M = NSAMPLE[g], 
                                 verbose = FALSE)
    time <- difftime(Sys.time(), startTime0, units=("secs"))[[1]]
    TIMES_RIDGE[k,g] <- time
    PROBS_RIDGE[k,g]  <- (orthant_p$logCDF)/log(2)
    cat(paste(k,"--",g,"\n"))
  }
}
saveRDS(PROBS_RIDGE, "RDS/Comparison_in_0/fung_RIDGE_log2P_256.RDS")
saveRDS(TIMES_RIDGE, "RDS/Comparison_in_0/fung_RIDGE_times_256.RDS")
saveRDS(log2(ep_prob), "RDS/Comparison_in_0/Dense_EPCHOL2_log2P_256.RDS")


# Cov50 ----------------------------------------------------------------

Rcpp::sourceCpp("cpp_source/ridgeway_cpp.cpp")
cov_rhos <- readRDS("RDS/00_all_covs_rhos.RDS")
cov256   <- cov_rhos[[4]][,,3]
n <- nrow(cov256)
NSIM <- 10

NSAMPLE <- c(5000,10000,20000,50000)
PROBS_GENZ <- TIMES_GENZ <- matrix(NA,NSIM,length(NSAMPLE))

for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime1 <- Sys.time()
    prob       <- log2(tlrmvnmvt::pmvn(lower = rep(-Inf,  n), 
                                       upper = rep(u,n), 
                                       sigma = cov256,
                                       algorithm = tlrmvnmvt::TLRQMC(N = NSAMPLE[g]/20,
                                                                     m =round(sqrt(n)))))
    time       <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
    TIMES_GENZ[k,g] <- time
    PROBS_GENZ[k,g]  <- prob
    cat(paste(k,"--",g,"\n"))
  }
}
saveRDS(PROBS_GENZ, "RDS/Comparison_in_0/rho05_TLRNK_log2P_256.RDS")
saveRDS(TIMES_GENZ, "RDS/Comparison_in_0/rho05_TLRNK_times_256.RDS")

PROBS_BOTEV <- TIMES_BOTEV <- matrix(NA,NSIM,length(NSAMPLE))
for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime1 <- Sys.time()
    prob       <-   log2(
      TruncatedNormal::pmvnorm(mu = rep(0, n), 
                               sigma = cov256, 
                               ub = rep(u,n),
                               B = NSAMPLE[g])
    )
    time       <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
    TIMES_BOTEV[k,g] <- time
    PROBS_BOTEV[k,g]  <- prob
    cat(paste(k,"--",g,"\n"))
  }
}
saveRDS(PROBS_BOTEV, "RDS/Comparison_in_0/rho05_BOTEV_log2P_256.RDS")
saveRDS(TIMES_BOTEV, "RDS/Comparison_in_0/rho05_BOTEV_times_256.RDS")


PROBS_RIDGE <- TIMES_RIDGE <- matrix(NA,NSIM,length(NSAMPLE))
for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime0 <- Sys.time()
    orthant_p <- ridgway_smc_cpp(x = rep(u,n),
                                 Sigma = cov256,
                                 ESS_fraction = .5,
                                 M = NSAMPLE[g], 
                                 verbose = FALSE)
    time <- difftime(Sys.time(), startTime0, units=("secs"))[[1]]
    TIMES_RIDGE[k,g] <- time
    PROBS_RIDGE[k,g]  <- (orthant_p$logCDF)/log(2)
    cat(paste(k,"--",g,"\n"))
  }
}
saveRDS(PROBS_RIDGE, "RDS/Comparison_in_0/rho05_RIDGE_log2P_256.RDS")
saveRDS(TIMES_RIDGE, "RDS/Comparison_in_0/rho05_RIDGE_times_256.RDS")



