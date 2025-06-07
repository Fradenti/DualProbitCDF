Rcpp::sourceCpp("cpp_source/ridgeway_cpp.cpp")
cov_dense <- readRDS("RDS/00_all_covs_dense.RDS")
cov256 <- cov_dense[[4]]
n <- nrow(cov256)
NSIM <- 10

NSAMPLE <- c(500,1000,5000,10000, 25000)
PROBS_GENZ <- TIMES_GENZ <- matrix(NA,NSIM,length(NSAMPLE))

for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime1 <- Sys.time()
    prob       <- log2(tlrmvnmvt::pmvn(lower = rep(-Inf,  n), 
                                upper = rep(0,n), 
                                sigma = cov256,
                                algorithm = tlrmvnmvt::GenzBretz(N = NSAMPLE[g]))[1])
    time       <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
    TIMES_GENZ[k,g] <- time
    PROBS_GENZ[k,g]  <- prob
    cat(paste(k,"--",g,"\n"))
}
}
saveRDS(TIMES_GENZ, "RDS/Comparison_in_0/GENZ_log2P_256.RDS")
saveRDS(TIMES_GENZ, "RDS/Comparison_in_0/GENZ_times_256.RDS")

PROBS_BOTEV <- TIMES_BOTEV <- matrix(NA,NSIM,length(NSAMPLE))
for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
        set.seed(123*k)
        startTime1 <- Sys.time()
        prob       <-   log2(
                            TruncatedNormal::pmvnorm(mu = rep(0, n), 
                                                     sigma = cov256, 
                                                     ub = rep(0,n),
                                                     B = NSAMPLE[g])
                          )
        time       <- difftime(Sys.time(), startTime1, units=("secs"))[[1]]
        TIMES_BOTEV[k,g] <- time
        PROBS_BOTEV[k,g]  <- prob
        cat(paste(k,"--",g,"\n"))
      }
}
saveRDS(PROBS_BOTEV, "RDS/Comparison_in_0/BOTEV_log2P_256.RDS")
saveRDS(TIMES_BOTEV, "RDS/Comparison_in_0/BOTEV_times_256.RDS")


PROBS_RIDGE <- TIMES_RIDGE <- matrix(NA,NSIM,length(NSAMPLE))
for(g in seq_along(NSAMPLE)){
  for(k in 1:NSIM){
    set.seed(123*k)
    startTime0 <- Sys.time()
    orthant_p <- ridgway_smc_cpp(x = rep(0,256),
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
saveRDS(PROBS_RIDGE, "RDS/Comparison_in_0/RIDGE_log2P_256.RDS")
saveRDS(TIMES_RIDGE, "RDS/Comparison_in_0/RIDGE_times_256.RDS")



## La matrice di vcov rimane la stessa across u e simulationa
colnames(PROBS_BOTEV) = paste("N =",NSAMPLE)
colnames(PROBS_GENZ) = paste("N =",NSAMPLE)
colnames(PROBS_RIDGE) = paste("N =",NSAMPLE)

boxplot(PROBS_BOTEV,ylim=c(-320,-300))
boxplot(PROBS_GENZ,add=T)
boxplot(PROBS_RIDGE,add=T,fill=2)
ep_prob  <- EPmvnCDF::FASt_cdf(x = rep(0,n),
                               Sigma = cov256,
                               log.p = FALSE,
                               eps = 100, tol = 1e-4,
                               algorithm = "B",
                               method = "chol")
abline(h = log2(ep_prob),col=2)



## La matrice di vcov rimane la stessa across u e simulationa
colnames(TIMES_BOTEV) = paste("N =",NSAMPLE)
colnames(TIMES_GENZ) = paste("N =",NSAMPLE)
colnames(TIMES_RIDGE) = paste("N =",NSAMPLE)

boxplot(TIMES_BOTEV, ylim=c(0,1000))
boxplot(TIMES_GENZ,add=T)
boxplot(TIMES_RIDGE,add=T,fill=2)
ep_prob  <- EPmvnCDF::FASt_cdf(x = rep(0,n),
                               Sigma = cov256,
                               log.p = FALSE,
                               eps = 100, tol = 1e-4,
                               algorithm = "B",
                               method = "chol")
abline(h = log2(ep_prob),col=2)
