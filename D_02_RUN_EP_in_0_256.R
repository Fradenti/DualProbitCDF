Rcpp::sourceCpp("cpp_source/ridgeway_cpp.cpp")
cov_dense <- readRDS("RDS/00_all_covs_dense.RDS")
cov256 <- cov_dense[[4]]
n <- nrow(cov256)
NSIM <- 10
u    <- 0

set.seed(123)
t1 = Sys.time()
ep_prob  <- EPmvnCDF::FASt_cdf(x = rep(u,n),
                               Sigma = cov256,
                               log.p = TRUE,
                               eps = 100, 
                               tol = 1e-4,
                               algorithm = "B",
                               method = "chol")
timeEP <- difftime(Sys.time(), t1, units=("secs"))[[1]]
timeEP
ep_prob

saveRDS((ep_prob)/log(2), "RDS/Comparison_in_0/Dense_EPCHOL2_log2P_256.RDS")
saveRDS(timeEP, "RDS/Comparison_in_0/Dense_EPCHOL2_times_256.RDS")

cov_fung <- readRDS("RDS/00_all_covs_fungi.RDS")
cov256 <- cov_fung[[4]]
n <- nrow(cov256)
NSIM <- 10
set.seed(123)
t1 = Sys.time()
ep_prob  <- EPmvnCDF::FASt_cdf(x = rep(u,n),
                               Sigma = cov256,
                               log.p = FALSE,
                               eps = 100, tol = 1e-4,
                               algorithm = "B",
                               method = "chol")
timeEP <- difftime(Sys.time(), t1, units=("secs"))[[1]]
timeEP
abline(h = log2(ep_prob),col=2)
saveRDS(log2(ep_prob), "RDS/Comparison_in_0/fung_EPCHOL2_log2P_256.RDS")
saveRDS(timeEP, "RDS/Comparison_in_0/fung_EPCHOL2_times_256.RDS")




Rcpp::sourceCpp("cpp_source/ridgeway_cpp.cpp")
cov_rhos <- readRDS("RDS/00_all_covs_rhos.RDS")
cov256   <- cov_rhos[[4]][,,3]
n <- nrow(cov256)
NSIM <- 10
set.seed(123)
t1 = Sys.time()
ep_prob  <- EPmvnCDF::FASt_cdf(x = rep(u,n),
                               Sigma = cov256,
                               log.p = FALSE,
                               eps = 100, tol = 1e-4,
                               algorithm = "B",
                               method = "chol")
timeEP <- difftime(Sys.time(), t1, units=("secs"))[[1]]
timeEP
abline(h = log2(ep_prob),col=2)
saveRDS(log2(ep_prob), "RDS/Comparison_in_0/rho05_EPCHOL2_log2P_256.RDS")
saveRDS(timeEP, "RDS/Comparison_in_0/rho05_EPCHOL2_times_256.RDS")
