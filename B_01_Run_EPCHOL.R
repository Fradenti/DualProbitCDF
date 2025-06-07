source("A_01_1_aux_functions.R")
# source("A_01_2_Create_COVmat.R")
source("A_02_algos.R")

Ns   <- c(16, 64, 128, 256, 512, 1024)
RHOs <- c(0, .25, .5, .75)
NSIM <- 10
Bs   <- seq(-2,2,length.out=20)


sink("RDS/EP_CHOL/progress_EPChol.txt")

fixedCOVs <- readRDS("RDS/00_all_covs_rhos.RDS")
# CONSTANT CASE -------------------------------------------------
for(j in 1:length(RHOs)){
  for(i in 1:length(Ns)){
    PAR <- c()
    COV <- fixedCOVs[[i]][,,j]
    cat("---> Starting with dimension",Ns[i]," and correlation",RHOs[j],"\n")  
    for(b in 1:length(Bs)){
      
      for(k in 1:NSIM){
        set.seed(123*k)
        cat(paste("Run:", k, "out of", NSIM," - ", 
                  (Sys.time()), "with b:", round(Bs[b],3),"\n"))
        PAR <- rbind(PAR,
                     EP_CHOL(covM = COV, 
                             b_vec = rep(Bs[b],Ns[i]), 
                             eps = 100, tol = 1e-4, 
                             type = "const")
        )
        cat(paste("--------------------------------------\n"))
      }
    }
    nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
                 Ns[i],"_rho_", RHOs[j],"_b_fixed_A2.RDS")
    saveRDS(PAR, nam)
    cat("---> Finished with dimension",Ns[i]," and correlation",RHOs[j],"\n")  
  }
}

# Fungible CASE -------------------------------------------------
fungiCOVs <- readRDS("RDS/00_all_covs_fungi.RDS")
for(i in 1:length(Ns)){
  PAR = c()
  COV <- fungiCOVs[[i]]
  cat("---> Starting with dimension",Ns[i]," and correlation",RHOs[j],"\n")  
  for(b in 1:length(Bs)){
    for(k in 1:NSIM){
      set.seed(123*k)
      cat(paste("Run:", k, "out of", NSIM," - ", 
                (Sys.time()), "with b:", round(Bs[b],3),"\n"))
      PAR <- rbind(PAR,
                   EP_CHOL(covM = COV, 
                          b_vec = rep(Bs[b],Ns[i]), 
                          eps = 100, tol = 1e-4, type = "fungible")
      )
      cat(paste("--------------------------------------\n"))
    }
  }
  nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
               Ns[i],"_fungible_b_fixed_A2.RDS")
  saveRDS(PAR, nam)
    cat("---> Finished with dimension",Ns[i],"\n")  
}


# Dense CASE -------------------------------------------------
denseCOVs <- readRDS("RDS/00_all_covs_dense.RDS")
for(i in 1:length(Ns)){
  PAR = c()
  COV <- denseCOVs[[i]]
  cat("---> Starting with dimension",Ns[i],"\n")  
  for(b in 1:length(Bs)){
    for(k in 1:NSIM){
      set.seed(123*k)
      cat(paste("Run:", k, "out of", NSIM," - ", 
                (Sys.time()), "with b:", round(Bs[b],3),"\n"))
      PAR <- rbind(PAR,
                   EP_CHOL(covM = COV, 
                          b_vec = rep(Bs[b],Ns[i]), 
                          eps = 100, tol = 1e-4, type = "dense")      )
      cat(paste("--------------------------------------\n"))
    }
  }
  nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
               Ns[i],"_dense_XtX_b_fixed_A2.RDS")
  saveRDS(PAR, nam)
  cat("---> Finished with dimension",Ns[i],"\n")  
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

source("A_01_1_aux_functions.R")
# source("A_01_2_Create_COVmat.R")
source("A_02_algos.R")

Ns   <- c(16, 64, 128, 256, 512, 1024)
RHOs <- c(0, .25, .5, .75)
NSIM <- 10
Bs   <- seq(-2,2,length.out=20)


fixedCOVs <- readRDS("RDS/00_all_covs_rhos.RDS")
# CONSTANT CASE -------------------------------------------------
for(j in 1:length(RHOs)){
  for(i in 1:length(Ns)){
    PAR <- c()
    COV <- fixedCOVs[[i]][,,j]
    cat("---> Starting with dimension",Ns[i]," and correlation",RHOs[j],"\n")  
    for(b in 1:length(Bs)){
      
      for(k in 1:NSIM){
        set.seed(123*k)
        cat(paste("Run:", k, "out of", NSIM," - ", 
                  (Sys.time()), "with b:", round(Bs[b],3),"\n"))
        PAR <- rbind(PAR,
                     EP_CHOL_algo1(covM = COV, 
                             b_vec = rep(Bs[b],Ns[i]), 
                             eps = 100, tol = 1e-4, 
                             type = "const")
        )
        cat(paste("--------------------------------------\n"))
      }
    }
    nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
                 Ns[i],"_rho_", RHOs[j],"_b_fixed_A1.RDS")
    saveRDS(PAR, nam)
    cat("---> Finished with dimension",Ns[i]," and correlation",RHOs[j],"\n")  
  }
}

fungiCOVs <- readRDS("RDS/00_all_covs_fungi.RDS")
for(i in 1:length(Ns)){
  PAR = c()
  COV <- fungiCOVs[[i]]
  cat(paste("Run:", k, "out of", NSIM," - ", 
            (Sys.time()), "with b:", round(Bs[b],3),"\n"))
  for(b in 1:length(Bs)){
    for(k in 1:NSIM){
      set.seed(123*k)
      PAR <- rbind(PAR,
                   EP_CHOL_algo1(covM = COV, 
                           b_vec = rep(Bs[b],Ns[i]), 
                           eps = 100, tol = 1e-4, type = "fungible")
      )
      cat(paste("--------------------------------------\n"))
    }
  }
  nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
               Ns[i],"_fungible_b_fixed_A1.RDS")
  saveRDS(PAR, nam)
  cat("---> Finished with dimension",Ns[i],"\n")  
}


denseCOVs <- readRDS("RDS/00_all_covs_dense.RDS")
for(i in 1:length(Ns)){
  PAR = c()
  COV <- denseCOVs[[i]]
  cat("---> Starting with dimension",Ns[i],"\n")  
  for(b in 1:length(Bs)){
    for(k in 1:NSIM){
      set.seed(123*k)
      PAR <- rbind(PAR,
                   EP_CHOL_algo1(covM = COV, 
                           b_vec = rep(Bs[b],Ns[i]), 
                           eps = 100, tol = 1e-4, type = "dense")      )
      cat(paste("--------------------------------------\n"))
    }
  }
  nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
               Ns[i],"_dense_XtX_b_fixed_A1.RDS")
  saveRDS(PAR, nam)
  cat("---> Finished with dimension",Ns[i],"\n")  
}

sink()
