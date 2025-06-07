
source("A_01_1_aux_functions.R")
# source("A_01_2_Create_COVmat.R")
source("A_02_algos.R")

Ns   <- c(16, 64, 128, 256, 512) #1024)
RHOs <- c(0, .25, .5, .75)
NSIM <- 10
Bs   <- seq(-2,2,length.out=20)

Ns[1]

sink("RDS/ORTH/progress_Orth.txt")

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
                     Orthants_Ridgeway(covM = COV, 
                                       b_vec = rep(Bs[b],Ns[i]), 
                                       ESS_fraction =  ifelse(Ns[i]>500, .1, .5),
                                       M = 1e4, type = "const"
                     )
        )
        cat(paste("--------------------------------------\n"))
      }
    }
    nam = paste0("RDS/ORTH/ORTH_dim_",
                 Ns[i],"_rho_", RHOs[j],"_b_fixed.RDS")
    saveRDS(PAR, nam)
    cat("---> Finished with dimension",Ns[i]," and correlation",RHOs[j],"\n")  
  }
}

fungiCOVs <- readRDS("RDS/00_all_covs_fungi.RDS")
for(i in 1:length(Ns)){
  PAR = c()
  COV <- fungiCOVs[[i]]
  cat("---> Starting with dimension",Ns[i],"\n")  
  for(b in 1:length(Bs)){
    for(k in 1:NSIM){
      set.seed(123*k)
      cat(paste("Run:", k, "out of", NSIM," - ", 
                (Sys.time()), "with b:", round(Bs[b],3),"\n"))
      PAR <- rbind(PAR,
                   Orthants_Ridgeway(covM = COV, 
                                     b_vec = rep(Bs[b],Ns[i]), 
                                     ESS_fraction =  ifelse(Ns[i]>500, .1, .5),
                                     M = 1e4, type = "fungible"
                   )
      )
      cat(paste("--------------------------------------\n"))
    }
  }
  nam = paste0("RDS/ORTH/ORTH_dim_",
               Ns[i],"_fungible_b_fixed.RDS")
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
      cat(paste("Run:", k, "out of", NSIM," - ", 
                (Sys.time()), "with b:", round(Bs[b],3),"\n"))
      PAR <- rbind(PAR,
                   Orthants_Ridgeway(covM = COV, 
                                     b_vec = rep(Bs[b],Ns[i]), 
                                     ESS_fraction = ifelse(Ns[i]>500, .1, .5),
                                     M = 1e4, type = "dense")
      )
      cat(paste("--------------------------------------\n"))
    }
  }
  nam = paste0("RDS/ORTH/ORTH_dim_",
               Ns[i],"_dense_XtX_b_fixed.RDS")
  saveRDS(PAR, nam)
  cat("---> Finished with dimension",Ns[i],"\n")  
}


sink()
