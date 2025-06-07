# TRUE PROB CONSTANT CASE -------------------------------------------------
TRUEp <- c()
for(j in 1:length(RHOs)){
  for(i in 1:length(Ns)){
    for(b in 1:length(Bs)){
      
      if( RHOs[j] == 0 ){
       prob <- Ns[i] * pnorm(Bs[b], 0, 1,log.p = TRUE)      
      }else{
       prob <- log(constRhoProb(b = rep(Bs[b],Ns[i]),rho = RHOs[j]))
      }
      TRUEp <- rbind(TRUEp,
                     c("rho" = RHOs[j], "Ns" = Ns[i], "b" = Bs[b], "Truep" = prob ))
      
      cat(paste("b",b,"--- i",i,"--- j",j,"\n"))
    }
  }
}
TRUEp <- data.frame(TRUEp) 
TRUEp %>% filter(Truep == 0)

plot(TRUEp$Truep)

# saveRDS(TRUEp, "RDS/00_trueLOGProbs_constantSigma.RDS")

fixedCOVs <- list()
for(i in 1:length(Ns)){
 covs <- array(NA,c(Ns[i],Ns[i],length(RHOs)))
  for(j in 1:length(RHOs)){
  covs[,,j]  <- create_covM(n = Ns[i], rho = RHOs[j])
}
 fixedCOVs[[i]] <- covs
}
fixedCOVs[[1]]
saveRDS(fixedCOVs, "RDS/00_all_covs_rhos.RDS")


fungiCOVs <- list()
for(i in 1:length(Ns)){
  fungiCOVs[[i]]  <- create_covM(n = Ns[i], seed = 2508*i, type = "fungible")
}
fixedCOVs[[1]]
saveRDS(fungiCOVs, "RDS/00_all_covs_fungi.RDS")

denseCOVs <- list()
for(i in 1:length(Ns)){
  denseCOVs[[i]]  <- create_covM(n = Ns[i], seed = 2508*i, type = "dense_XtX")
}
denseCOVs[[1]]
saveRDS(denseCOVs, "RDS/00_all_covs_dense.RDS")