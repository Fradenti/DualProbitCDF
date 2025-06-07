library(tidyverse)
Ns   = c(16, 64, 128, 256, 512, 1024)
RHOs = c(0,.25,.5,.75)
NSIM  = 10
Bs = seq(-2,2,length.out=20)

theme_set(theme_bw())

i=j=1

# Constant COR matrix ----------------------------------------------------------

D_TRUNCN = c()
for(i in seq_along(Ns)){
  for(j in seq_along(RHOs)){
  nam = paste0("RDS/TRUNCN/TRUNCN_dim_",
                Ns[i],"_rho_", RHOs[j],"_b_fixed.RDS")
  dat <-  readRDS(nam) %>% mutate(rho = RHOs[j])
  D_TRUNCN <- D_TRUNCN %>% bind_rows(dat)
  }
  }
D_TRUNCN  

D_TLRGenz = c()
for(i in seq_along(Ns)){
  for(j in seq_along(RHOs)){
    nam = paste0("RDS/TLRGenz/TLRGenz_dim_",
                 Ns[i],"_rho_", RHOs[j],"_b_fixed.RDS")
    dat <-  readRDS(nam) %>% mutate(rho = RHOs[j], dim = Ns[i])
    D_TLRGenz <- D_TLRGenz %>% bind_rows(dat)
  }
}
D_TLRGenz  

D_EPeig = c()
for(i in seq_along(Ns)){
  for(j in seq_along(RHOs)){
    nam = paste0("RDS/EP_EIG/EP_EIG_dim_",
                 Ns[i],"_rho_", RHOs[j],"_b_fixed_A2.RDS")
    dat <-  readRDS(nam) %>% mutate(rho = RHOs[j])
    D_EPeig <- D_EPeig %>% bind_rows(dat)
  }
}
D_EPeig


D_EPeig_A1 = c()
for(i in seq_along(Ns)){
  for(j in seq_along(RHOs)){
    nam = paste0("RDS/EP_EIG/EP_EIG_dim_",
                 Ns[i],"_rho_", RHOs[j],"_b_fixed_A1.RDS")
    dat <-  readRDS(nam) %>% mutate(rho = RHOs[j])
    D_EPeig_A1 <- D_EPeig_A1 %>% bind_rows(dat)
  }
}
D_EPeig_A1


D_EPchol = c()
for(i in seq_along(Ns)){
  for(j in seq_along(RHOs)){
    nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
                 Ns[i],"_rho_", RHOs[j],"_b_fixed_A2.RDS")
    dat <-  readRDS(nam) %>% mutate(rho = RHOs[j])
    D_EPchol <- D_EPchol %>% bind_rows(dat)
  }
}
D_EPchol


D_EPchol_A1 = c()
for(i in seq_along(Ns)){
  for(j in seq_along(RHOs)){
    nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
                 Ns[i],"_rho_", RHOs[j],"_b_fixed_A1.RDS")
    dat <-  readRDS(nam) %>% mutate( rho = RHOs[j])
    D_EPchol_A1 <- D_EPchol_A1 %>% bind_rows(dat)
  }
}
D_EPchol_A1

D_ORTH = c()
for(i in seq_along(Ns)[-6]){
  for(j in seq_along(RHOs)){
    nam = paste0("RDS/ORTH/ORTH_dim_",
                 Ns[i],"_rho_",RHOs[j],"_b_fixed.RDS")
    dat <-  readRDS(nam) %>% mutate(dim=Ns[i], info= "none", rho = RHOs[j])
    D_ORTH <- D_ORTH %>% bind_rows(dat)
  }
}
D_ORTH

D_TRUNCN$info <- as.character(D_TRUNCN$info)
D_rhofix <- bind_rows(D_TRUNCN,D_TLRGenz,D_EPeig,D_EPeig_A1,D_EPchol,D_EPchol_A1,D_ORTH)
saveRDS(D_rhofix, "RDS/ALL_rhofixed_2025.RDS")


# FUNGIBLE ----------------------------------------------------------------
D_TRUNCN = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/TRUNCN/TRUNCN_dim_",
                 Ns[i],"_fungible_b_fixed.RDS")
    dat <-  readRDS(nam) 
    D_TRUNCN <- D_TRUNCN %>% bind_rows(dat)
  }
D_TRUNCN  

D_TLRGenz = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/TLRGenz/TLRGenz_dim_",
                 Ns[i],"_fungible_b_fixed.RDS")
    dat <-  readRDS(nam) 
    D_TLRGenz <- D_TLRGenz %>% bind_rows(dat)
  }
D_TLRGenz  

D_EPeig = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/EP_EIG/EP_EIG_dim_",
                 Ns[i],"_fungible_b_fixed_A2.RDS")
    dat <-  readRDS(nam) 
    D_EPeig <- D_EPeig %>% bind_rows(dat)
}
D_EPeig


D_EPeig_A1 = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/EP_EIG/EP_EIG_dim_",
                 Ns[i],"_fungible_b_fixed_A1.RDS")
    dat <-  readRDS(nam) 
    D_EPeig_A1 <- D_EPeig_A1 %>% bind_rows(dat)
}
D_EPeig_A1


D_EPchol = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
                 Ns[i],"_fungible_b_fixed_A2.RDS")
    dat <-  readRDS(nam)
    D_EPchol <- D_EPchol %>% bind_rows(dat)
}
D_EPchol


D_EPchol_A1 = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
                 Ns[i],"_fungible_b_fixed_A1.RDS")
    dat <-  readRDS(nam) 
    D_EPchol_A1 <- D_EPchol_A1 %>% bind_rows(dat)
  }
D_EPchol_A1

D_ORTH = c()
for(i in seq_along(Ns)[-6]){
    nam = paste0("RDS/ORTH/ORTH_dim_",
                 Ns[i],"_fungible_b_fixed.RDS")
    dat <-  readRDS(nam) %>% mutate(dim=Ns[i], info= "none")
    D_ORTH <- D_ORTH %>% bind_rows(dat)
}
D_ORTH

D_TRUNCN$info <- as.character(D_TRUNCN$info)
D_fungible <- bind_rows(D_TRUNCN,D_TLRGenz,D_EPeig,D_EPeig_A1,D_EPchol,D_EPchol_A1,D_ORTH)
saveRDS(D_fungible, "RDS/ALL_fungible_2025.RDS")


# DENSE XtX ----------------------------------------------------------------
D_TRUNCN = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/TRUNCN/TRUNCN_dim_",
                 Ns[i],"_dense_XtX_b_fixed.RDS")
    dat <-  readRDS(nam)
    D_TRUNCN <- D_TRUNCN %>% bind_rows(dat)
}
D_TRUNCN  

D_TLRGenz = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/TLRGenz/TLRGenz_dim_",
                 Ns[i],"_dense_XtX_b_fixed.RDS")
    dat <-  readRDS(nam)
    D_TLRGenz <- D_TLRGenz %>% bind_rows(dat)
}
D_TLRGenz  

D_EPeig = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/EP_EIG/EP_EIG_dim_",
                 Ns[i],"_dense_XtX_b_fixed_A2.RDS")
    dat <-  readRDS(nam) 
    D_EPeig <- D_EPeig %>% bind_rows(dat)
  }
D_EPeig


D_EPeig_A1 = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/EP_EIG/EP_EIG_dim_",
                 Ns[i],"_dense_XtX_b_fixed_A1.RDS")
    dat <-  readRDS(nam)
    D_EPeig_A1 <- D_EPeig_A1 %>% bind_rows(dat)
  }
D_EPeig_A1


D_EPchol = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
                 Ns[i],"_dense_XtX_b_fixed_A2.RDS")
    dat <-  readRDS(nam)
    D_EPchol <- D_EPchol %>% bind_rows(dat)
  }
D_EPchol


D_EPchol_A1 = c()
for(i in seq_along(Ns)){
    nam = paste0("RDS/EP_CHOL/EP_CHOL_dim_",
                 Ns[i],"_dense_XtX_b_fixed_A1.RDS")
    dat <-  readRDS(nam) 
    D_EPchol_A1 <- D_EPchol_A1 %>% bind_rows(dat)
  }
D_EPchol_A1


D_ORTH = c()
for(i in seq_along(Ns)[-6]){
    nam = paste0("RDS/ORTH/ORTH_dim_",
                 Ns[i],"_dense_XtX_b_fixed.RDS")
    dat <-  readRDS(nam) %>% mutate(dim=Ns[i], info= "none")
    D_ORTH <- D_ORTH %>% bind_rows(dat)
  }
D_ORTH


D_TRUNCN$info <- as.character(D_TRUNCN$info)
D_dense_XtX <- bind_rows(D_TRUNCN,D_TLRGenz,D_EPeig,D_EPeig_A1,D_EPchol,D_EPchol_A1,D_ORTH)
saveRDS(D_dense_XtX, "RDS/ALL_dense_XtX_2025.RDS")


