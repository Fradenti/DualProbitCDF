library(tidyverse)
# extract_results ---------------------------------------------------------



NSAMPLE <- paste("N =",c(500,1000,5000,10000, 25000))
PROBS_BOTEV <- readRDS("RDS/Comparison_in_0/rho05_BOTEV_log2P_256.RDS")
PROBS_GENZ  <- readRDS("RDS/Comparison_in_0/rho05_TLRNK_log2P_256.RDS")
PROBS_RIDGE <- readRDS("RDS/Comparison_in_0/rho05_RIDGE_log2P_256.RDS")

colnames(PROBS_BOTEV) <- colnames(PROBS_GENZ) <- colnames(PROBS_RIDGE) <- NSAMPLE

a1 = PROBS_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = PROBS_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = PROBS_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_rho <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Const.")
ep_rho <- readRDS("RDS/Comparison_in_0/rho05_EPCHOL2_log2P_256.RDS")


# -------------------------------------------------------------------------


NSAMPLE <- paste("N =",c(500,1000,5000,10000, 25000))
PROBS_BOTEV <- readRDS("RDS/Comparison_in_0/fung_BOTEV_log2P_256.RDS")
PROBS_GENZ  <- readRDS("RDS/Comparison_in_0/fung_GENZ_log2P_256.RDS")
PROBS_RIDGE <- readRDS("RDS/Comparison_in_0/fung_RIDGE_log2P_256.RDS")

colnames(PROBS_BOTEV) <- colnames(PROBS_GENZ) <- colnames(PROBS_RIDGE) <- NSAMPLE

a1 = PROBS_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = PROBS_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = PROBS_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_fung <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Fung")
ep_fung <- readRDS("RDS/Comparison_in_0/fung_EPCHOL2_log2P_256.RDS")


# -------------------------------------------------------------------------

NSAMPLE <- paste("N =",c(500,1000,5000,10000, 25000))
PROBS_BOTEV <- readRDS("RDS/Comparison_in_0/Dense_BOTEV_log2P_256.RDS")
PROBS_GENZ  <- readRDS("RDS/Comparison_in_0/Dense_GENZ_log2P_256.RDS")
PROBS_RIDGE <- readRDS("RDS/Comparison_in_0/Dense_RIDGE_log2P_256.RDS")

colnames(PROBS_BOTEV) <- colnames(PROBS_GENZ) <- colnames(PROBS_RIDGE) <- NSAMPLE

a1 = PROBS_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = PROBS_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = PROBS_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_dense <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Dense")
ep_dense <- readRDS("RDS/Comparison_in_0/Dense_EPCHOL2_log2P_256.RDS")


EP <- data.frame(ep = ep_dense, cov = "Dense") %>% 
  bind_rows(data.frame(ep= ep_fung, cov= "Fung") )%>% 
  bind_rows(data.frame(ep= ep_rho, cov= "Const."))
# -------------------------------------------------------------------------

A <- rbind(A_rho,A_fung, A_dense)
colnames(A)

ggplot()+
  geom_boxplot(data=A,aes(x=Var2,y=value,col=algo))+
  facet_wrap(~cov, scale="free_y")+
  theme_bw()+
  geom_hline(data=EP, aes(yintercept = ep))






library(tidyverse)
# extract_results ---------------------------------------------------------



NSAMPLE <- paste("N =",c(500,1000,5000,10000, 25000))
TIMES_BOTEV <- readRDS("RDS/Comparison_in_0/rho05_BOTEV_times_256.RDS")
TIMES_GENZ  <- readRDS("RDS/Comparison_in_0/rho05_TLRNK_times_256.RDS")
TIMES_RIDGE <- readRDS("RDS/Comparison_in_0/rho05_RIDGE_times_256.RDS")

colnames(TIMES_BOTEV) <- colnames(TIMES_GENZ) <- colnames(TIMES_RIDGE) <- NSAMPLE

a1 = TIMES_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = TIMES_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = TIMES_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_rho <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Const.")
ep_rho <- readRDS("RDS/Comparison_in_0/rho05_EPCHOL2_times_256.RDS")


# -------------------------------------------------------------------------


NSAMPLE <- paste("N =",c(500,1000,5000,10000, 25000))
TIMES_BOTEV <- readRDS("RDS/Comparison_in_0/fung_BOTEV_times_256.RDS")
TIMES_GENZ  <- readRDS("RDS/Comparison_in_0/fung_GENZ_times_256.RDS")
TIMES_RIDGE <- readRDS("RDS/Comparison_in_0/fung_RIDGE_times_256.RDS")

colnames(TIMES_BOTEV) <- colnames(TIMES_GENZ) <- colnames(TIMES_RIDGE) <- NSAMPLE

a1 = TIMES_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = TIMES_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = TIMES_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_fung <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Fung")
ep_fung <- readRDS("RDS/Comparison_in_0/fung_EPCHOL2_times_256.RDS")


# -------------------------------------------------------------------------

NSAMPLE <- paste("N =",c(500,1000,5000,10000, 25000))
TIMES_BOTEV <- readRDS("RDS/Comparison_in_0/Dense_BOTEV_times_256.RDS")
TIMES_GENZ  <- readRDS("RDS/Comparison_in_0/Dense_GENZ_times_256.RDS")
TIMES_RIDGE <- readRDS("RDS/Comparison_in_0/Dense_RIDGE_times_256.RDS")

colnames(TIMES_BOTEV) <- colnames(TIMES_GENZ) <- colnames(TIMES_RIDGE) <- NSAMPLE

a1 = TIMES_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = TIMES_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = TIMES_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_dense <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Dense")
ep_dense <- readRDS("RDS/Comparison_in_0/Dense_EPCHOL2_times_256.RDS")


EP <- data.frame(ep = ep_dense, cov = "Dense") %>% 
  bind_rows(data.frame(ep= ep_fung, cov= "Fung") )%>% 
  bind_rows(data.frame(ep= ep_rho, cov= "Const."))
# -------------------------------------------------------------------------

A <- rbind(A_rho,A_fung, A_dense)
colnames(A)

ggplot()+
  geom_boxplot(data=A,aes(x=Var2,y=value,col=algo))+
  facet_wrap(~cov, scale="free_y")+scale_y_log10()+
  theme_bw()+
  geom_hline(data=EP, aes(yintercept = ep))
