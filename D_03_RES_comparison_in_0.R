library(tidyverse)
library(reshape2)
library(latex2exp)
# extract_results ---------------------------------------------------------



NSAMPLE <- paste(c(5000,10000,20000,50000))
PROBS_BOTEV <- readRDS("RDS/Comparison_in_0/rho05_BOTEV_log2P_256.RDS")
PROBS_GENZ  <- readRDS("RDS/Comparison_in_0/rho05_TLRNK_log2P_256.RDS")
PROBS_RIDGE <- readRDS("RDS/Comparison_in_0/rho05_RIDGE_log2P_256.RDS")

colnames(PROBS_BOTEV) <- colnames(PROBS_GENZ) <- colnames(PROBS_RIDGE) <- NSAMPLE

a1 = PROBS_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = PROBS_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = PROBS_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_rho <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = factor("Case~(ii)*','~rho==0.5") )
ep_rho <- readRDS("RDS/Comparison_in_0/rho05_EPCHOL2_log2P_256.RDS")

# -------------------------------------------------------------------------

PROBS_BOTEV <- readRDS("RDS/Comparison_in_0/fung_BOTEV_log2P_256.RDS")
PROBS_GENZ  <- readRDS("RDS/Comparison_in_0/fung_GENZ_log2P_256.RDS")
PROBS_RIDGE <- readRDS("RDS/Comparison_in_0/fung_RIDGE_log2P_256.RDS")

colnames(PROBS_BOTEV) <- colnames(PROBS_GENZ) <- colnames(PROBS_RIDGE) <- NSAMPLE

a1 = PROBS_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = PROBS_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = PROBS_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_fung <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Case (i)")
ep_fung <- readRDS("RDS/Comparison_in_0/fung_EPCHOL2_log2P_256.RDS")


# -------------------------------------------------------------------------

PROBS_BOTEV <- readRDS("RDS/Comparison_in_0/Dense_BOTEV_log2P_256.RDS")
PROBS_GENZ  <- readRDS("RDS/Comparison_in_0/Dense_GENZ_log2P_256.RDS")
PROBS_RIDGE <- readRDS("RDS/Comparison_in_0/Dense_RIDGE_log2P_256.RDS")

colnames(PROBS_BOTEV) <- colnames(PROBS_GENZ) <- colnames(PROBS_RIDGE) <- NSAMPLE

a1 = PROBS_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = PROBS_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = PROBS_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_dense <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Case (iii)")
ep_dense <- readRDS("RDS/Comparison_in_0/Dense_EPCHOL2_log2P_256.RDS")


EP <- data.frame(ep = ep_dense, cov = "Case (iii)") %>% 
  bind_rows(data.frame(ep= ep_fung, cov= "Case (i)") )%>% 
  bind_rows(data.frame(ep= ep_rho, cov= factor("Case~(ii)*','~rho==0.5")))
# -------------------------------------------------------------------------

EP$cov <- factor(EP$cov, levels = c("Case (i)", "Case~(ii)*','~rho==0.5", "Case (iii)"))

A <- rbind(A_rho,A_fung, A_dense)
colnames(A)

ggplot()+
  geom_hline(data=EP, aes(yintercept = ep), col= 2, lwd=1.5, lty=1)+
  geom_boxplot(data=A,aes(x=factor(Var2),y=value,col=algo))+
  scale_color_manual(values = c(1,4,3))+
  facet_wrap(~cov, scale="free_y", label = label_parsed)+
  theme_bw()+
  xlab("N")+ylab(TeX("Estimated $log_2$ probabilities"))+
  theme(legend.position = "bottom", text = element_text(size=18))+
  ggview::canvas(h=5,w=15)


ggplot()+
  geom_hline(data=EP, aes(yintercept = ep), col= 2, lwd=1.5, lty=1)+
  geom_boxplot(data=A,aes(x=factor(Var2),y=value,col=algo))+
  scale_color_manual(values = c(1,4,3))+
  facet_wrap(~cov, scale="free_y", label = label_parsed)+
  theme_bw()+
  xlab("N")+ylab(TeX("Estimated $log_2$ probabilities"))+
  theme(legend.position = "none", text = element_text(size=18))+
  ggview::canvas(h=5,w=15)
ggsave("NewFigures/log2p_comparisonin0_m256.pdf",h=6,w=18)




library(tidyverse)
# extract_results ---------------------------------------------------------

TIMES_BOTEV <- readRDS("RDS/Comparison_in_0/rho05_BOTEV_times_256.RDS")
TIMES_GENZ  <- readRDS("RDS/Comparison_in_0/rho05_TLRNK_times_256.RDS")
TIMES_RIDGE <- readRDS("RDS/Comparison_in_0/rho05_RIDGE_times_256.RDS")

colnames(TIMES_BOTEV) <- colnames(TIMES_GENZ) <- colnames(TIMES_RIDGE) <- NSAMPLE

a1 = TIMES_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = TIMES_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = TIMES_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_rho <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = factor("Case~(ii)*','~rho==0.5"))
ep_rho <- readRDS("RDS/Comparison_in_0/rho05_EPCHOL2_times_256.RDS")


# -------------------------------------------------------------------------

TIMES_BOTEV <- readRDS("RDS/Comparison_in_0/fung_BOTEV_times_256.RDS")
TIMES_GENZ  <- readRDS("RDS/Comparison_in_0/fung_GENZ_times_256.RDS")
TIMES_RIDGE <- readRDS("RDS/Comparison_in_0/fung_RIDGE_times_256.RDS")

colnames(TIMES_BOTEV) <- colnames(TIMES_GENZ) <- colnames(TIMES_RIDGE) <- NSAMPLE

a1 = TIMES_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = TIMES_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = TIMES_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_fung <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Case (i)")
ep_fung <- readRDS("RDS/Comparison_in_0/fung_EPCHOL2_times_256.RDS")


# -------------------------------------------------------------------------

TIMES_BOTEV <- readRDS("RDS/Comparison_in_0/Dense_BOTEV_times_256.RDS")
TIMES_GENZ  <- readRDS("RDS/Comparison_in_0/Dense_GENZ_times_256.RDS")
TIMES_RIDGE <- readRDS("RDS/Comparison_in_0/Dense_RIDGE_times_256.RDS")

colnames(TIMES_BOTEV) <- colnames(TIMES_GENZ) <- colnames(TIMES_RIDGE) <- NSAMPLE

a1 = TIMES_BOTEV %>% reshape2::melt() %>% mutate(algo = "Botev")
a2 = TIMES_GENZ %>% reshape2::melt() %>% mutate(algo = "Genz/TLRank")
a3 = TIMES_RIDGE %>% reshape2::melt() %>% mutate(algo = "Ridgeway")

A_dense <- a1 %>% bind_rows(a2,a3) %>% mutate(cov = "Case (iii)")
ep_dense <- readRDS("RDS/Comparison_in_0/Dense_EPCHOL2_times_256.RDS")


EP <- data.frame(ep = ep_dense, cov = "Case (iii)") %>% 
  bind_rows(data.frame(ep= ep_fung, cov= "Case (i)") )%>% 
  bind_rows(data.frame(ep= ep_rho, cov= factor("Case~(ii)*','~rho==0.5")))
# -------------------------------------------------------------------------
EP$cov <- factor(EP$cov, levels = c("Case (i)", "Case~(ii)*','~rho==0.5", "Case (iii)"))

A <- rbind(A_rho,A_fung, A_dense)
colnames(A)

ggplot()+
  geom_hline(data=EP, aes(yintercept = ep), col= 2, lwd=1.5, lty=1)+
  geom_boxplot(data=A,aes(x=factor(Var2),y=value,col=algo))+
  scale_color_manual(values = c(1,4,3))+
  facet_wrap(~cov, scale="free_y", label = label_parsed)+
  theme_bw()+scale_y_log10()+
  xlab("N")+ylab(TeX("Estimated $log_2$ probabilities"))+
  theme(legend.position = "bottom", text = element_text(size=18))+
  ggview::canvas(h=5,w=15)


ggplot()+
  geom_hline(data=EP, aes(yintercept = ep), col= 2, lwd=1.5, lty=1)+
  geom_boxplot(data=A,aes(x=factor(Var2),y=value,col=algo))+
  scale_color_manual(values = c(1,4,3))+
  facet_wrap(~cov, scale="free_y", label = label_parsed)+
  theme_bw()+scale_y_log10()+
  xlab("N")+ylab(TeX("Computational time in seconds ($log_{10} scale$)"))+
  theme(legend.position = "none", text = element_text(size=18))+
  ggview::canvas(h=6,w=18)

ggsave("NewFigures/times_comparisonin0_m256.pdf",h=6,w=18)
