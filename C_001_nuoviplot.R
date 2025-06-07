library(tidyverse)
theme_set(theme_bw())
dense <- readRDS("RDS/ALL_dense_XtX_2025.RDS")

table(dense$algo)
eca <- dense %>% filter(algo == "EP_Chol_algo")
ecb <- dense %>% filter(algo == "EP_Chol_algoB")
eea <- dense %>% filter(algo == "EP_Eigen_algoA")
eeb <- dense %>% filter(algo == "EP_Eigen_algoB")

pairs(cbind(eca$log2prob, ecb$log2prob, eea$log2prob, eeb$log2prob))


pairs(cbind((eca$log2prob-eca$log2prob), (ecb$log2prob-eca$log2prob)/eca$log2prob, 
                          (eea$log2prob-eca$log2prob)/eca$log2prob, 
                          (eeb$log2prob-eca$log2prob)/eca$log2prob ))

boxplot(cbind((eca$log2prob-eca$log2prob), (ecb$log2prob-eca$log2prob)/eca$log2prob, 
            (eea$log2prob-eca$log2prob)/eca$log2prob, 
            (eeb$log2prob-eca$log2prob)/eca$log2prob ))



dense <- readRDS("RDS/ALL_dense_XtX_2025.RDS")
dense_epcB <- dense %>% filter(algo != "EP_Eigen_algoA" , 
                 algo != "EP_Eigen_algoB" , 
                 algo != "EP_Chol_algo")


table(dense_epcB$algo)
long_epcb <- dense_epcB %>% filter(algo == "EP_Chol_algoB") %>% pull(log2prob) %>% rep(4)
long_epcb2 <- dense_epcB %>% filter(algo == "Botev") %>% pull(info) %>% as.numeric %>% rep(4)

d1 <- dense_epcB %>% mutate(log2ep = long_epcb[1:4600]) %>% 
  mutate(log2relative = (log2prob-log2ep)/log2ep, log2delta = (log2prob-log2ep)) %>% 
  filter(algo != "EP_Chol_algoB")

d1$log2prob-d1$log2ep
#View(d1)
ggplot(d1)+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = log2relative),outlier.shape = NA)+
  facet_grid(algo~dim)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(d1)+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = abs(log2delta)))+
  facet_grid(algo~dim)

d1 <- dense_epcB %>% mutate(log2ep = long_epcb[1:4600], err = long_epcb2[1:4600]) %>%
  mutate(log2relative = (log2prob-log2ep)/log2ep, log2delta = (log2prob-log2ep),
         lowerr= (log2(2^log2prob*(1-err))-log2ep)/log2ep,
         upperr= (log2(2^log2prob*(1+err))-log2ep)/log2ep) %>% 
  group_by(dim, algo) %>% mutate(mlowerr=mean(lowerr, na.rm=TRUE),
                                 mupperr=mean(upperr, na.rm=TRUE)) %>% ungroup() %>% 
  filter(algo != "EP_Chol_algoB") %>% mutate(mlowerr2 = case_when(algo != "Botev"~NA,
                                                                  algo == "Botev"~mlowerr),
                                             mupperr2 = case_when(algo != "Botev"~NA,
                                                                  algo == "Botev"~mupperr))

ggplot(d1)+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_hline(aes(yintercept = c(mlowerr2)),col=4)+
  geom_hline(aes(yintercept = c(mupperr2)),col=4)+  
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = log2relative))+
  facet_grid(algo~dim)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# -------------------------------------------------------------------------



library(tidyverse)
theme_set(theme_bw())
fung <- readRDS("RDS/ALL_fungible_2025.RDS")

table(fung$algo)
eca <- fung %>% filter(algo == "EP_Chol_algo")
ecb <- fung %>% filter(algo == "EP_Chol_algoB")
eea <- fung %>% filter(algo == "EP_Eigen_algoA")
eeb <- fung %>% filter(algo == "EP_Eigen_algoB")

pairs(cbind(eca$log2prob, ecb$log2prob, eea$log2prob, eeb$log2prob))


pairs(cbind((eca$log2prob-eca$log2prob), (ecb$log2prob-eca$log2prob)/eca$log2prob, 
            (eea$log2prob-eca$log2prob)/eca$log2prob, 
            (eeb$log2prob-eca$log2prob)/eca$log2prob ))

boxplot(cbind((eca$log2prob-eca$log2prob), (ecb$log2prob-eca$log2prob)/eca$log2prob, 
              (eea$log2prob-eca$log2prob)/eca$log2prob, 
              (eeb$log2prob-eca$log2prob)/eca$log2prob ))



fung_epcB <- fung %>% filter(algo != "EP_Eigen_algoA" , 
                               algo != "EP_Eigen_algoB" , 
                               algo != "EP_Chol_algo")


table(fung_epcB$algo)
long_epcb <- fung_epcB %>% filter(algo == "EP_Chol_algoB") %>% pull(log2prob) %>% rep(4)
long_epcb2 <- fung_epcB %>% filter(algo == "Botev") %>% pull(info) %>% as.numeric %>% rep(4)

d1 <- fung_epcB %>% mutate(log2ep = long_epcb[1:4600]) %>% 
  mutate(log2relative = (log2prob-log2ep)/log2ep, log2delta = (log2prob-log2ep)) %>% 
  filter(algo != "EP_Chol_algoB")

d1$log2prob-d1$log2ep
ggplot(d1)+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = log2relative))+
  facet_grid(algo~dim)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(d1)+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = abs(log2delta)))+
  facet_grid(algo~dim)


#mean((log2EP - log2(pBotev*(1-err)))/log2EP) e mean((log2EP - log2(pBotev*(1+err)))/log2EP)

d1 <- fung_epcB %>% mutate(log2ep = long_epcb[1:4600], err = long_epcb2[1:4600]) %>%
  mutate(log2relative = (log2prob-log2ep)/log2ep, log2delta = (log2prob-log2ep),
         lowerr= (log2(2^log2prob*(1-err))-log2ep)/log2ep,
         upperr= (log2(2^log2prob*(1+err))-log2ep)/log2ep) %>% 
  group_by(dim, algo) %>% mutate(mlowerr=mean(lowerr, na.rm=TRUE),
                                 mupperr=mean(upperr, na.rm=TRUE)) %>% ungroup() %>% 
  filter(algo != "EP_Chol_algoB") %>% mutate(mlowerr2 = case_when(algo != "Botev"~NA,
                                                                  algo == "Botev"~mlowerr),
                                             mupperr2 = case_when(algo != "Botev"~NA,
                                                                  algo == "Botev"~mupperr))

ggplot(d1)+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_hline(aes(yintercept = c(mlowerr2)),col=4)+
  geom_hline(aes(yintercept = c(mupperr2)),col=4)+  
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = log2relative))+
  facet_grid(algo~dim)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# -------------------------------------------------------------------------




library(tidyverse)
theme_set(theme_bw())
truth <- readRDS("RDS/00_trueLOGProbs_constantSigma.RDS")
rhosf <- readRDS("RDS/ALL_rhofixed_2025.RDS")

eca <- rhosf %>% filter(algo == "EP_Chol_algo")
ecb <- rhosf %>% filter(algo == "EP_Chol_algoB")
eea <- rhosf %>% filter(algo == "EP_Eigen_algoA")
eeb <- rhosf %>% filter(algo == "EP_Eigen_algoB")

pairs(cbind(eca$log2prob, ecb$log2prob, eea$log2prob, eeb$log2prob))


pairs(cbind((eca$log2prob-eca$log2prob), (ecb$log2prob-eca$log2prob)/eca$log2prob, 
            (eea$log2prob-eca$log2prob)/eca$log2prob, 
            (eeb$log2prob-eca$log2prob)/eca$log2prob ))

boxplot(cbind((eca$log2prob-eca$log2prob), (ecb$log2prob-eca$log2prob)/eca$log2prob, 
              (eea$log2prob-eca$log2prob)/eca$log2prob, 
              (eeb$log2prob-eca$log2prob)/eca$log2prob ))



rhosf_epcB <- rhosf %>% filter(algo != "EP_Eigen_algoA" , 
                               algo != "EP_Eigen_algoB" , 
                               algo != "EP_Chol_algo")


table(rhosf_epcB$algo)

rhosf_epcB <- rhosf_epcB %>% left_join(truth %>% rename(dim = Ns),by = c("b","dim","rho")) %>% mutate(log2true = Truep / log(2))
rhosf_epcb2 <- rhosf_epcB %>% filter(algo == "Botev") %>% pull(info) %>% as.numeric %>% rep(4)
d1 <- rhosf_epcB %>% 
  mutate(log2relative = (2^log2prob-2^log2true)/2^log2true, log2delta = (2^log2prob-2^log2true))
View(d1)

ggplot(d1 %>% filter(rho == 0))+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = log2relative))+
  facet_grid(algo~dim)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("rho = 0")


ggplot(d1 %>% filter(rho == 0.25))+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = log2relative))+
  facet_grid(algo~dim)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("rho = .25")

ggplot(d1 %>% filter(rho == 0.5))+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = log2relative))+
  facet_grid(algo~dim,scale = "free_y")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("rho = .5")

ggplot(d1 %>% filter(rho == 0.75))+
  geom_hline(aes(yintercept = 0),col=2)+
  geom_boxplot(aes(x = factor(round(b,2)),
                   y = log2relative))+
  facet_grid(algo~dim)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("rho = .75")



