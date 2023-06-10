##########################################################

# IOI stability models:
#H1: testing hybrid instability

##########################################################

#R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"

rm(list = ls())

library(rstan)
library(arm)
library(MuMIn)
library(dplyr)
library(broom)
library(coda)
library(grid) 
library(bayesplot)
library(brms)
library(merTools)
library(tidybayes)
library(parallel)
library(ggplot2)
library(patchwork)
library(ggridges)
library(afex)
library(rstan)
library(openair)
library(sjPlot)
require(effects)
require(sjstats)
require(cmdstanr)
require(bayestestR)
require(readxl)
require(xlsx)
require(qdapTools)
require(readxl)
require(posterior)

#### Import data ####
d<-as.data.frame(read_excel("../Tinkerbird_stability/MS_StabilityMaster_Tinkerbirds_SingleRecs_10Mar23.xlsx")) %>% 
  drop_na(ringno)

d <- subset(d, d$ringno != "AS28707") #remove juvenile individual
str(d)

d$bird_ID<-as.factor(d$bird_ID)
d$song_id<-as.factor(d$song_id)
d$ringno<-as.factor(d$ringno)
d$prop_RF_elai<-as.numeric(d$prop_RF_elai)
d$prop_YF_elai<-as.numeric(d$prop_YF_elai)

# Scale variables
d$IOISta<-scale(d$IOI_spec1, scale = T)
d$prop_RF_elai_sta<-scale(d$prop_RF_elai, scale = T)
d$IOI_nrxn1_pp_sta<-scale(d$IOI_nrxn1_pp, scale = T)
d$IOI_coq8a_pp_sta<-scale(d$IOI_coq8a_pp, scale = T)
d$IOI_msh2_pp_sta<-scale(d$IOI_msh2_pp, scale = T)
d$IOI_enah_pp_sta<-scale(d$IOI_enah_pp, scale = T)
d$IOI_coq8a_pp_fastr_sta<-scale(d$IOI_coq8a_pp_fastr, scale = T)
d$IOI_nrxn1_pp_fastr_sta<-scale(d$IOI_nrxn1_pp_fastr, scale = T)
d$IOI_msh2_pp_fastr_sta<-scale(d$IOI_msh2_pp_fastr, scale = T)
d$IOI_enah_pp_fastr_sta<-scale(d$IOI_enah_pp_fastr, scale = T)

#keep only contact zone individuals
sympatric<-subset(d, d$Location == "Dinedor" |
                    d$Location == "Dombeya" |
                    d$Location == "Malindza" |
                    d$Location == "Mpaka" |
                    d$Location == "Mpofu" |
                    d$Location == "Skhuphe" |
                    d$Location == "Tshaneni" |
                    d$Location == "Volinde")

#### Let's model - testing hybrid instability in the contact zone - All individuals ####
DHGLMTinkerbird_cz_all=bf(IOISta ~ 1 + prop_RF_elai_sta + (1|ringno), 
                    sigma ~ I(prop_RF_elai_sta^2) + prop_RF_elai_sta + (1|ringno)) 

m1_dhglm_symp_hinst_all<- brm(DHGLMTinkerbird_cz_all, data = sympatric, 
                        warmup = 500,
                        iter = 5500,
                        chains = 5, 
                        control = list(max_treedepth = 15),
                        cores = 5,
                        seed = 12345)

m1_dhglm_symp_hinst_all <- add_criterion(m1_dhglm_symp_hinst_all, "waic")

waic_m3_all <- waic(m1_dhglm_symp_hinst_all)

plot(m1_dhglm_symp_hinst_all)

summary(m1_dhglm_symp_hinst_all)

plot_model(m1_dhglm_symp_hinst_all) +
  theme_classic()

# Remove Linear Effect
DHGLMTinkerbird_cz_all1=bf(IOISta ~ 1 + prop_RF_elai_sta + (1|ringno), 
                          sigma ~ I(prop_RF_elai_sta^2) + (1|ringno)) 

m1_dhglm_symp_hinst_all1<- brm(DHGLMTinkerbird_cz_all1, data = sympatric, 
                              warmup = 500,
                              iter = 5500,
                              chains = 5, 
                              control = list(max_treedepth = 15),
                              cores = 5,
                              seed = 12345)

m1_dhglm_symp_hinst_all1 <- add_criterion(m1_dhglm_symp_hinst_all1, "waic")

waic_m3_all1 <- waic(m1_dhglm_symp_hinst_all1)

plot(m1_dhglm_symp_hinst_all1)

summary(m1_dhglm_symp_hinst_all1)

plot_model(m1_dhglm_symp_hinst_all1) +
  theme_classic()

# Remove Quadratic Effect
DHGLMTinkerbird_cz_all2=bf(IOISta ~ 1 + prop_RF_elai_sta + (1|ringno),
                          sigma ~ prop_RF_elai_sta + (1|ringno)) 

m1_dhglm_symp_hinst_all2<- brm(DHGLMTinkerbird_cz_all2, data = sympatric, 
                              warmup = 500,
                              iter = 5500,
                              chains = 5, 
                              control = list(max_treedepth = 15),
                              cores = 5,
                              seed = 12345)

m1_dhglm_symp_hinst_all2 <- add_criterion(m1_dhglm_symp_hinst_all2, "waic")

plot(m1_dhglm_symp_hinst_all2)
summary(m1_dhglm_symp_hinst_all2)
plot_model(m1_dhglm_symp_hinst_all2) +
  theme_classic()

#### Let's model - testing hybrid instability in the contact zone - All individuals - pp NRXN1 ####
DHGLMTinkerbird_cz_all_nrxn1=bf(IOISta ~ 1 + IOI_nrxn1_pp_sta + (1|ringno),
                                sigma ~ I(IOI_nrxn1_pp_sta^2) + IOI_nrxn1_pp_sta + (1|ringno)) 

m1_dhglm_symp_hinst_all_nrxn1<- brm(DHGLMTinkerbird_cz_all_nrxn1, data = sympatric,
                                    warmup = 500,
                                    iter = 5500,
                                    chains = 5, 
                                    control = list(max_treedepth = 15),
                                    cores = 5,
                                    seed = 12345)

m1_dhglm_symp_hinst_all_nrxn1 <- add_criterion(m1_dhglm_symp_hinst_all_nrxn1, "waic")

plot(m1_dhglm_symp_hinst_all_nrxn1)
summary(m1_dhglm_symp_hinst_all_nrxn1)
plot_model(m1_dhglm_symp_hinst_all_nrxn1) +
  theme_classic()

#posterior distribution of model estimates
nrxn1_elai_draws <- brms::as_draws(m1_dhglm_symp_hinst_all_nrxn1)

sigma_intercept <- cbind(nrxn1_elai_draws$`1`$b_sigma_Intercept, 
                         nrxn1_elai_draws$`2`$b_sigma_Intercept, 
                         nrxn1_elai_draws$`3`$b_sigma_Intercept,
                         nrxn1_elai_draws$`4`$b_sigma_Intercept,
                         nrxn1_elai_draws$`5`$b_sigma_Intercept)

sigma_linear <- cbind(nrxn1_elai_draws$`1`$b_sigma_IOI_nrxn1_pp_sta, 
                         nrxn1_elai_draws$`2`$b_sigma_IOI_nrxn1_pp_sta, 
                         nrxn1_elai_draws$`3`$b_sigma_IOI_nrxn1_pp_sta,
                         nrxn1_elai_draws$`4`$b_sigma_IOI_nrxn1_pp_sta,
                         nrxn1_elai_draws$`5`$b_sigma_IOI_nrxn1_pp_sta)

sigma_quadratic <- cbind(nrxn1_elai_draws$`1`$b_sigma_IIOI_nrxn1_pp_staE2, 
                      nrxn1_elai_draws$`2`$b_sigma_IIOI_nrxn1_pp_staE2, 
                      nrxn1_elai_draws$`3`$b_sigma_IIOI_nrxn1_pp_staE2,
                      nrxn1_elai_draws$`4`$b_sigma_IIOI_nrxn1_pp_staE2,
                      nrxn1_elai_draws$`5`$b_sigma_IIOI_nrxn1_pp_staE2)

# Calculate the traits' value at peak performance and store in table
ypeak_nrxn1<-sigma_intercept-((sigma_linear^2/(4*sigma_quadratic)))

coef_ypeak_nrxn1=quantile (ypeak_nrxn1, prob=c(0.025, 0.5, 0.975))

# IOI valaue at peak
#2.5%        50%      97.5% 
#-1.2910088 -1.0549380 -0.8105339  # Backscaled is IOI = 0.47 (0.45 - 0.48)

coefYpeakTable = as.data.frame(cbind(coefName="Peak y", as.data.frame(t(as.data.frame(coef_ypeak_nrxn1)))))

# Calculate the ancestry at peak performance and store in table
xpeak_nrxn1= -(sigma_linear)/(2*sigma_quadratic)

coef_xpeak_mod2=quantile (xpeak_nrxn1, prob=c(0.025, 0.5, 0.975))

coefXpeakTable = as.data.frame(cbind(coefName="Peak x", as.data.frame(t(as.data.frame(coef_xpeak_mod2)))))

coefXpeakTable

#coefName      2.5%       50%      97.5%
#coef_xpeak_mod2   Peak x -1.138128 -0.499821 -0.2493761

## Post hoc test
# Assign a value of 0 (or no) to post-peak observations and of 1 (yes) to pre-peak observations
sympatric$PrePeak_nrxn1[sympatric$IOI_nrxn1_pp_sta< -1.05]<-"yes"

sympatric$PrePeak_nrxn1[sympatric$IOI_nrxn1_pp_sta> -1.05]<-"no"

nrxn1_posthoc <- lmer(IOISta ~ 1 + IOI_nrxn1_pp_fastr_sta * PrePeak_nrxn1 + (1|ringno),
                      data = sympatric)

summary(nrxn1_posthoc)
plot(simulateResiduals(nrxn1_posthoc))
check_model(nrxn1_posthoc)

#Simulate posterior distribution of model estimates
n.sim <- 5000 

bsim_nrxn1 <- sim(nrxn1_posthoc, n.sim=n.sim) 

colnames(bsim_nrxn1@fixef) <- names(fixef(nrxn1_posthoc))

head(bsim_nrxn1@fixef)

# Retrieve mean and credible intervals for fixed effects (effect of ancestry = post-peak effect of ancestry)
coef_fixed_nrxn1_posthoc = apply(bsim_nrxn1@fixef, 2, quantile, prob=c(0.025,0.5, 0.975)) 

head(coef_fixed_nrxn1_posthoc)

# Recalculate to retrieve pre-peak age effect from the same model
pre_peak_anc_nrxn1_posthoc= (bsim_nrxn1@fixef[,2])+(bsim_nrxn1@fixef[,4]) 

head(pre_peak_anc_nrxn1_posthoc)

coef_pre_peak_anc_nrxn1_posthoc=quantile (pre_peak_anc_nrxn1_posthoc, prob=c(0.025, 0.5, 0.975))

head(coef_pre_peak_anc_nrxn1_posthoc)

# Store pre- and post-peak age effects in table
coefPrePeakAncTable_nrxn1 = as.data.frame(cbind(coefName="Pre-peak ancestry", as.data.frame(t(as.data.frame(coef_pre_peak_anc_nrxn1_posthoc)))))

coefPrePeakAncTable_nrxn1

coefPostPeakAncTable_nrxn1 = as.data.frame(t(as.data.frame(coef_fixed_nrxn1_posthoc))) %>% tibble::rownames_to_column("coefName") %>%  filter(coefName=="IOI_nrxn1_pp_fastr_sta")

coefPostPeakAncTable_nrxn1

#### Let's model - testing hybrid instability in the contact zone - All individuals - pp COQ8A ####
DHGLMTinkerbird_cz_all_coq8a=bf(IOISta ~ 1 + IOI_coq8a_pp_sta + (1|ringno),
                                sigma ~ I(IOI_coq8a_pp_sta^2) + IOI_coq8a_pp_sta + (1|ringno)) 

m1_dhglm_symp_hinst_all_coq8a<- brm(DHGLMTinkerbird_cz_all_coq8a, data = sympatric,
                                    warmup = 500,
                                    iter = 5500,
                                    chains = 5, 
                                    control = list(max_treedepth = 15),
                                    cores = 5,
                                    seed = 56789)

m1_dhglm_symp_hinst_all_coq8a <- add_criterion(m1_dhglm_symp_hinst_all_coq8a, "waic")

plot(m1_dhglm_symp_hinst_all_coq8a)

summary(m1_dhglm_symp_hinst_all_coq8a)

plot_model(m1_dhglm_symp_hinst_all_coq8a) +
  theme_classic()

#posterior distribution of model estimates
coq8a_elai_draws <- as_draws(m1_dhglm_symp_hinst_all_coq8a)

coq8a_sigma_intercept <- cbind(coq8a_elai_draws$`1`$b_sigma_Intercept, 
                         coq8a_elai_draws$`2`$b_sigma_Intercept, 
                         coq8a_elai_draws$`3`$b_sigma_Intercept,
                         coq8a_elai_draws$`4`$b_sigma_Intercept,
                         coq8a_elai_draws$`5`$b_sigma_Intercept)

coq8a_sigma_linear <- cbind(coq8a_elai_draws$`1`$b_sigma_IOI_coq8a_pp_sta, 
                      coq8a_elai_draws$`2`$b_sigma_IOI_coq8a_pp_sta, 
                      coq8a_elai_draws$`3`$b_sigma_IOI_coq8a_pp_sta,
                      coq8a_elai_draws$`4`$b_sigma_IOI_coq8a_pp_sta,
                      coq8a_elai_draws$`5`$b_sigma_IOI_coq8a_pp_sta)

coq8a_sigma_quadratic <- cbind(coq8a_elai_draws$`1`$b_sigma_IIOI_coq8a_pp_staE2, 
                         coq8a_elai_draws$`2`$b_sigma_IIOI_coq8a_pp_staE2, 
                         coq8a_elai_draws$`3`$b_sigma_IIOI_coq8a_pp_staE2,
                         coq8a_elai_draws$`4`$b_sigma_IIOI_coq8a_pp_staE2,
                         coq8a_elai_draws$`5`$b_sigma_IIOI_coq8a_pp_staE2)

# Calculate the traits' value at peak performance and store in table
ypeak_coq8a<-coq8a_sigma_intercept-((coq8a_sigma_linear^2/(4*coq8a_sigma_quadratic)))

coef_ypeak_coq8a=quantile (ypeak_coq8a, prob=c(0.025, 0.5, 0.975))

coef_ypeak_coq8a

#2.5%       50%     97.5% 
#-1.320964 -1.059615 -0.812429   # Backscaled is IOI = 0.47 (0.45 - 0.48)

coefYpeakTable = as.data.frame(cbind(coefName="Peak y", as.data.frame(t(as.data.frame(coef_ypeak_coq8a)))))

coefYpeakTable

# Calculate the ancestry at peak performance and store in table
xpeak_coq8a= -(coq8a_sigma_linear)/(2*coq8a_sigma_quadratic)

coef_xpeak_coq8a=quantile (xpeak_coq8a, prob=c(0.025, 0.5, 0.975))

coefXpeakTable = as.data.frame(cbind(coefName="Peak x", as.data.frame(t(as.data.frame(coef_xpeak_coq8a)))))

coefXpeakTable

#coefName      2.5%       50%      97.5%
#coef_xpeak_mod2   Peak x -1.138128 -0.499821 -0.2493761

## Post hoc test
# Assign a value of 0 (or no) to post-peak observations and of 1 (yes) to pre-peak observations
sympatric$PrePeak_coq8a[sympatric$IOI_coq8a_pp_sta< -1.05]<-"yes"

sympatric$PrePeak_coq8a[sympatric$IOI_coq8a_pp_sta> -1.05]<-"no"

ympatric$PrePeak_coq8a[sympatric$IOI_coq8a_pp_sta< -1.05]<-"yes"
sympatric$PrePeak_coq8a[sympatric$IOI_coq8a_pp_sta> -1.05]<-"no"


coq8a_posthoc <- lmer(IOISta ~ 1 + IOI_coq8a_pp_fastr_sta * PrePeak_coq8a + (1|ringno),
                      data = sympatric)

summary(coq8a_posthoc)

plot(simulateResiduals(coq8a_posthoc))

check_model(coq8a_posthoc)

plot(allEffects(coq8a_posthoc))

#Simulate posterior distribution of model estimates
n.sim <- 5000 

bsim_coq8a <- sim(coq8a_posthoc, n.sim=n.sim) 

colnames(bsim_coq8a@fixef) <- names(fixef(coq8a_posthoc))

head(bsim_coq8a@fixef)

# Retrieve mean and credible intervals for fixed effects (effect of ancestry = post-peak effect of ancestry)
coef_fixed_coq8a_posthoc = apply(bsim_coq8a@fixef, 2, quantile, prob=c(0.025,0.5, 0.975)) 

head(coef_fixed_coq8a_posthoc)

# Recalculate to retrieve pre-peak ancestry effect from the same model
pre_peak_anc_coq8a_posthoc= (bsim_coq8a@fixef[,2])+(bsim_coq8a@fixef[,4]) 

head(pre_peak_anc_coq8a_posthoc)

coef_pre_peak_anc_coq8a_posthoc=quantile (pre_peak_anc_coq8a_posthoc, prob=c(0.025, 0.5, 0.975))

head(coef_pre_peak_anc_coq8a_posthoc)

# Store pre- and post-peak ancestry effects in table
coefPrePeakAncTable = as.data.frame(cbind(coefName="Pre-peak ancestry", as.data.frame(t(as.data.frame(coef_pre_peak_anc_coq8a_posthoc)))))

coefPrePeakAncTable

coefPostPeakAncTable = as.data.frame(t(as.data.frame(coef_fixed_coq8a_posthoc))) %>% tibble::rownames_to_column("coefName") %>%  filter(coefName=="IOI_coq8a_pp_fastr_sta")

coefPostPeakAncTable

#### Let's model - testing hybrid instability in the contact zone - All individuals - pp MSH2 ####
DHGLMTinkerbird_cz_all_msh2=bf(IOISta ~ 1 + IOI_msh2_pp_sta + (1|ringno),
                                sigma ~ I(IOI_msh2_pp_sta^2) + IOI_msh2_pp_sta + (1|ringno)) 

m1_dhglm_symp_hinst_all_msh2<- brm(DHGLMTinkerbird_cz_all_msh2, data = sympatric,
                                    warmup = 500,
                                    iter = 5500,
                                    chains = 5, 
                                    control = list(max_treedepth = 15),
                                    cores = 5,
                                    seed = 12345)

m1_dhglm_symp_hinst_all_msh2 <- add_criterion(m1_dhglm_symp_hinst_all_msh2, "waic")

plot(m1_dhglm_symp_hinst_all_msh2)

summary(m1_dhglm_symp_hinst_all_msh2)

plot_model(m1_dhglm_symp_hinst_all_msh2) +
  theme_classic()

#### Let's model - testing hybrid instability in the contact zone - All individuals - pp ENAH ####
DHGLMTinkerbird_cz_all_enah=bf(IOISta ~ 1 + IOI_enah_pp_sta + (1|ringno),
                               sigma ~ I(IOI_enah_pp_sta^2) + IOI_enah_pp_sta + (1|ringno)) 

m1_dhglm_symp_hinst_all_msh2<- brm(DHGLMTinkerbird_cz_all_enah, data = sympatric,
                                   warmup = 500,
                                   iter = 5500,
                                   chains = 5, 
                                   control = list(max_treedepth = 15),
                                   cores = 5,
                                   seed = 12345)

m1_dhglm_symp_hinst_all_msh2 <- add_criterion(m1_dhglm_symp_hinst_all_msh2, "waic")

plot(m1_dhglm_symp_hinst_all_msh2)

summary(m1_dhglm_symp_hinst_all_msh2)

plot_model(m1_dhglm_symp_hinst_all_msh2) +
  theme_classic()