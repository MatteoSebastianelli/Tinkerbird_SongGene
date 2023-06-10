##########################################################

#IOI stability models:
#H2: testing character displacement

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
allopatric<-as.data.frame(read_excel("../Tinkerbird_stability/MS_StabilityMaster_Tinkerbirds_SingleRecs_10Mar23.xlsx")) %>% 
  filter(pop == "allopatric")

allopatric$bird_ID<-as.factor(allopatric$bird_ID)
allopatric$song_id<-as.factor(allopatric$song_id)
allopatric$ringno<-as.factor(allopatric$ringno)
allopatric$prop_RF_elai<-as.numeric(allopatric$prop_RF_elai)
allopatric$prop_YF_elai<-as.numeric(allopatric$prop_YF_elai)
allopatric$prop_RF_fastr<-as.numeric(allopatric$prop_RF_fastr)
allopatric$prop_YF_fastr<-as.numeric(allopatric$prop_YF_fastr)


allopatric$IOISta<-scale(allopatric$IOI_spec1, scale = T)
allopatric$prop_RF_elai_sta<-scale(allopatric$prop_RF_elai, scale = T)

# Merge ringed individuals from the contact zone with all allopatric individuals
purged <- rbind(sympatric, allopatric) # The sympatric dataset was generated in the previous script
head(purged)

# Get only >= 99 % pure individuals
pure_filtered_rf_99 <- purged %>% filter(prop_RF_fastr >= 0.99)

pure_filtered_yf_99 <- purged %>% filter(prop_YF_fastr >= 0.99)

pure_filtered_99 <- rbind(pure_filtered_rf_99, pure_filtered_yf_99)

head(pure_filtered_99)

#### Let's model - testing character displacement in the cz with pure filtered 99 ####
DHGLMTinkerbird_cd1=bf(IOISta ~ 1 + Species*pop + (1|bird_ID), 
                      sigma ~ Species*pop + (1|bird_ID))

m2_dhglm_pure99CD<- brm(DHGLMTinkerbird_cd1, data = pure_filtered_99, # prior = prior,
                      warmup = 500,
                      iter = 7000,
                      chains = 5, 
                      #thin = 8,
                      control = list(max_treedepth = 15),
                      cores = 5,
                      seed = 12345)

m2_dhglm_pure99CD <- add_criterion(m2_dhglm_pure99CD, "waic")

plot(m2_dhglm_pure99CD)
summary(m2_dhglm_pure99CD)
plot_model(m2_dhglm_pure99CD) +
  theme_classic()