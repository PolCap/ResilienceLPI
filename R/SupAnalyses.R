# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupAnalyses.R
# - DATE:        20/11/2020
# - DESCRIPTION: Supplementary analyses using the lambda trends to analyse the 
#                other factors determining the poulation change
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list = ls(all = TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(brms)
library(data.table)
library(zoo)
library(dplyr)
library(expss)
library(MASS)
library(tidybayes)
library(bayesplot)
library(bayestestR)

# Set working directory

path <- gsub("LPI/R", "LPI/",
             dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path, "R")
DataPath <-  paste0(path, "Data")
ResultPath <-  paste0(path, "Results")

# Resistance ###################################################################

# Load data

load(paste0(DataPath, "/ResistanceChange.RData"))
load(paste0(DataPath, "/RecoveryChange.RData"))

# Read csv with body sizes

setwd(DataPath)
bm <- read.csv("Body_Mass_means.csv")

# Match body mass data

bm$Binomial <- gsub("_", " ", bm$Binomial)

res_data <- left_join(res_data, bm[, c("Binomial", "bs_g")],
            by = c("SpeciesName" = "Binomial"))

rec_data <- left_join(rec_data, bm[, c("Binomial", "bs_g")],
            by = c("SpeciesName" = "Binomial"))

# Transform number of threats in factor

res_data <- res_data %>%
  mutate(n.threat = as.factor(n.threat),
         n.threat = factor(n.threat, levels = c("0", "1", "2", "3")))

rec_data <- rec_data %>%
  mutate(n.threat = as.factor(n.threat),
         n.threat = factor(n.threat, levels = c("0", "1", "2", "3")))

# Classify taxonomic groups into major ones

rec_data <- rec_data %>% 
  mutate(Taxon= ifelse(Class=="Holocephali"|Class=="Elasmobranchii" | 
                         Class=="Myxini"|Class=="Cephalaspidomorphi"|
                         Class=="Actinopterygii"|Class=="Sarcopterygii",
                       "Fish", 
                       ifelse(Class=="Aves", "Birds",
                              ifelse(Class=="Mammalia", 
                                     "Mammals",
                                     ifelse(Class=="Amphibia", 
                                            "Amphibians",
                                            ifelse(Class=="Reptilia",
                                                   "Reptiles"))))))

res_data <- res_data %>% 
  mutate(Taxon= ifelse(Class=="Holocephali"|Class=="Elasmobranchii" | 
                         Class=="Myxini"|Class=="Cephalaspidomorphi"|
                         Class=="Actinopterygii"|Class=="Sarcopterygii",
                       "Fish", 
                       ifelse(Class=="Aves", "Birds",
                              ifelse(Class=="Mammalia", 
                                     "Mammals",
                                     ifelse(Class=="Amphibia", 
                                            "Amphibians",
                                            ifelse(Class=="Reptilia",
                                                   "Reptiles"))))))

# Transform to a data frame

res_data <- as.data.frame(res_data)
rec_data <- as.data.frame(rec_data)

#Change tau for the analyses

colnames(res_data)[15] <- "tau"
colnames(rec_data)[15] <- "tau"

# Set modelling parameters

iter = 10000
thin = 0.0005 * iter
warmup = 0.1 * iter
prior <- c(prior(normal(0, 1), class = b),
           prior(exponential(1), class = sigma))

# Resistance ###################################################################

# Body size --------------------------------------------------------------------

bs_res <- brm(mu | se(tau, sigma = TRUE) ~ scale(bs_g)*Taxon + (1|SpeciesName),
              iter = iter, thin = thin, warmup = warmup,
              prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
              data = subset(res_data, !is.na(bs_g)),
              family = gaussian,
              cores = 20)

# Length of the time series ----------------------------------------------------

duration_res <- brm(mu | se(tau, sigma = TRUE) ~ scale(Duration) + (1|SpeciesName),
                    iter = iter, thin = thin,  warmup = warmup,
                    prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
                    data = res_data,
                    family = gaussian,
                    cores = 20)


# Protection status ------------------------------------------------------------

protect_res <- brm(mu | se(tau, sigma = TRUE) ~ Protected_status - 1 + (1|SpeciesName) ,
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, 
                   control = list(adapt_delta = .975, max_treedepth = 20),
                   data = res_data,
                   family = gaussian,
                   cores = 20)

# Managed ----------------------------------------------------------------------

managed_res <- brm(mu | se(tau, sigma = TRUE) ~ Managed - 1 + (1|SpeciesName) ,
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, 
                   control = list(adapt_delta = .975, max_treedepth = 20),
                   data = res_data,
                   family = gaussian,
                   cores = 20)

# System taxa ------------------------------------------------------------------

mst_res <- brm(mu | se(tau, sigma = TRUE) ~ System:Taxon - 1 + (1|SpeciesName),
               iter = iter, thin = thin, warmup = warmup,
               prior= prior, 
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = res_data,
               family = gaussian,
               cores = 20)

# Time sensitivity -------------------------------------------------------------

mts_res <- brm(mu | se(tau, sigma = TRUE) ~ Decade - 1 + (1|SpeciesName) ,
               iter = iter, thin = thin, warmup = warmup,
               prior= prior, 
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = res_data,
               family = gaussian,
               cores = 20)

# Recovery #####################################################################

# Body size --------------------------------------------------------------------

bs_rec <- brm(mu | se(tau, sigma = TRUE) ~ scale(bs_g)*Taxon + (1|SpeciesName) ,
              iter = iter, thin = thin, warmup = warmup,
              prior= prior, 
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = subset(rec_data, !is.na(bs_g)),
              family = gaussian,
              cores = 20)

# Length of the time series ----------------------------------------------------

duration_rec <- brm(mu | se(tau, sigma = TRUE) ~ scale(Duration) + (1|SpeciesName) ,
                    iter = iter, thin = thin, warmup = warmup,
                    prior= prior, 
                    control = list(adapt_delta = .975, max_treedepth = 20),
                    data = rec_data,
                    family = gaussian,
                    cores = 20)

# Protection status ------------------------------------------------------------

protect_rec <- brm(mu | se(tau, sigma = TRUE) ~ Protected_status - 1 + (1|SpeciesName) ,
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, 
                   control = list(adapt_delta = .975, max_treedepth = 20),
                   data = rec_data,
                   family = gaussian,
                   cores = 20)

# Managed  ---------------------------------------------------------------------

managed_rec<- brm(mu | se(tau, sigma = TRUE) ~ Managed - 1 + (1|SpeciesName) ,
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, 
                  control = list(adapt_delta = .975, max_treedepth = 20),
                  data = rec_data,
                  family = gaussian,
                  cores = 20)

# Time sensitivity ------------------------------------------------------------- 

mts_rec <- brm(mu | se(tau, sigma = TRUE) ~ Decade - 1 + (1|SpeciesName) ,
               iter = iter, thin = thin, warmup = warmup,
               prior= prior, 
               control = list(adapt_delta = .975, max_treedepth = 20),
               data = rec_data,
               family = gaussian,
               cores = 20)

# System taxon -----------------------------------------------------------------

mst_rec <- brm(mu | se(tau, sigma = TRUE) ~ System:Taxon - 1 + (1|SpeciesName) ,
               iter = iter, thin = thin, warmup = warmup,
               prior= prior, 
               control = list(adapt_delta = .975, max_treedepth = 20),
               data =rec_data,
               family = gaussian,
               cores = 20)

# Save the models ##############################################################

setwd(ResultPath)

save(bs_res, duration_res, protect_res, managed_res, 
     bs_rec, duration_rec, protect_rec, managed_rec,
     mts_res, mts_rec, mst_rec, mst_rec, 
     file = "SupMod.RData")

# Re-analyses without the outliers #############################################

# Function to remove outliers 

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Remove outliers from the data

rec_data$mu <- remove_outliers(rec_data$mu)
rec_data <- rec_data %>% filter(!is.na(mu))
res_data$mu <- remove_outliers(res_data$mu)
res_data <- res_data %>% filter(!is.na(mu))

# Recovery #####################################################################

# Number of threats ------------------------------------------------------------

mo1 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat-1 + (1|SpeciesName), 
           iter = iter, thin = thin, warmup = warmup,
           prior = prior,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = rec_data, family = gaussian,
           cores = 20)

# System -----------------------------------------------------------------------

mo2 <- brm(mu ~ System-1 + (1|SpeciesName), 
           iter = iter, thin = thin, warmup = warmup,
           prior = prior,
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = rec_data, family = gaussian,
           cores = 20)

# Taxa -------------------------------------------------------------------------

mo3 <- brm(mu | se(tau, sigma = TRUE) ~ Taxon-1 + (1|SpeciesName), 
           iter = iter, thin = thin, warmup = warmup,
           prior, control = list(adapt_delta = .975, max_treedepth = 20),
           data = rec_data, family = gaussian,
           cores = 20)

# Number of threats and system -------------------------------------------------

mio1 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:System - 1 + (1|SpeciesName),
            iter = iter, thin = thin,  warmup = warmup,
            prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
            data = rec_data,
            family = gaussian,
            cores = 20)

# Number of threats and Taxa ---------------------------------------------------

mio2 <-  brm(mu | se(tau, sigma = TRUE) ~ n.threat:Taxon - 1 + (1|SpeciesName),
      iter = iter, thin = thin,  warmup = warmup,
      prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
      data = rec_data,
      family = gaussian,
      cores = 20)

# Resistance ###################################################################

# Number of threats ------------------------------------------------------------

mo4 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat-1 + (1|SpeciesName), 
           iter = iter, thin = thin, warmup = warmup,
           prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
           data = res_data, family = gaussian,
           cores = 20)


# System -----------------------------------------------------------------------

mo5 <- brm(mu | se(tau, sigma = TRUE) ~ System-1 + (1|SpeciesName), 
           iter = iter, thin = thin, warmup = warmup,
           prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
           data = res_data, family = gaussian,
           cores = 20)

# Taxa -------------------------------------------------------------------------

mo6 <- brm(mu | se(tau, sigma = TRUE) ~ Taxon-1 + (1|SpeciesName), 
           iter = iter, thin = thin, warmup = warmup,
           prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
           data = res_data, family = gaussian,
           cores = 20)

# Number of threats and system -------------------------------------------------

mio3 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:System - 1 + (1|SpeciesName),
            iter = iter, thin = thin,  warmup = warmup,
            prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
            data = res_data,
            family = gaussian,
            cores = 20)

# Number of threats and taxa ---------------------------------------------------

mio4 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:Taxon - 1 + (1|SpeciesName),
            iter = iter,
            thin = thin,
            warmup = warmup,
            prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
            data = res_data,
            family = gaussian,
            cores = 20)

# Save outlier removed models ##################################################

setwd(DataPath)
save(file = "Outliers.RData", 
     mo1,mo2,mo3,mo4,mo5,mo6,
     mio1, mio2, mio3, mio4)

# Models for only declining populations ########################################

# Load population data

load(paste0(DataPath,"/PopChangeData.RData"))

# We now subset resistance and recovery data that are declining 
# First we subset populations that are declining

dec_pop <- pops_data %>% filter(mu<0)

# Now subset using the IDs for each 

rec_data_dec <- rec_data %>% filter(ID%in%dec_pop$ID) # Recovery
res_data_dec <- res_neg %>% filter(ID%in%dec_pop$ID) # Resistance

# Number of threats on recovery ------------------------------------------------

m1 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat-1 + (1|SpeciesName),
          iter = iter, thin = thin, warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = rec_data_dec, family = gaussian,
          cores = 20)

# System -----------------------------------------------------------------------

m2 <- brm(mu | se(tau, sigma = TRUE) ~ System-1 + (1|SpeciesName),
          iter = iter, thin = thin, warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = rec_data_dec, family = gaussian,
          cores = 20)

# Taxa -------------------------------------------------------------------------

m3 <- brm(mu | se(tau, sigma = TRUE) ~ Taxon-1 + (1|SpeciesName),
          iter = iter, thin = thin, warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = rec_data_dec, family = gaussian,
          cores = 20)

# Number of threats and system -------------------------------------------------

mi1 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:System - 1 + (1 | SpeciesName),
           iter = iter, thin = thin,  warmup = warmup,
           prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
           data = rec_data_dec,
           family = gaussian,
           cores = 20)

# Number of threats and taxa ---------------------------------------------------

mi2 <-  brm(mu | se(tau, sigma = TRUE) ~ n.threat:Taxon - 1 + (1 | SpeciesName),
            iter = iter,
            thin = thin,
            warmup = warmup,
            prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
            data = rec_data_dec,
            family = gaussian,
            cores = 20)

# Number of threats ------------------------------------------------------------

m4 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat-1 + (1|SpeciesName) ,
          iter = iter, thin = thin, warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = res_data,
          family = gaussian,
          cores = 20)

# System -----------------------------------------------------------------------

m5 <- brm(mu | se(tau, sigma = TRUE) ~ System - 1 + (1|SpeciesName) ,
          iter = iter,  thin = thin,  warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = res_data_dec,
          family = gaussian,
          cores = 20)

# Taxa -------------------------------------------------------------------------

m6 <- brm(mu | se(tau, sigma = TRUE) ~ Taxon - 1 + (1|SpeciesName) ,
          iter = iter, thin = thin,  warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = res_data_dec,
          family = gaussian,
          cores = 20)

# Number of threats and system -------------------------------------------------

mi3 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:System - 1 + (1|SpeciesName) ,
           iter = iter, thin = thin,  warmup = warmup,
           prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
           data = res_data_dec,
           family = gaussian,
           cores = 20)

# Number of threats and Taxa ---------------------------------------------------

mi4 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:Taxon - 1 + (1|SpeciesName) ,
           iter = iter, thin = thin, warmup = warmup,
           prior= prior, 
           control = list(adapt_delta = .975, max_treedepth = 20),
           data = res_data_dec,
           family = gaussian,
           cores = 20)

# Save the models ##############################################################

setwd(DataPath)
save(file = "SuppAnal2.RData", 
     m1,m2,m3,m4,m5,m6,
     mi1, mi2, mi3,mi4, 
     res_data_dec, rec_data_dec)