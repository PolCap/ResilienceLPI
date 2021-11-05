# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupAnalyses.R
# - DATE:        20/11/2020
# - DESCRIPTION: Supplementary analyses using the population growth (r) trends to 
#                analyse the other factors determining the resilience change
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
library(bayestestR)
library(phytools)
library(geiger)


# Set working directory

path <- gsub("LPI/Code", "LPI/", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path, "Code")
DataPath <-  paste0(path, "Data")
ResultPath <-  paste0(path, "Results")

# Load data

load(paste0(DataPath, "/ResistanceChange.RData"))
load(paste0(DataPath, "/RecoveryChange.RData"))

# Resistance ###################################################################
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
res_data_dec <- res_data %>% filter(ID%in%dec_pop$ID) # Resistance

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

# Phylogenetic Signal ##########################################################

# Load resistance data

load(paste0(DataPath, "/ResistanceChange.RData"))

# Load recovery data

load(paste0(DataPath, "/RecoveryChange.RData"))

# Read the trees 

setwd(DataPath)
amp_tree <- read.nexus("amphibians.nex")
bird_tree <- read.nexus("birds.nex")
mam_tree <- read.nexus("mammals.nex")
rep_tree <- read.nexus("reptiles.nex")
fish_tree <- read.tree("fishtree.tre")

# Subset the first tree from the list (run just once)

amp_tree <- amp_tree[[1]]
bird_tree <- bird_tree[[1]]
mam_tree <- mam_tree[[1]]
rep_tree <- rep_tree[[1]]

# Correct species names 

amp_tree$tip.label <-  gsub("_", " ", amp_tree$tip.label)
bird_tree$tip.label <-  gsub("_", " ", bird_tree$tip.label)
mam_tree$tip.label <-  gsub("_", " ", mam_tree$tip.label)
rep_tree$tip.label <-  gsub("_", " ", rep_tree$tip.label)
fish_tree$tip.label <-  gsub("_", " ", fish_tree$tip.label)

# Subset species names for both datasets (resistance and recovery)

# Amphibians

amp_res <- res_data %>% filter(Taxon=="Amphibians")
amp_rec <- rec_data %>% filter(Taxon=="Amphibians")

# Birds 

bird_res <- res_data %>% filter(Taxon=="Birds") 
bird_rec <- rec_data %>% filter(Taxon=="Birds") 

# Mammals 

mam_res <- res_data %>% filter(Taxon=="Mammals")
mam_rec <- rec_data %>% filter(Taxon=="Mammals")

# Reptiles

rep_res <- res_data %>% filter(Taxon=="Reptiles")
rep_rec <- rec_data %>% filter(Taxon=="Reptiles") 

# Fish

fish_res <- res_data %>% filter(Taxon=="Fish")
fish_rec <- rec_data %>% filter(Taxon=="Fish")

# We prune the trees and the datasets

# Amphibians
# Resistance
amp_sp_res <- unique(amp_res$SpeciesName)
names(amp_sp_res) <- amp_sp_res
(chk<- name.check(amp_tree, amp_sp_res))
amp_tree_res <- drop.tip(amp_tree, chk$tree_not_data)
amp_res <- amp_res[amp_res$SpeciesName%in%amp_tree$tip.label,]

#Recovery
amp_sp_rec <- unique(amp_rec$SpeciesName)
names(amp_sp_rec) <- amp_sp_rec
(chk<- name.check(amp_tree, amp_sp_rec))
amp_tree_rec <- drop.tip(amp_tree, chk$tree_not_data)
amp_rec <- amp_rec[amp_rec$SpeciesName%in%amp_tree$tip.label,]

# Birds
# Resistance
bird_sp_res <- unique(bird_res$SpeciesName)
names(bird_sp_res) <- bird_sp_res
(chk<- name.check(bird_tree, bird_sp_res))
bird_tree_res <- drop.tip(bird_tree, chk$tree_not_data)
bird_res <- bird_res[bird_res$SpeciesName%in%bird_tree$tip.label,]

#Recovery
bird_sp_rec <- unique(bird_rec$SpeciesName)
names(bird_sp_rec) <- bird_sp_rec
(chk<- name.check(bird_tree, bird_sp_rec))
bird_tree_rec <- drop.tip(bird_tree, chk$tree_not_data)
bird_rec <- bird_rec[bird_rec$SpeciesName%in%bird_tree$tip.label,]

# Mammals
# Resistance
mam_sp_res <- unique(mam_res$SpeciesName)
names(mam_sp_res) <- mam_sp_res
(chk<- name.check(mam_tree, mam_sp_res))
mam_tree_res <- drop.tip(mam_tree, chk$tree_not_data)
mam_res <- mam_res[mam_res$SpeciesName%in%mam_tree$tip.label,]

#Recovery
mam_sp_rec <- unique(mam_rec$SpeciesName)
names(mam_sp_rec) <- mam_sp_rec
(chk<- name.check(mam_tree, mam_sp_rec))
mam_tree_rec <- drop.tip(mam_tree, chk$tree_not_data)
mam_rec <- mam_rec[mam_rec$SpeciesName%in%mam_tree$tip.label,]

# Reptiles 
# Resistance
rep_sp_res <- unique(rep_res$SpeciesName)
names(rep_sp_res) <- rep_sp_res
(chk<- name.check(rep_tree, rep_sp_res))
rep_tree_res <- drop.tip(rep_tree, chk$tree_not_data)
rep_res <- rep_res[rep_res$SpeciesName%in%rep_tree$tip.label,]

#Recovery
rep_sp_rec <- unique(rep_rec$SpeciesName)
names(rep_sp_rec) <- rep_sp_rec
(chk<- name.check(rep_tree, rep_sp_rec))
rep_tree_rec <- drop.tip(rep_tree, chk$tree_not_data)
rep_rec <- rep_rec[rep_rec$SpeciesName%in%rep_tree$tip.label,]

# Fishes
# Resistance
fish_sp_res <- unique(fish_res$SpeciesName)
names(fish_sp_res) <- fish_sp_res
(chk<- name.check(fish_tree, fish_sp_res))
fish_tree_res <- drop.tip(fish_tree, chk$tree_not_data)
fish_res <- fish_res[fish_res$SpeciesName%in%fish_tree$tip.label,]

#Recovery
fish_sp_rec <- unique(fish_rec$SpeciesName)
names(fish_sp_rec) <- fish_sp_rec
(chk<- name.check(fish_tree, fish_sp_rec))
fish_tree_rec <- drop.tip(fish_tree, chk$tree_not_data)
fish_rec <- fish_rec[fish_rec$SpeciesName%in%fish_tree$tip.label,]

# Make each tree ultrametric

amp_tree_res <-  chronos(amp_tree_res)
amp_tree_rec <-  chronos(amp_tree_rec)
bird_tree_res <-  chronos(bird_tree_res)
bird_tree_rec <-  chronos(bird_tree_rec)
mam_tree_res <-  chronos(mam_tree_res)
mam_tree_rec <-  chronos(mam_tree_rec)
rep_tree_res <-  chronos(rep_tree_res)
rep_tree_rec <-  chronos(rep_tree_rec)
fish_tree_res <-  chronos(fish_tree_res)
fish_tree_rec <-  chronos(fish_tree_rec)

# Add animal as column

amp_res <- amp_res %>% mutate(animal=SpeciesName)
amp_rec <- amp_rec %>% mutate(animal=SpeciesName)
bird_res <- bird_res %>% mutate(animal=SpeciesName)
bird_rec <- bird_rec %>% mutate(animal=SpeciesName)
mam_res <- mam_res %>% mutate(animal=SpeciesName)
mam_rec <- mam_rec %>% mutate(animal=SpeciesName)
rep_res <- rep_res %>% mutate(animal=SpeciesName)
rep_rec <- rep_rec %>% mutate(animal=SpeciesName)
fish_res <- fish_res %>% mutate(animal=SpeciesName)
fish_rec <- fish_rec %>% mutate(animal=SpeciesName)

# Set modelling parameters

iter <- 10000
thin <- 0.0005 * iter
warmup <- 0.1 * iter

# Set priors

priors <- c(prior(normal(0, 10), class = Intercept),
            prior(exponential(1), class = sd))

# Amphibians -------------------------------------------------------------------

# Resistance

# Define the correlation structure

A <- ape::vcv.phylo(amp_tree_res)

# Model

ps_amp_res <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                  data = amp_res,
                  family = gaussian(), 
                  data2 = list(A = A),
                  cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_res_amp <- hypothesis(ps_amp_res, hyp, class = NULL))


# Recovery

# Define the correlation structure

A <- ape::vcv.phylo(amp_tree_rec)

# Model

ps_amp_rec <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                  data = amp_rec,
                  family = gaussian(), 
                  data2 = list(A = A),
                  cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_rec_amp <- hypothesis(ps_amp_rec, hyp, class = NULL))

# Birds ------------------------------------------------------------------------

# Resistance

# Define the correlation structure

A <- ape::vcv.phylo(bird_tree_res)

# Model

ps_bird_res <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                   data = bird_res,
                   family = gaussian(), 
                   data2 = list(A = A),
                   cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_res_bird <- hypothesis(ps_bird_res, hyp, class = NULL))


# Recovery

# Define the correlation structure

A <- ape::vcv.phylo(bird_tree_rec)

# Model

ps_bird_rec <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                   data = bird_rec,
                   family = gaussian(), 
                   data2 = list(A = A),
                   cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_rec_bird <- hypothesis(ps_bird_rec, hyp, class = NULL))

# Mammals ----------------------------------------------------------------------

# Resistance

# Define the correlation structure

A <- ape::vcv.phylo(mam_tree_res)

# Model

ps_mam_res <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                  data = mam_res,
                  family = gaussian(), 
                  data2 = list(A = A),
                  cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_res_mam <- hypothesis(ps_mam_res, hyp, class = NULL))


# Recovery

# Define the correlation structure

A <- ape::vcv.phylo(mam_tree_rec)

# Model

ps_mam_rec <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                  data = mam_rec,
                  family = gaussian(), 
                  data2 = list(A = A),
                  cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_rec_mam <- hypothesis(ps_mam_rec, hyp, class = NULL))

# Reptiles ----------------------------------------------------------------------

# Resistance

# Define the correlation structure

A <- ape::vcv.phylo(rep_tree_res)

# Model

ps_rep_res <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                  data = rep_res,
                  family = gaussian(), 
                  data2 = list(A = A),
                  cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_res_rep <- hypothesis(ps_rep_res, hyp, class = NULL))


# Recovery

# Define the correlation structure

A <- ape::vcv.phylo(rep_tree_rec)

# Model

ps_rep_rec <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                  data = rep_rec,
                  family = gaussian(), 
                  data2 = list(A = A),
                  cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_rec_rep <- hypothesis(ps_rep_rec, hyp, class = NULL))

# Fishes -----------------------------------------------------------------------

# Resistance

# Define the correlation structure

A <- ape::vcv.phylo(fish_tree_res)

# Model

ps_fish_res <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                   data = fish_res,
                   family = gaussian(), 
                   data2 = list(A = A),
                   cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_res_fish <- hypothesis(ps_fish_res, hyp, class = NULL))


# Recovery

# Define the correlation structure

A <- ape::vcv.phylo(fish_tree_rec)

# Model

ps_fish_rec <- brm(mu ~ 1 + (1|gr(animal, cov = A)) + (1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
                   data = fish_rec,
                   family = gaussian(), 
                   data2 = list(A = A),
                   cores = 20)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_SpeciesName__Intercept^2 + sigma^2) = 0")

(psign_rec_fish <- hypothesis(ps_fish_rec, hyp, class = NULL))

# Save the results -------------------------------------------------------------

setwd(ResultPath)
save(ps_amp_rec, ps_amp_res,ps_bird_rec,ps_bird_res, 
     ps_mam_rec,ps_mam_res,ps_rep_rec,ps_rep_res, 
     ps_fish_rec, ps_fish_res, file="PSModels.RData")
