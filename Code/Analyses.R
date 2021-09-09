# -----------------------------------------------------------------------------#
# - FILE NAME:   Analyses.R
# - DATE:        21/09/2020
# - DESCRIPTION: Use the population growth (r) trends to analyse the factors 
#               determining the resistance and recovery trend over time.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -----------------------------------------------------------------------------#

rm(list = ls(all = TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(brms)
library(data.table)
library(zoo)
library(dplyr)
library(expss)
library(MASS)

# Set working directory

path <- gsub("LPI/Code", "LPI/", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path, "Code")
DataPath <-  paste0(path, "Data")
ResultPath <-  paste0(path, "Results")

# Recovery #####################################################################

# Load data 

load(paste0(DataPath, "/RecoveryChange.RData"))

#Change tau for the analyses

colnames(rec_data)[15] <- "tau"

# SECTION 1: EXPLORE THE DATA #####################

# Number of populations #3347
length(unique(rec_data$ID))

length(unique(rec_data$ID[!is.na(rec_data$bs_g)]))

# Number of species #1341

length(unique(rec_data$SpeciesName))
length(unique(rec_data$SpeciesName[!is.na(rec_data$bs_g)]))

# Proportion of data per groups

table1 <-rec_data %>% 
  distinct(ID, .keep_all=T) %>% 
  group_by(System, Taxon, n.threat) %>% 
  tally() %>% 
  mutate(per= (n/sum(n))*100)

# Proportion increasing and decreasing

rec_data %>% 
  mutate(incr= ifelse(mu>0,
                      "Increase",
                      ifelse(mu<0, "Decrease", "Stable"))) %>% 
  group_by(incr) %>% 
  tally() %>% 
  mutate(per= (n/sum(n))*100)

# SECTION 2: MODELLING ######################################

#Set modelling parameters

iter = 10000
thin = 0.0005*iter
warmup = 0.1*iter

# Set weakly informed prior

prior <- c(prior(normal(0, 1), class = b),
           prior(exponential(1), class = sigma))

# Transform to a data frame

rec_data <- as.data.frame(rec_data)

# Transform number of threats in factor

rec_data <- rec_data %>% 
  mutate(n.threat=as.factor(n.threat),
         n.threat = factor(n.threat, levels = c("0", "1","2","3")),
         System=as.factor(System))

# Number of threats ------------------------------------------------------------

m1 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat-1 + (1|SpeciesName),
          iter = iter, thin = thin, warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = rec_data, family = gaussian,
          cores = 20)


# System -----------------------------------------------------------------------

m2 <- brm(mu | se(tau, sigma = TRUE) ~ System-1 + (1|SpeciesName),
          iter = iter, thin = thin, warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = rec_data, family = gaussian,
          cores = 20)

# Taxonomic group --------------------------------------------------------------

m3 <- brm(mu | se(tau, sigma = TRUE) ~ Taxon-1 + (1|SpeciesName),
          iter = iter, thin = thin, warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = rec_data, family = gaussian,
          cores = 20)

# Interactive models recovery ##################################################

# Number of threats and System -------------------------------------------------

mi1 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:System - 1 + (1 | SpeciesName),
  iter = iter,
  thin = thin,
  warmup = warmup,
  prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
  data = rec_data,
  family = gaussian,
  cores = 20)

# Number of threats and taxon --------------------------------------------------

mi2 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:Taxon - 1 + (1 | SpeciesName),
    iter = iter,
    thin = thin,
    warmup = warmup,
    prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
    data = rec_data,
    family = gaussian,
    cores = 20)

# Resistance ###################################################################

# Load data

load(paste0(DataPath, "/ResistanceChange.RData"))

#Change tau for the analyses

colnames(res_data)[15] <- "tau"

# SECTION 1: EXPLORE THE DATA #####################

# Number of populations #3347

length(unique(res_data$ID))

# Number of species #1341

length(unique(res_data$SpeciesName))

# Proportion of data per groups

table2 <- res_data %>% 
  distinct(ID, .keep_all=T) %>% 
  group_by(System, Taxon, n.threat) %>% 
  tally() %>% 
  mutate(per= (n/sum(n))*100)

# Join with previous table

table <- rbind(table1, table2)

# Save into table S1

setwd(ResultPath)
write.csv(table, "TableS1.csv")

# Proportion of increase/decrease

res_data %>% 
  mutate(incr= ifelse(mu>0,
                      "Increase",
                      ifelse(mu<0, "Decrease", "Stable"))) %>% 
  group_by(incr) %>% 
  tally() %>% 
  mutate(per= (n/sum(n))*100)

# SECTION 2: MODELLING ######################################

#Set modelling parameters

iter = 10000
thin = 0.0005 * iter
warmup = 0.1 * iter

# Set weakly informed prior

prior <- c(prior(normal(0, 1), class = b),
           prior(exponential(1), class = sigma))

# Prepare the data
# Transform to a data frame

res_data <- as.data.frame(res_data)

# Transform number of threats in factor

res_data <- res_data %>%
  mutate(n.threat = as.factor(n.threat),
    n.threat = factor(n.threat, levels = c("0", "1", "2", "3")),
    System = as.factor(System))

# Number of threats ------------------------------------------------------------

m4 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat - 1 + (1|SpeciesName),
          iter = iter,
          thin = thin,
          warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = res_data,
          family = gaussian,
          cores = 20)

# System -----------------------------------------------------------------------

m5 <- brm(mu | se(tau, sigma = TRUE) ~ System - 1 + (1|SpeciesName) ,
          iter = iter,
          thin = thin,
          warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = res_data,
          family = gaussian,
          cores = 20)

# Taxon ------------------------------------------------------------------------

# fit the model

m6 <- brm(mu | se(tau, sigma = TRUE) ~ Taxon - 1 + (1|SpeciesName),
          iter = iter,
          thin = thin,
          warmup = warmup,
          prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
          data = res_data,
          family = gaussian,
          cores = 20)

# Interactive models ###########################################################

# Number of threats and System -------------------------------------------------

mi3 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:System - 1 + (1|SpeciesName) ,
           iter = iter,
           thin = thin,
           warmup = warmup,
           prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
           data = res_data,
           family = gaussian,
           cores = 20)

# Number of threats and Taxon --------------------------------------------------

mi4 <- brm(mu | se(tau, sigma = TRUE) ~ n.threat:Taxon - 1 + (1|SpeciesName) ,
           iter = iter,
           thin = thin,
           warmup = warmup,
           prior= prior, control = list(adapt_delta = .975, max_treedepth = 20),
           data = res_data,
           family = gaussian,
           cores = 20)

# Save all the models ##########################################################

setwd(ResultsPath)

save(file = "Factors.RData",
     m1,m2,m3,m4, m5, m6, 
     mi1, mi2, mi3, mi4)
