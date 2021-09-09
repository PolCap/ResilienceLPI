# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupFigures.R         
# - DATE:        03/08/2021
# - DESCRIPTION: Supplementary figures 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(tidybayes)
library(tidyverse)
library(rstan)
library(rstanarm)
library(brms)
library(ggplot2) #Now I load the latest version of ggplot2
library(bayesplot)
library(bayestestR)
library(ggdist)
library(magrittr)
library(forcats)
library(modelr)
library(emmeans)
library(cowplot)
library(rphylopic)
library(RCurl)
library(png)
library(mapdata)

# Set default ggplot theme

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=15, margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                  axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color="black", size = 12),
                  axis.text.y = element_text(color="black", size = 12),
                  strip.text.x = element_text(size = 12),
                  axis.ticks = element_line(color="black")))


# Set working directory

path <- gsub("LPI/Code", "LPI/", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path, "Code")
DataPath <-  paste0(path, "Data")
ResultPath <-  paste0(path, "Results")

# Figure S2: General model outcomes ############################################

# Load models

load(paste0(ResultPath, "/SuppAnal2.RData"))

# Panel a: negative lambda vs number of stressors ------------------------------

# Sample size

sample_size <- res_data_dec %>% 
  group_by(n.threat) %>% 
  summarise(n=n()) %>% 
  mutate(n.threat=factor(n.threat))

med <- m4 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_n.threat", "", .variable)) %>% 
  group_by(.variable) %>% 
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med, by=c("n.threat"=".variable"))

# Plot intervals 

dat <- m4 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_n.threat", "", .variable))

(g1 <- m4 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(.variable = gsub("b_n.threat", "", .variable)) %>% 
    median_qi(.width = c(.95, .8, .5)) %>%
    ggplot(aes(y = .variable, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, alpha=0.2, aes(fill=.variable)) +
    scale_fill_manual(values=c("#EDA30C", "#F7910C", "#E15F00", "#F74C0C"))+
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=n.threat, 
                                      label = paste0("n=", n)),
              vjust   = -1.5, hjust=-.1)+
    xlim(-0.03, 0.03)+
    labs(x="Posterior estimate", y = "Number of stressors") +
    theme(legend.position = "none"))

# Panel d: positive lambda vs number of stressors ------------------------------

# Sample size

sample_size <- rec_data_dec %>% 
  group_by(n.threat) %>% 
  summarise(n=n()) %>% 
  mutate(n.threat=factor(n.threat))

med <- m1 %>%
  gather_draws(b_n.threat0, b_n.threat1, b_n.threat2, b_n.threat3) %>%
  mutate(.variable = gsub("b_n.threat", "", .variable)) %>% 
  group_by(.variable) %>% 
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med, by=c("n.threat"=".variable"))

# Plot intervals 

dat <- m1 %>%
  gather_draws(b_n.threat0, b_n.threat1, b_n.threat2, b_n.threat3) %>%
  mutate(.variable = gsub("b_n.threat", "", .variable))

(g4 <- m1 %>%
    gather_draws(b_n.threat0, b_n.threat1, b_n.threat2, b_n.threat3) %>%
    mutate(.variable = gsub("b_n.threat", "", .variable)) %>% 
    median_qi(.width = c(.95, .8, .5)) %>%
    ggplot(aes(y = .variable, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, alpha=0.2, aes(fill=.variable)) +
    scale_fill_manual(values=c("#EDA30C", "#F7910C", "#E15F00", "#F74C0C"))+
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=n.threat, 
                                      label = paste0("n=", n)),
              vjust   = -1.5, hjust=-.1)+
    xlim(-0.03, 0.03)+
    labs(x="Posterior estimate", y = "Number of stressors") +
    theme(legend.position = "none"))

# Panel b: System effects on negative lambda ----------------------------------- 
# Sample size

sample_size <- res_data_dec %>% 
  group_by(System) %>% 
  summarise(n=n()) 

med <- m5 %>%
  gather_draws(b_SystemFreshwater, b_SystemMarine, b_SystemTerrestrial) %>%
  mutate(.variable = gsub("b_System", "", .variable)) %>% 
  group_by(.variable) %>% 
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med, by=c("System"=".variable"))

# Plot intervals 

dat <- m5 %>%
  gather_draws(b_SystemFreshwater, b_SystemMarine, b_SystemTerrestrial) %>%
  mutate(.variable = gsub("b_System", "", .variable)) 

(g2 <- m5 %>%
    gather_draws(b_SystemFreshwater, b_SystemMarine, b_SystemTerrestrial) %>%
    mutate(.variable = gsub("b_System", "", .variable)) %>% 
    median_qi(.width = c(.95, .8, .5)) %>%
    ggplot(aes(y = .variable, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, aes(fill=.variable), alpha=0.2) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=System, 
                                      label = paste0("n=", n)),
              vjust   = -1.5, hjust =-.1)+
    xlim(-0.020,0.020)+
    scale_fill_manual("", values = c("#A1D6E2","#336B87", "#CC954E"))+
    labs(x="Posterior estimate", y = "System") + 
    theme(legend.position = "none"))

# Panel e: System effects on positive lambda ----------------------------------- 

# Sample size

sample_size <- rec_data_dec %>% 
  group_by(System) %>% 
  summarise(n=n()) 

med <- m2 %>%
  gather_draws(b_SystemFreshwater, b_SystemMarine, b_SystemTerrestrial) %>%
  mutate(.variable = gsub("b_System", "", .variable)) %>% 
  group_by(.variable) %>% 
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med, by=c("System"=".variable"))

# Plot intervals 

dat <- m2 %>%
  gather_draws(b_SystemFreshwater, b_SystemMarine, b_SystemTerrestrial) %>%
  mutate(.variable = gsub("b_System", "", .variable)) 

(g5 <- m2 %>%
    gather_draws(b_SystemFreshwater, b_SystemMarine, b_SystemTerrestrial) %>%
    mutate(.variable = gsub("b_System", "", .variable)) %>% 
    median_qi(.width = c(.95, .8, .5)) %>%
    ggplot(aes(y = .variable, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, aes(fill=.variable), alpha=0.2) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=System, 
                                      label = paste0("n=", n)),
              vjust   = -1.5, hjust=-0.1)+
    xlim(-0.026,0.026)+
    scale_fill_manual("", values = c("#A1D6E2","#336B87", "#CC954E"))+
    labs(x="Posterior estimate", y = "System") + 
    theme(legend.position = "none"))

# Panel c: Taxon effects on negative lambda ----------------------------------- 

# Sample size

sample_size <- res_data_dec %>% 
  group_by(Taxon) %>% 
  summarise(n=n()) 

med <- m6 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_Taxon", "", .variable)) %>%  
  group_by(.variable) %>% 
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med, by=c("Taxon"=".variable"))

# Plot intervals 

dat <- m6 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_Taxon", "", .variable)) 

(g3 <- m6 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(.variable = gsub("b_Taxon", "", .variable))  %>%    
    median_qi(.width = c(.95, .8, .5)) %>%
    ggplot(aes(y = .variable, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, aes(fill=.variable), alpha=0.2) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=Taxon, 
                                      label = paste0("n=", n)),
              vjust   = -1.5,
              hjust=-0.1)+
    xlim(-0.048,0.048)+
    scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    labs(x="Posterior estimate", y = "Taxon") + 
    theme(legend.position = "none"))

# Panel f: Taxon effects on positive lambda ----------------------------------- 

# Sample size

sample_size <- rec_data_dec %>% 
  group_by(Taxon) %>% 
  summarise(n=n()) 

med <- m3 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_Taxon", "", .variable)) %>%  
  group_by(.variable) %>% 
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med, by=c("Taxon"=".variable"))

# Plot intervals 

dat <- m3 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_Taxon", "", .variable)) 

(g6 <- m3 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(.variable = gsub("b_Taxon", "", .variable))  %>%    
    median_qi(.width = c(.95, .8, .5)) %>%
    ggplot(aes(y = .variable, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, aes(fill=.variable), alpha=0.2) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=Taxon, 
                                      label = paste0("n=", n)),
              vjust   = -1.5,
              hjust=-0.1)+
    xlim(-0.048,0.048)+
    scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    labs(x="Posterior estimate", y = "Taxon") + 
    theme(legend.position = "none"))

# Combined Figure --------------------------------------------------------------

ggf <- plot_grid(g2+ xlab("")+ 
                   theme(plot.margin =  margin(0.5, 0, 0, 0, "cm"),
                         axis.title.y = element_text(margin = margin(0, -5, 0, 0))),
                 g3 + xlab("")+ 
                   theme(plot.margin =  margin(0.5, 0, 0, 0, "cm"),
                         axis.title.y = element_text(margin = margin(0, -5, 0, 0))),
                 g1+ xlab("")+ 
                   theme(plot.margin =  margin(0.5, 0, 0, 0, "cm"),
                         axis.title.y = element_text(margin = margin(0, 5, 0, 5))), 
                 
                 g5 + theme(plot.margin =  margin(0.5, 0, 0, 0, "cm"),
                            axis.title.y = element_text(margin = margin(0, -5, 0, 0))),
                 g6 + theme(plot.margin =  margin(0.5, 0, 0, 0, "cm"),
                            axis.title.y = element_text(margin = margin(0, -5, 0, 0))),
                 labels = "auto",
                 g4+ theme(plot.margin =  margin(0.5, 0, 0, 0, "cm"),
                           axis.title.y = element_text(margin = margin(0, 5, 0, 5))),
                 nrow = 2,
                 align = "h", axis = "b") 

# Add our customised titles 

t1 <- ggdraw() + draw_label("Resistance",angle = -90, 
                            fontface='bold', size=18)
t2 <- ggdraw() + draw_label("Recovery",angle = -90, 
                            fontface='bold',size = 18)
titles <- plot_grid(t1,t2, ncol = 1)

(figureS2<- plot_grid(ggf, titles,  
                      ncol=2, rel_widths =c(1, 0.02)))


ggsave("FigureS2.pdf", figureS2, 
       width = 12, height = 10,
       path = ResultPath)

# Figure S3: Interactions System ###############################################

# Sample size

sample_size <- res_data_dec %>% 
  group_by(n.threat, System) %>% 
  summarise(n=n()) %>% 
  mutate(n.threat=as.character(n.threat))

med <- mi3 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         System = gsub("System", "", .variable), 
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         System =gsub(":", "", System),
         System =gsub("b_n.threat", "", System),
         System = gsub('[[:digit:]]+', '', System)) %>%   
  group_by(n.threat, System) %>% 
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med)

# Resistance vs system

dat <- mi3 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         System = gsub("System", "", .variable), 
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         System =gsub(":", "", System),
         System =gsub("b_n.threat", "", System),
         System = gsub('[[:digit:]]+', '', System))   

(g1 <- mi3 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.width = c(.95, .8, .5)) %>%
    mutate(n.threat= gsub("[^0-9.-]", "", .variable),
           System = gsub("System", "", .variable), 
           n.threat= gsub("[^[:alnum:]]", "", n.threat),
           System =gsub(":", "", System),
           System =gsub("b_n.threat", "", System),
           System = gsub('[[:digit:]]+', '', System)) %>%   
    ggplot(aes(y = n.threat, x = .value)) +
    geom_vline(xintercept = 0, 
               linetype = "dashed", 
               colour="grey50") +
    stat_slab(data=dat, aes(fill=System, alpha=n.threat)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=n.threat, 
                                      label = paste0("n=", n)),
              vjust   = -1.5,
              hjust=-0.1)+
    scale_fill_manual(values =c("#A1D6E2", 
                                "#387594", 
                                "#CC954E"))+ 
    scale_alpha_discrete(range = c(0.2,0.7))+
    labs(x="Posterior estimate", y = "Number of stressors") +
    facet_wrap(.~System) +
    xlim(-0.045, 0.045) +
    theme(legend.position = "none"))

# Recovery vs system

# Sample size

sample_size <- rec_data_dec %>% 
  group_by(n.threat, System) %>% 
  summarise(n=n()) %>% 
  mutate(n.threat=as.character(n.threat))

med <- mi1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         System = gsub("System", "", .variable), 
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         System =gsub(":", "", System),
         System =gsub("b_n.threat", "", System),
         System = gsub('[[:digit:]]+', '', System)) %>%   
  group_by(n.threat, System) %>% 
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med)

# Resistance vs system

dat <- mi1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         System = gsub("System", "", .variable), 
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         System =gsub(":", "", System),
         System =gsub("b_n.threat", "", System),
         System = gsub('[[:digit:]]+', '', System))   

(g2 <- mi1 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.width = c(.95, .8, .5)) %>%
    mutate(n.threat= gsub("[^0-9.-]", "", .variable),
           System = gsub("System", "", .variable), 
           n.threat= gsub("[^[:alnum:]]", "", n.threat),
           System =gsub(":", "", System),
           System =gsub("b_n.threat", "", System),
           System = gsub('[[:digit:]]+', '', System)) %>%   
    ggplot(aes(y = n.threat, x = .value)) +
    geom_vline(xintercept = 0, 
               linetype = "dashed", 
               colour="grey50") +
    stat_slab(data=dat, aes(fill=System, alpha=n.threat)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=n.threat, 
                                      label = paste0("n=", n)),
              vjust   = -1.5,
              hjust=-0.1)+
    scale_fill_manual(values =c("#A1D6E2", 
                                "#387594", 
                                "#CC954E"))+ 
    scale_alpha_discrete(range = c(0.2,0.7))+
    labs(x="Posterior estimate", y = "Number of stressors") +
    facet_wrap(.~System) +
    xlim(-0.052, 0.052) +
    theme(legend.position = "none"))

# Combined Figure --------------------------------------------------------------

ggf <- plot_grid(g1+xlab(""), g2 + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()), 
  labels = "auto",
  nrow = 2,
  align = "h", axis = "b") 

# Add our customised titles 

t1 <- ggdraw() + draw_label("Resistance",angle = -90, 
                            fontface='bold', size=18)
t2 <- ggdraw() + draw_label("Recovery",angle = -90, 
                            fontface='bold',size = 18)
titles <- plot_grid(t1,t2, ncol = 1)

(figure3<- plot_grid(ggf, titles,  
                     ncol=2, rel_widths =c(1, 0.02)))

ggsave("FigureS3.pdf", figure3, 
       width = 12, height = 10,
       path = ResultPath)

# Figure S4: Interaction with Taxa #############################################

# Here we define the limits of the x scale for each faced

limits_data <- data.frame(Taxon = c("Amphibians", "Amphibians", 
                                    "Birds", "Birds", 
                                    "Fish", "Fish", 
                                    "Mammals", "Mammals",
                                    "Reptiles", "Reptiles"), 
                          y = "0", 
                          x = c(-0.15,0.15, 
                                -0.045, 0.045, 
                                -0.08, 0.08,
                                -0.08, 0.08,
                                -0.08, 0.08))

# Sample size

sample_size <- res_data_dec %>% 
  group_by(n.threat, Taxon) %>% 
  summarise(n=n()) %>% 
  mutate(n.threat=as.character(n.threat))

med <- mi2 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         Taxon = gsub("Taxon", "", .variable), 
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_n.threat", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon))  %>%   
  group_by(n.threat, Taxon) %>% 
  summarise(med=quantile(.value,probs = 0.9))

sample_size<- sample_size %>% 
  left_join(med)

# Plot intervals 

dat <- mi2 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         Taxon = gsub("Taxon", "", .variable), 
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_n.threat", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon)) 

(g1 <- mi2 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(n.threat= gsub("[^0-9.-]", "", .variable),
           n.threat= gsub("[^[:alnum:]]", "", n.threat),
           Taxon = gsub("Taxon", "", .variable), 
           Taxon =gsub(":", "", Taxon),
           Taxon =gsub("b_n.threat", "", Taxon),
           Taxon = gsub('[[:digit:]]+', '', Taxon)) %>%
    ggplot(aes(y = n.threat, x = .value)) +
    facet_wrap(Taxon~., nrow = 1,scales = "free_x")+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, aes(fill=Taxon, alpha=n.threat)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=n.threat, 
                                      label = paste0("n=", n)),
              vjust   = -1.5,
              hjust=-0.1)+
    scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    scale_alpha_discrete(range = c(0.2,0.7))+
    geom_blank(data = limits_data, aes(x = x, y = y)) + 
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    #xlim(-0.15,0.15)+
    labs(x="Posterior estimate", y = "Number of stressors") + 
    theme(legend.position = "none"))

# Recovery vs Taxon

# Sample size

sample_size <- rec_data_dec %>% 
  group_by(n.threat, Taxon) %>% 
  summarise(n=n()) %>% 
  mutate(n.threat=as.character(n.threat))

med <- mi4 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         Taxon = gsub("Taxon", "", .variable), 
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_n.threat", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon))  %>%   
  group_by(n.threat, Taxon) %>% 
  summarise(med=quantile(.value,probs = 0.9))

sample_size<- sample_size %>% 
  left_join(med)

# Plot it

dat <- mi4 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         Taxon = gsub("Taxon", "", .variable), 
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_n.threat", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon)) 

(g2 <- mi4 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(n.threat= gsub("[^0-9.-]", "", .variable),
           n.threat= gsub("[^[:alnum:]]", "", n.threat),
           Taxon = gsub("Taxon", "", .variable), 
           Taxon =gsub(":", "", Taxon),
           Taxon =gsub("b_n.threat", "", Taxon),
           Taxon = gsub('[[:digit:]]+', '', Taxon)) %>%
    ggplot(aes(y = n.threat, x = .value)) +
    facet_wrap(Taxon~., nrow = 1, scales = "free_x")+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, aes(fill=Taxon, alpha=n.threat)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=n.threat, 
                                      label = paste0("n=", n)),
              vjust   = -1.5,
              hjust=-0.1)+
    scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    scale_alpha_discrete(range = c(0.2,0.7))+
    geom_blank(data = limits_data, aes(x = x, y = y)) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    #xlim(-0.15,0.15)+
    labs(x="Posterior estimate", y = "Number of stressors") + 
    theme(legend.position = "none"))


# Combined Figure --------------------------------------------------------------

ggf <- plot_grid(g1+xlab(""), g2 + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()), 
  labels = "auto",
  nrow = 2,
  align = "h", axis = "b") 

# Add our customised titles 

t1 <- ggdraw() + draw_label("Resistance",angle = -90, 
                            fontface='bold', size=18)
t2 <- ggdraw() + draw_label("Recovery",angle = -90, 
                            fontface='bold',size = 18)
titles <- plot_grid(t1,t2, ncol = 1)

(ggf4<- plot_grid(ggf, titles,  
                  ncol=2, rel_widths =c(1, 0.02)))

# Add the silhouettes 
# Import phylopics

mammals <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/05f87521-20d4-4a05-8ac6-aa0bab7f1394.512.png"))
birds <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/a6ef5684-5683-4a46-aa2d-d39c0d134ba7.512.png"))
amphibians <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/35237d88-e323-44dc-8349-b70fe078c5e7.512.png"))
reptiles <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/a359e147-c088-4c18-a2f1-49abfb2b9325.512.png"))
fish <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/8a9bc7b5-360c-4f2a-9487-723dcd3deb8b.256.png"))

# Add the figures

(figure4 <- ggdraw(ggf4) + 
    draw_image(amphibians, x = 0.14, 
               y = 0.075, 
               width = 0.1, height = 0.05)+
    draw_image(birds, x = 0.33, 
               y = 0.075, 
               width = 0.1, height = 0.07)+
    draw_image(fish, x = 0.53, 
               y = 0.075, 
               width = 0.1, height = 0.04)+
    draw_image(mammals, x = 0.70, 
               y = 0.075, 
               width = 0.1, height = 0.07)+
    draw_image(reptiles, 
               x = 0.92, 
               y = 0.055, 
               width = 0.1, height = 0.08))

# Save it

ggsave("FigureS4.pdf", figure4, 
       width = 14.5, height = 8,
       path = ResultPath)

# Further supplementary figures ################################################

load(paste0(DataPath, "/ResistanceChange.RData"))
load(paste0(DataPath, "/RecoveryChange.RData"))
load(paste0(ResultPath,"/SupMod.RData"))

# Figure S5: Body size effects #################################################

(gs5a <- bs_res %>%
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.value, .width = c(.95, .8, .5)) %>% 
   mutate(stress =gsub("b_", "", .variable),
          stress =gsub("scale", "", stress),
          stress =gsub("bs_g", "Body mass", stress),
          stress =gsub("Body mass:Taxon", "Body mass:", stress)) %>%
   filter(stress!="TaxonReptiles", 
          stress!="TaxonMammals", 
          stress!="TaxonFish",
          stress!="TaxonBirds", 
          stress!="Intercept") %>% 
   ggplot(aes(y = stress, x = .value)) +
   geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
   geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                      interval_size_range = c(0.5, 2),
                      colour="grey40") +
   labs(x="Posterior estimate", y = "Body mass (g) and Taxonomic group") + 
   theme(legend.position = "none"))

(gs5b <- bs_rec %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(stress =gsub("b_", "", .variable),
           stress =gsub("scale", "", stress),
           stress =gsub("bs_g", "Body mass", stress),
           stress =gsub("Body mass:Taxon", "Body mass:", stress)) %>%
    filter(stress!="TaxonReptiles", 
           stress!="TaxonMammals", 
           stress!="TaxonFish",
           stress!="TaxonBirds", 
           stress!="Intercept") %>%  
    ggplot(aes(y = stress, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    labs(x="Posterior estimate", y = "") + 
    theme(legend.position = "none"))


(gs5 <- plot_grid(gs5a, gs5b, 
                  labels = "auto", 
                  nrow = 1,
                  align = "h", 
                  axis = "b") )

ggsave(gs5, filename = "Figure S5.pdf",
       path = ResultPath, 
       width = 10, height = 6)

# Figure S6: Phylogenetic signal ###############################################
# Import data 
load(paste0(ResultPath, "/PSModels.RData"))

# Import phylopics

mammals <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/05f87521-20d4-4a05-8ac6-aa0bab7f1394.512.png"))
birds <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/a6ef5684-5683-4a46-aa2d-d39c0d134ba7.512.png"))
amphibians <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/35237d88-e323-44dc-8349-b70fe078c5e7.512.png"))
reptiles <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/a359e147-c088-4c18-a2f1-49abfb2b9325.512.png"))
fish <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/8a9bc7b5-360c-4f2a-9487-723dcd3deb8b.256.png"))

# Create a common data frame

# Resistance 

psign_res <- rbind(data.frame(taxon= "Amhibians",dist=psign_res_amp$samples$H1), 
                   data.frame(taxon= "Birds", dist=psign_res_bird$samples$H1),
                   data.frame(taxon= "Mammals",dist=psign_res_mam$samples$H1),
                   data.frame(taxon= "Reptiles",dist=psign_res_rep$samples$H1),
                   data.frame(taxon= "Fish",dist=psign_res_fish$samples$H1))

# Calculate the median 

res_median <- psign_res %>%
  group_by(taxon) %>% 
  summarise(median=median(dist))

# Plot 

(ga <- psign_res %>% 
    ggplot(aes(y=taxon, x=dist, color=taxon)) +
    stat_halfeye(aes(color = taxon,
                     fill=after_scale(colorspace::lighten(color, .3))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8,
                 adjust = .5) +
    geom_text(data=res_median,
              aes(x = median, label = format(round(median, 2), nsmall = 2)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    scale_color_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    geom_vline(xintercept = 0.5, linetype = "dashed", colour="grey50") +
    labs(x="Phylogenetic signal", y="")+
    theme(legend.position = "none")+
    xlim(0,1))

# Add the silhouettes 

(ga <- ggdraw(ga) + 
    draw_image(amphibians, x = 0.7, 
               y =0.2, 
               width = 0.1, height = 0.05)+  
    draw_image(birds, x = 0.7, 
               y = 0.35, 
               width = 0.1, height = 0.07)+
    draw_image(fish, x = 0.7, 
               y = 0.52, 
               width = 0.1, height = 0.04)+
    draw_image(mammals, x = 0.70, 
               y = 0.67, 
               width = 0.1, height = 0.07)+
    draw_image(reptiles, 
               x = 0.7, 
               y = 0.8, 
               width = 0.1, height = 0.07))
# Recovery

psign_rec <- rbind(data.frame(taxon= "Amhibians",dist=psign_rec_amp$samples$H1), 
                   data.frame(taxon= "Birds", dist=psign_rec_bird$samples$H1),
                   data.frame(taxon= "Mammals",dist=psign_rec_mam$samples$H1),
                   data.frame(taxon= "Reptiles",dist=psign_rec_rep$samples$H1),
                   data.frame(taxon= "Fish",dist=psign_rec_fish$samples$H1))

# Calculate the median 

rec_median <- psign_rec %>%
  group_by(taxon) %>% 
  summarise(median=median(dist))

# Plot 

(gb <- psign_rec %>% 
    ggplot(aes(y=taxon, x=dist, color=taxon)) +
    stat_halfeye(aes(color = taxon,
                     fill=after_scale(colorspace::lighten(color, .3))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8,
                 adjust = .5) +
    geom_text(data=rec_median,
              aes(x = median, label = format(round(median, 2), nsmall = 2)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    scale_color_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    geom_vline(xintercept = 0.5, linetype = "dashed", colour="grey50") +
    labs(x="Phylogenetic signal", y="")+
    theme(legend.position = "none")+
    xlim(0,1))

# Add the silhouettes 

(gb <- ggdraw(gb) + 
    draw_image(amphibians, x = 0.7, 
               y =0.2, 
               width = 0.1, height = 0.05)+  
    draw_image(birds, x = 0.7, 
               y = 0.35, 
               width = 0.1, height = 0.07)+
    draw_image(fish, x = 0.7, 
               y = 0.52, 
               width = 0.1, height = 0.04)+
    draw_image(mammals, x = 0.70, 
               y = 0.67, 
               width = 0.1, height = 0.07)+
    draw_image(reptiles, 
               x = 0.7, 
               y = 0.8, 
               width = 0.1, height = 0.07))

# Combine figures

(figureS6 <- plot_grid(ga, gb, 
                       labels = "auto",
                       nrow = 1,
                       align = "hv", axis = "b")) 

# Save it

ggsave("FigureS6.pdf", figureS6, 
       width = 10, height = 8,
       path = ResultPath)

# Figure S7: Duration ##########################################################

(gS7a <- duration_res %>%
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.value, .width = c(.95, .8, .5)) %>% 
   mutate(stress =gsub("b_", "", .variable),
          stress =gsub("scale", "", stress)) %>%    
   ggplot(aes(y = stress, x = .value)) +
   geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
   geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                      interval_size_range = c(0.5, 2),
                      colour="grey40") +
   labs(x="Posterior estimate", y = "Duration") + 
   theme(legend.position = "none"))

(gS7b <- duration_rec %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(stress =gsub("b_", "", .variable),
           stress =gsub("scale", "", stress)) %>%    
    ggplot(aes(y = stress, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    labs(x="Posterior estimate", y = "") + 
    theme(legend.position = "none"))

(gS7 <- plot_grid(gS7a, gS7b, 
                  labels = "auto", 
                  nrow = 1,
                  align = "h", 
                  axis = "b") )

ggsave(gS7, filename = "Figure S7.pdf",
       path = ResultPath, 
       width = 10, height = 6)

# Figure S8: Duration vs start year ############################################

# Plot it 

(gS8a <- ggplot(res_data, 
                aes(x=Start,y=Duration, 
                    group=Start)) +
    geom_boxplot(fill="#C4AF8F") +
    scale_x_continuous(breaks = (seq(1950,2015, by=10)))+
    labs(x=""))

(gS8b <- ggplot(rec_data, 
                aes(x=Start,y=Duration, group=Start)) +
    geom_boxplot(fill="#90AFC5") +
    scale_x_continuous(breaks = (seq(1950,2015, by=10)))+
    labs(x="Starting year of the time series"))

(gS8 <- plot_grid(gS8a, gS8b, 
                  labels = "auto", 
                  nrow = 2,
                  align = "v", 
                  axis = "b") )

ggsave(gS8, filename = "Figure S8.pdf",
       path = ResultPath, 
       width = 11, height = 8)

# Figure S9: Decades ###########################################################

(gS9a <- mts_res %>%
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.value, .width = c(.95, .8, .5)) %>% 
   mutate(stress= gsub(":.*", "", .variable),
          stress= gsub("[^0-9.-]", "", stress),
          stress= gsub("[^[:alnum:]]", "", stress),
          Decade = gsub("Decade", "", .variable),
          Decade = gsub(".*:","", Decade),
          Decade =gsub(":", "", Decade),
          Decade =gsub("b_n.threat", "", Decade),
          Decade =gsub("b_", "", Decade)) %>%    
   filter(Decade!="2010") %>% 
   ggplot(aes(x = Decade, y = .value)) +
   geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
   geom_pointinterval(aes(ymin = .lower, ymax = .upper),
                      interval_size_range = c(0.5, 2),
                      colour="grey40") +
   labs(x="Decade", y = "Posterior estimate"))

(gS9b <- mts_rec %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(stress= gsub(":.*", "", .variable),
           stress= gsub("[^0-9.-]", "", stress),
           stress= gsub("[^[:alnum:]]", "", stress),
           Decade = gsub("Decade", "", .variable),
           Decade = gsub(".*:","", Decade),
           Decade =gsub(":", "", Decade),
           Decade =gsub("b_n.threat", "", Decade),
           Decade =gsub("b_", "", Decade)) %>%
    filter(Decade!="2010") %>% 
    ggplot(aes(x = Decade, y = .value)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour="grey50") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    labs(x="Decade", y = "Posterior estimate"))


(gS9 <- plot_grid(gS9a+labs(x=""), gS9b, 
                  labels = "auto", 
                  nrow = 2,
                  align = "v", 
                  axis = "b") )

ggsave(gS9, filename = "Figure S9.pdf",
       path = ResultPath, 
       width = 12, height = 8)

# Figure S10: System taxa #######################################################

# Resistance

# Sample size

sample_size <-res_data %>% 
  group_by(System, Taxon) %>% 
  summarise(n=n())

med <- mst_res %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Taxon = gsub("Taxon", "", .variable),
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_System", "", Taxon),
         Taxon =gsub("Freshwater", "", Taxon),
         Taxon =gsub("Marine", "", Taxon),
         Taxon =gsub("Terrestrial", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon),
         System = gsub("b_System", "", .variable), 
         System = gsub("Taxon", "", System), 
         System = gsub("Mammals", "", System), 
         System = gsub("Birds", "", System), 
         System = gsub("Fish", "", System), 
         System = gsub("Reptiles", "", System), 
         System = gsub("Amphibians", "", System), 
         System =gsub("Taxon", "", System),
         System =gsub(":", "", System),
         System = gsub('[[:digit:]]+', '', System))  %>%   
  group_by(System, Taxon) %>% 
  summarise(med=quantile(.value,probs = 0.9))

sample_size<- sample_size %>% 
  left_join(med)

# Plot intervals 

dat <- mst_res %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Taxon = gsub("Taxon", "", .variable),
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_System", "", Taxon),
         Taxon =gsub("Freshwater", "", Taxon),
         Taxon =gsub("Marine", "", Taxon),
         Taxon =gsub("Terrestrial", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon),
         System = gsub("b_System", "", .variable), 
         System = gsub("Taxon", "", System), 
         System = gsub("Mammals", "", System), 
         System = gsub("Birds", "", System), 
         System = gsub("Fish", "", System), 
         System = gsub("Reptiles", "", System), 
         System = gsub("Amphibians", "", System), 
         System =gsub("Taxon", "", System),
         System =gsub(":", "", System),
         System = gsub('[[:digit:]]+', '', System)) %>% 
  filter(Taxon!="Amphibians"| System!="Marine",
         Taxon!="Fish"|System!="Terrestrial")

(gS10a <- mst_res %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Taxon = gsub("Taxon", "", .variable),
           Taxon =gsub(":", "", Taxon),
           Taxon =gsub("b_System", "", Taxon),
           Taxon =gsub("Freshwater", "", Taxon),
           Taxon =gsub("Marine", "", Taxon),
           Taxon =gsub("Terrestrial", "", Taxon),
           Taxon = gsub('[[:digit:]]+', '', Taxon),
           System = gsub("b_System", "", .variable), 
           System = gsub("Taxon", "", System), 
           System = gsub("Mammals", "", System), 
           System = gsub("Birds", "", System), 
           System = gsub("Fish", "", System), 
           System = gsub("Reptiles", "", System), 
           System = gsub("Amphibians", "", System), 
           System =gsub("Taxon", "", System),
           System =gsub(":", "", System),
           System = gsub('[[:digit:]]+', '', System)) %>%
    filter(Taxon!="Amphibians"|System!="Marine",
           Taxon!="Fish"|System!="Terrestrial")%>% 
    ggplot(aes(y = Taxon, x = .value)) +
    facet_wrap(System~., nrow = 1,scales = "free_x")+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, aes(fill=Taxon), alpha=0.2) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=Taxon, 
                                      label = paste0("n=", n)),
              vjust   = -1.5,
              hjust=-0.1)+
    scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    #xlim(-0.15,0.15)+
    labs(x="Posterior estimate", y = "Taxon") + 
    theme(legend.position = "none"))

# Recovery

# Sample size

sample_size <- rec_data %>% 
  group_by(System, Taxon) %>% 
  summarise(n=n())

med <- mst_rec %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Taxon = gsub("Taxon", "", .variable),
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_System", "", Taxon),
         Taxon =gsub("Freshwater", "", Taxon),
         Taxon =gsub("Marine", "", Taxon),
         Taxon =gsub("Terrestrial", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon),
         System = gsub("b_System", "", .variable), 
         System = gsub("Taxon", "", System), 
         System = gsub("Mammals", "", System), 
         System = gsub("Birds", "", System), 
         System = gsub("Fish", "", System), 
         System = gsub("Reptiles", "", System), 
         System = gsub("Amphibians", "", System), 
         System =gsub("Taxon", "", System),
         System =gsub(":", "", System),
         System = gsub('[[:digit:]]+', '', System))  %>%   
  group_by(System, Taxon) %>% 
  summarise(med=quantile(.value,probs = 0.9))

sample_size<- sample_size %>% 
  left_join(med)

# Plot intervals 

dat <- mst_rec %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Taxon = gsub("Taxon", "", .variable),
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_System", "", Taxon),
         Taxon =gsub("Freshwater", "", Taxon),
         Taxon =gsub("Marine", "", Taxon),
         Taxon =gsub("Terrestrial", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon),
         System = gsub("b_System", "", .variable), 
         System = gsub("Taxon", "", System), 
         System = gsub("Mammals", "", System), 
         System = gsub("Birds", "", System), 
         System = gsub("Fish", "", System), 
         System = gsub("Reptiles", "", System), 
         System = gsub("Amphibians", "", System), 
         System =gsub("Taxon", "", System),
         System =gsub(":", "", System),
         System = gsub('[[:digit:]]+', '', System)) %>% 
  filter(Taxon!="Amphibians"| System!="Marine",
         Taxon!="Fish"|System!="Terrestrial")

(g9b <- mst_rec %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Taxon = gsub("Taxon", "", .variable),
           Taxon =gsub(":", "", Taxon),
           Taxon =gsub("b_System", "", Taxon),
           Taxon =gsub("Freshwater", "", Taxon),
           Taxon =gsub("Marine", "", Taxon),
           Taxon =gsub("Terrestrial", "", Taxon),
           Taxon = gsub('[[:digit:]]+', '', Taxon),
           System = gsub("b_System", "", .variable), 
           System = gsub("Taxon", "", System), 
           System = gsub("Mammals", "", System), 
           System = gsub("Birds", "", System), 
           System = gsub("Fish", "", System), 
           System = gsub("Reptiles", "", System), 
           System = gsub("Amphibians", "", System), 
           System =gsub("Taxon", "", System),
           System =gsub(":", "", System),
           System = gsub('[[:digit:]]+', '', System)) %>%
    filter(Taxon!="Amphibians"|System!="Marine",
           Taxon!="Fish"|System!="Terrestrial")%>% 
    ggplot(aes(y = Taxon, x = .value)) +
    facet_wrap(System~., nrow = 1,scales = "free_x")+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    stat_slab(data=dat, aes(fill=Taxon), alpha=0.2) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=Taxon, 
                                      label = paste0("n=", n)),
              vjust   = -1.5,
              hjust=-0.1)+
    scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    #xlim(-0.15,0.15)+
    labs(x="Posterior estimate", y = "Taxon") + 
    theme(legend.position = "none"))

# Combined figure 

ggf <- plot_grid(g9a+ xlab("")+ 
                   theme(plot.margin =  margin(0.5, 0, 0, 0, "cm"),
                         axis.title.y = element_text(margin = margin(0, -5, 0, 0))),
                 g9b + 
                   theme(plot.margin =  margin(0.5, 0, 0, 0, "cm"),
                         axis.title.y = element_text(margin = margin(0, -5, 0, 0))),
                 nrow = 2, labels = "auto",
                 align = "h", axis = "b") 

# Add our customised titles 

t1 <- ggdraw() + draw_label("Resistance",angle = -90, 
                            fontface='bold', size=18)
t2 <- ggdraw() + draw_label("Recovery",angle = -90, 
                            fontface='bold',size = 18)
titles <- plot_grid(t1,t2, ncol = 1)

(figureS10<- plot_grid(ggf, titles,  
                      ncol=2, rel_widths =c(1, 0.02)))

ggsave("FigureS10.pdf", figureS10, 
       width = 12, height = 10,
       path = ResultPath)

