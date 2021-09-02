# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figures.R         
# - DATE:        15/09/2020
# - DESCRIPTION: Use models from PosLambdaAnalyses.R to produce figures
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(ggplot2, lib.loc = "C:/Users/fg20342/Documents/R/") #Note that I use an old version
# I use an older version given a bug with the use of the function coord_proj
library(dplyr)
library(cowplot)
library(broom)
library(mapdata)
library(maps)
library(RColorBrewer)
library(tidyr)
library(ggalt)
library(proj4)
library(sf)       
library(lwgeom)   
library(ggpubr)

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

path <- gsub("LPI/R", "LPI/", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"R")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(DataPath, "/RecoveryChange.RData"))
load(paste0(DataPath, "/ResistanceChange.RData"))

# Figure 1: Maps ###############################################################

# Set the world map
 
world <- map_data("world")

# Create plot1

(p1 <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray70", fill = "gray70", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() + 
    geom_point(data = res_data, 
               aes(x = Longitude, y = Latitude, 
                   colour = mu, size=n.threat),
               alpha=0.6) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(size=guide_legend(title = "Number of stressors",
                             title.position = "top"))+
    scale_color_gradient2("Resistance trend",
                          midpoint = 0, 
                          low = "#c1666b", 
                          mid = "#e4dfda",
                          high = "#4281a4",
                          guide =guide_colourbar(nbin=100,
                                                 barwidth = 9, 
                                                 title.position = "top"),
                          breaks=c(-0.35, 0.0, 0.35, 0.75))+
    expand_limits(x = 0, y = 0)+
    theme(legend.position = "none",
          legend.direction = "horizontal",
          legend.title = element_text(size = 12, vjust= -1),
          legend.text = element_text(size = 10),
          plot.margin = unit(c(0,0,1.5,6), units = , "cm")))

# Create distribution 1

gd1 <- ggplot(res_data, aes(x = mu)) +
  geom_density(adjust = 3, alpha = .6,
               position = "stack", 
               fill = "grey40", color = "grey40") +
               #fill = "#C4AF8F", color = "#C4AF8F") + 
  geom_vline(aes(xintercept = 0),
             linetype = "longdash", size=0.5) +
  labs(x="Trend", y="")

# Create legend1

legend1 <- ggplot() +
  geom_map(map = world, data = world,
           aes(long, lat, map_id = region), 
           color = "gray70", fill = "gray70", size = 0.3) +
  coord_proj("+proj=wintri") +
  theme_map() + 
  geom_point(data = res_data, 
             aes(x = Longitude, y = Latitude, 
                 colour = mu),
             alpha=0.6) +
  scale_y_continuous(limits = c(-80, 80)) +
  scale_color_gradient2("Trend",
                        midpoint = 0, 
                        low = "#c1666b", 
                        mid = "#e4dfda",
                        high = "#4281a4",
                        guide =guide_colourbar(nbin=100,
                                               barwidth = 9),
                        breaks=c(-0.35, 0.0, 0.35, 0.75))+
  expand_limits(x = 0, y = 0)+
  theme(legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.margin = unit(c(0,0,0,6), units = , "cm"))

# Create legend2

legend2 <- ggplot() +
  geom_map(map = world, data = world,
           aes(long, lat, map_id = region), 
           color = "gray70", fill = "gray70", size = 0.3) +
  coord_proj("+proj=wintri") +
  theme_map() + 
  geom_point(data = res_data, 
             aes(x = Longitude, y = Latitude, 
                 size=n.threat),
             alpha=0.6) +
  scale_y_continuous(limits = c(-80, 80)) +
  guides(size=guide_legend(title = "Number of stressors",
                           title.position = "top"))+
  expand_limits(x = 0, y = 0)+
  theme(legend.direction = "horizontal",
        legend.title = element_text(size = 12, vjust= -1),
        legend.text = element_text(size = 10),
        plot.margin = unit(c(0,0,0,6), units = , "cm"))

# Extract "legends only" from ggplot object

legend1 <- get_legend(legend1)
legend2 <- get_legend(legend2)

# Plot 1 ----

gg1 <- ggdraw(p1) +
  draw_plot(gd1+theme(axis.title.x = element_text(size = 12,
                                                  margin = margin(0, 0, 0, 0))),
            x = -0.3, y = -0.2, scale = 0.45)+
  draw_plot(legend1, 
            x=0.4, y=-0.43, scale = 0.8)+ 
  draw_plot(legend2, x=0.65, y=-0.42, scale = 0.8)

# Create plot2 ---------------

(p2 <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray70", fill = "gray70", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() + 
    geom_point(data = rec_data, 
               aes(x = Longitude, y = Latitude, 
                   colour = mu, size=n.threat),
               alpha=0.6) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(size=guide_legend(title = "Number of stressors",
                             title.position = "top"))+
    scale_color_gradient2("Resistance trend",
                          midpoint = 0, 
                          low = "#c1666b", 
                          mid = "#e4dfda",
                          high = "#4281a4",
                          guide =guide_colourbar(nbin=100,
                                                 barwidth = 9, 
                                                 title.position = "top"),
                          breaks=c(-0.35, 0.0, 0.35, 0.75))+
    expand_limits(x = 0, y = 0)+
    theme(legend.position = "none",
          legend.direction = "horizontal",
          legend.title = element_text(size = 12, vjust= -1),
          legend.text = element_text(size = 10),
          plot.margin = unit(c(0,0,1.5,6), units = , "cm")))

# Create distribution 2

gd2 <- ggplot(rec_data, aes(x = mu)) +
  geom_density(adjust = 3, alpha = .6,
               position = "stack", 
               fill="grey40", color = "grey40") +
  geom_vline(aes(xintercept = 0),
             linetype = "longdash", size=0.5) +
  labs(x="Trend", y="")

# Create legend3

legend3 <- ggplot() +
  geom_map(map = world, data = world,
           aes(long, lat, map_id = region), 
           color = "gray70", fill = "gray70", size = 0.3) +
  coord_proj("+proj=wintri") +
  theme_map() + 
  geom_point(data = rec_data, 
             aes(x = Longitude, y = Latitude, 
                 colour = mu),
             alpha=0.6) +
  scale_y_continuous(limits = c(-80, 80)) +
  scale_color_gradient2("Trend", 
                        midpoint = 0, 
                        low = "#c1666b", 
                        mid = "#e4dfda",
                        high = "#4281a4",
                        guide =guide_colourbar(nbin=100,
                                               barwidth = 9),
                        breaks=c(-0.35, 0.0, 0.35, 0.75))+
  expand_limits(x = 0, y = 0)+
  theme(legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.margin = unit(c(0,0,0,6), units = , "cm"))

# Create legend4

legend4 <- ggplot() +
  geom_map(map = world, data = world,
           aes(long, lat, map_id = region), 
           color = "gray70", fill = "gray70", size = 0.3) +
  coord_proj("+proj=wintri") +
  theme_map() + 
  geom_point(data = rec_data, 
             aes(x = Longitude, y = Latitude, 
                 size=n.threat),
             alpha=0.6) +
  scale_y_continuous(limits = c(-80, 80)) +
  guides(size=guide_legend(title = "Number of stressors",
                           title.position = "top"))+
  expand_limits(x = 0, y = 0)+
  theme(legend.direction = "horizontal",
        legend.title = element_text(size = 12, vjust= -1),
        legend.text = element_text(size = 10),
        plot.margin = unit(c(0,0,0,6), units = , "cm"))

# Extract "legends only" from ggplot object

legend3 <- get_legend(legend3)
legend4 <- get_legend(legend4)

# Plot2 ----

(gg2 <- ggdraw(p2) +
    draw_plot(gd2+theme(axis.title.x = element_text(size = 12,
                        margin = margin(0, 0, 0, 0))),
              x = -0.3, y = -0.2, scale = 0.45)+
    draw_plot(legend3, 
              x=0.4, y=-0.43, scale = 0.8)+ 
    draw_plot(legend4, x=0.65, y=-0.42, scale = 0.8))
    

# Put the maps together ------------

figure1 <- ggarrange(gg1, 
                 gg2,
                 labels = "auto",
                 nrow = 2,
                 align = "v")

# Save it

ggsave("Figure1.pdf", figure1, 
       width = 10, height = 10,
       path = ResultPath)

# Section 2 ####################################################################
# Load packages 
detach("package:ggpubr", unload=TRUE)
detach("package:ggalt", unload=TRUE)
detach("package:cowplot", unload=TRUE)
detach("package:ggplot2", unload=TRUE)

library(tidybayes)
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
library(performance)

# Load model data 

load(paste0(DataPath, "/Factors.RData"))

# Set the theme for the plots

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=15, margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                  axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color="black", size = 12),
                  axis.text.y = element_text(color="black", size = 12),
                  strip.text = element_text(size = 12),
                  axis.ticks = element_line(color="black"),
                  plot.margin = unit(c(0,0,0,0), units = , "cm")))


# Figure 3: General model outcomes #############################################

# Panel a: negative lambda vs number of stressors ------------------------------

# Sample size

sample_size <- res_data %>% 
  group_by(n.threat) %>% 
  summarise(n=n()) %>% 
  mutate(n.threat=factor(n.threat))

med <- m4 %>%
  gather_draws(b_n.threat0, b_n.threat1, b_n.threat2, b_n.threat3) %>%
  mutate(.variable = gsub("b_n.threat", "", .variable)) %>% 
  group_by(.variable) %>% 
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med, by=c("n.threat"=".variable"))

# Plot intervals 

dat <- m4 %>%
  gather_draws(b_n.threat0, b_n.threat1, b_n.threat2, b_n.threat3) %>%
  mutate(.variable = gsub("b_n.threat", "", .variable))

(g1 <- m4 %>%
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

# Compute the R2 of the model 

r2_bayes(m4, robust = T, ci=0.95)

# Panel d: positive lambda vs number of stressors --------------------------------

# Sample size

sample_size <- rec_data %>% 
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

# Compute the R2 of the model 

r2_bayes(m1, robust = T, ci=0.95)

# Panel b: System effects on negative lambda ----------------------------------- 
# Sample size

sample_size <- res_data %>% 
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

# Compute the R2 of the model 

r2_bayes(m5, robust = T, ci=0.95)

# Panel e: System effects on positive lambda ----------------------------------- 

# Sample size

sample_size <- rec_data %>% 
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

# Compute the R2 of the model 

r2_bayes(m2, robust = T, ci=0.95)

# Panel c: Taxon effects on negative lambda ----------------------------------- 

# Sample size

sample_size <- res_data %>% 
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

# Compute the R2 of the model 

r2_bayes(m6, robust = T, ci=0.95)

# Panel f: Taxon effects on positive lambda ----------------------------------- 

# Sample size

sample_size <- rec_data %>% 
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

# Compute the R2 of the model 

r2_bayes(m3, robust = T, ci=0.95)


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

(figure2<- plot_grid(ggf, titles,  
                     ncol=2, rel_widths =c(1, 0.02)))

ggsave("Figure3.pdf", figure2, 
       width = 12, height = 10,
       path = ResultPath)

# Figure 4: Interactions System ################################################

# Resistance -------------------------------------------------------------------

# Sample size

sample_size <- res_data %>% 
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

# Plot

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

# Compute the R2 of the model 

r2_bayes(mi3, robust = T, ci=0.95)

# Recovery ---------------------------------------------------------------------

# Sample size

sample_size <- rec_data %>% 
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

# Plot

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

# Compute the R2 of the model 

r2_bayes(mi1, robust = T, ci=0.95)

# Combined Figure --------------------------------------------------------------

ggf <- plot_grid(g1+xlab(""), g2 + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()), 
                 labels = "auto",
                 nrow = 2,
                 align = "h", axis = "b") 

# Add titles 

t1 <- ggdraw() + draw_label("Resistance",angle = -90, 
                            fontface='bold', size=18)
t2 <- ggdraw() + draw_label("Recovery",angle = -90, 
                            fontface='bold',size = 18)
titles <- plot_grid(t1,t2, ncol = 1)

# Combine again

(figure3<- plot_grid(ggf, titles,  
                     ncol=2, rel_widths =c(1, 0.02)))

# Save the plot

ggsave("Figure4.pdf", figure3, 
       width = 12, height = 10,
       path = ResultPath)

# Figure 5: Interaction Taxon ##################################################

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

# Resistance -------------------------------------------------------------------

# Sample size

sample_size <- res_data %>% 
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

# Plot 

dat <- mi4 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         Taxon = gsub("Taxon", "", .variable), 
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_n.threat", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon)) 

(g1 <- mi4 %>%
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
    labs(x="Posterior estimate", y = "Number of stressors") + 
    theme(legend.position = "none"))

# Compute the R2 of the model 

r2_bayes(mi4, robust = T, ci=0.95)

# Recovery ---------------------------------------------------------------------

# Sample size

sample_size <- rec_data %>% 
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

# Plot it

dat <- mi2 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(n.threat= gsub("[^0-9.-]", "", .variable),
         n.threat= gsub("[^[:alnum:]]", "", n.threat),
         Taxon = gsub("Taxon", "", .variable), 
         Taxon =gsub(":", "", Taxon),
         Taxon =gsub("b_n.threat", "", Taxon),
         Taxon = gsub('[[:digit:]]+', '', Taxon)) 

(g2 <- mi2 %>%
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
    labs(x="Posterior estimate", y = "Number of stressors") + 
    theme(legend.position = "none"))

# Compute the R2 of the model 

r2_bayes(mi3, robust = T, ci=0.95)

# Combined Figure --------------------------------------------------------------

ggf <- plot_grid(g1+xlab(""), g2 + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()), 
  labels = "auto",
  nrow = 2,
  align = "h", axis = "b") 

# Add titles 

t1 <- ggdraw() + draw_label("Resistance",angle = -90, 
                            fontface='bold', size=18)
t2 <- ggdraw() + draw_label("Recovery",angle = -90, 
                            fontface='bold',size = 18)
titles <- plot_grid(t1,t2, ncol = 1)

# Combine again 

(ggf4<- plot_grid(ggf, titles,  
                     ncol=2, rel_widths =c(1, 0.02)))

# Import phylopics

mammals <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/05f87521-20d4-4a05-8ac6-aa0bab7f1394.512.png"))
birds <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/a6ef5684-5683-4a46-aa2d-d39c0d134ba7.512.png"))
amphibians <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/35237d88-e323-44dc-8349-b70fe078c5e7.512.png"))
reptiles <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/a359e147-c088-4c18-a2f1-49abfb2b9325.512.png"))
fish <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/8a9bc7b5-360c-4f2a-9487-723dcd3deb8b.256.png"))

# Add the silhouettes 

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

ggsave("Figure5.pdf", figure4, 
       width = 14.5, height = 8,
       path = ResultPath)

