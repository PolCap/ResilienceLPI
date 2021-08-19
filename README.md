# Global patterns of resilience decline in vertebrate populations

This repository contains the code for the statistical analyses of the manuscript "Global patterns of resilience decline in vertebrate populations".

Pol Capdevila<sup>1</sup>*, Nicola Noviello<sup>1</sup>, Louise McRae<sup>2</sup>, Robin Freeman<sup>2</sup>, Christopher F Clements<sup>1</sup>

<sup>1</sup>School of Biological Sciences, University of Bristol, 24 Tyndall Ave, BS8 1TQ, Bristol, UK. 

<sup>2</sup>Institute of Zoology, Zoological Society of London, Regent’s Park, London NW1 4RY, UK.

#### Contact: pcapdevila.pc[at]gmail.com

# Aim
_In this work we quantify the temporal trends of two key components of resilience – resistance and recovery – in >2,000 population time-series of >1,000 vertebrate species globally._

# Data

- __`LambdaChangeData.RData`__: contains the outcomes of the state-space models fit on the recovery time-series. 
- __`NegativeLambdaChangeData.RData`__: contains the outcomes of the state-space models fit on the recovery time-series.
- __`Body_Mass_means.csv`__: contains the mean values of the adult body weight (g) of the studied species according to different sources.


# Code

To run the statistical analyses we used different R scripts: 

- __`Analyses.R`__: code to analyse the factors influencing resistance and recovery loss. This analyses require a lot of computing power, so be minfult of that.
- __`Figures.R`__: code to create the figures and tables of the study. 
- __`SupAnalyses.R`__: code to do the supplementary analyses. 
- __`SupFigures.R`__: code to create the suplementary figures. 

_For a detailed explanation of the methods, please see the methods section of the main manuscript._

# Software
_R version 3.6.1 or greater_

To download `R`, go to https://www.r-project.org and for `RStudio`, go to https://www.rstudio.com/products/rstudio/download/ .

