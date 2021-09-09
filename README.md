# Global patterns of resilience decline in vertebrate populations

Pol Capdevila<sup>1</sup>*, Nicola Noviello<sup>1</sup>, Louise McRae<sup>2</sup>, Robin Freeman<sup>2</sup>, Christopher F Clements<sup>1</sup>

<sup>1</sup>School of Biological Sciences, University of Bristol, 24 Tyndall Ave, BS8 1TQ, Bristol, UK. 

<sup>2</sup>Institute of Zoology, Zoological Society of London, Regent’s Park, London NW1 4RY, UK.

#### Contact: pcapdevila.pc[at]gmail.com

---

## Aim of the study

_Maintaining the resilience of natural populations, their ability to resist and recover from disturbance, is crucial to prevent biodiversity loss. However, the lack of appropriate data and quantitative tools has hampered our understanding of the factors determining resilience on a global scale. Here, we quantified the temporal trends of two key components of resilience – resistance and recovery – in >2,000 population time-series of >1,000 vertebrate species globally. We show that the number of threats to which a population is exposed is the main driver of resilience decline in vertebrate populations. Such declines are driven by a non-uniform loss of different components of resilience (i.e. resistance and recovery). Increased anthropogenic threats accelerating resilience loss through a decline in the recovery ability – but not resistance – of vertebrate populations. These findings suggest we may be underestimating the impacts of global change, highlighting the need to account for the multiple components of resilience in global biodiversity assessments._

---

## Data

- __`RecoveryChange.RData`__: contains the outcomes of the state-space models fit on the recovery time-series. 
- __`ResistanceChange.RData`__: contains the outcomes of the state-space models fit on the recovery time-series.
- __`Body_Mass_means.csv`__: contains the mean values of the adult body weight (g) of the studied species, please check [Noviello et al. 2021](https://doi.org/10.1101/2020.12.17.423192) for further details on how the data was collected.
- __`PopChangeData.RData`__: contains the outcomes of the state-space models from the population abundance time-series from the Living Planet Database.  
- __`amphibians.nex`__: phylogeny of amphibians obtained from [vertlife.org](https://vertlife.org/).  
- __`birds.nex`__: phylogeny of amphibians obtained from [vertlife.org](https://vertlife.org/).  
- __`mammals.nex`__: phylogeny of amphibians obtained from [vertlife.org](https://vertlife.org/).
- __`reptiles.nex`__: phylogeny of amphibians obtained from [vertlife.org](https://vertlife.org/).
- __`fishtree.tre`__: phylogeny of fishes obtained from [Rabosky et al. 2014](https://www.nature.com/articles/s41586-018-0273-1).

---

# Code

To run the statistical analyses we used different R scripts: 

- __`Analyses.R`__: code to analyse the factors influencing resistance and recovery loss. This analyses require a lot of computing power, so be minfult of that.
- __`Figures.R`__: code to create the figures and tables of the study. 
- __`SupAnalyses.R`__: code to do the supplementary analyses. 
- __`SupFigures.R`__: code to create the suplementary figures. 

---

# Software

_R version 4.0.2 or greater_

To download `R`, go to https://www.r-project.org and for `RStudio`, go to https://www.rstudio.com/products/rstudio/download/ .

