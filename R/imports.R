
# load packages

library(spatstat)
library(ggplot2)
library(reshape2) # for melt() in my function plot_sample()
library(dplyr)
library(readr)
library(ape) # for Moran.I() in my function var_estim_vw()
library(spsurvey)
library(spdep) # for geary() in my function var_estim_vstr2()
library(sp) # for SpatialPoints() in my function var_estim_vstr2()
