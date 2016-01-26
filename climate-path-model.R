#' # Path Modeling of Climatic Influence on Spillover
#' I'm using the `piecewiseSEM package, which operationalized the methods from Shipley 2009 & 2013 in particular for fitting and assessing path models that enable applying hierarchical modeling methods.
#' 
#+ load data and packages
load("pathmodel_data_20160112.RData")
library(plyr); library(dplyr); library(piecewiseSEM)
library(lmerTest); library(nlme)