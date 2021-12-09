##############################
# REQUIRE PACKAGES WITH PACMAN
##############################

remotes::install_github("gavinsimpson/ggvegan", force = T)

if (!require("OpenStreetMap")) install.packages("OpenStreetMap")
library(OpenStreetMap)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,reshape2,zoo,vegan,readxl,lme4,
               viridis,patchwork,ggpubr,visdat,labdsv,
               reshape2,ggvegan,janitor,MuMIn,
               sf,tigris,ggrepel,osmdata,rgdal,ggmap)

source("./R/functions/ScreePlotFunction.R")

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}
normalize <- function(x){(x-min(x))/(max(x)-min(x))}
standard_error <- function(x) {sd(x) / sqrt(length(x))}
percent_change <- function(x) {((x - lead(x))/(x))*100}
zscore <- function(x) {(x - mean(x))/sd(x)}
slope <- function(x, y){
  mean_x <- mean(x)
  mean_y <- mean(y)
  nom <- sum((x - mean_x)*(y-mean_y))
  denom <- sum((x - mean_x)^2)
  m <- nom / denom
  return(m)
}

