##################
# REQUIRE PACKAGES
##################

set.seed(1)
if (!require("ggvegan")) remotes::install_github("gavinsimpson/ggvegan", force = T)
if (!require("OpenStreetMap")) install.packages("OpenStreetMap")
if (!require("pacman")) install.packages("pacman")


pacman::p_load(tidyverse,zoo,vegan,readxl,lme4,OpenStreetMap,
               viridis,patchwork,ggpubr,visdat,labdsv,psych,
               reshape2,ggvegan,janitor,MuMIn,kableExtra,
               sf,tigris,ggrepel,osmdata,rgdal,ggmap,corrplot,
               MCMCvis, lubridate,tidybayes,ncdf4,imputeTS,devtools,
               scales,rjags,coda,R2jags,gridExtra,PerformanceAnalytics,
               broom.mixed)

rmse = function(m, o){
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

