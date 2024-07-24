
## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
rm(list = ls())


## read in libraries for 'format_data.R'
library(tidyverse)
library(stats)
library(gridExtra)
library(gtable)
library(grid)
library(deSolve)
library(rstudioapi)
library(reshape2)
library(hexbin)
library(janitor)
library(diffdf)

## read in libraries for 'fit_stan_model.R' & 'analyze_visualize_stan_model.R'
library(cmdstanr)
library(rstan)
library(StanHeaders)
library(shinystan)
library(posterior)
library(bayesplot)

## read in libraries for 'analyze_visualize_stan_model.R'
library(egg)



## check wd is appropriate
getwd()


## hardcode relative file paths
code <- "../code"
data_input <- "../data_input"
data_output <- "../data_output"
figs <- "../figs"


## end start ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##


source('format_data.R')
source('fit_stan_model.R')
source('analyze_visualize_stan_model.R')
