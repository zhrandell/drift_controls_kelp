
## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
rm(list = ls())

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Choose the model [1] or [2]
sel.model <- c('Logistic', 'vanLeeuwen')[1]
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Libraries for 'format_data.R'
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

## Libraries for 'fit_stan_model.R' & 'analyze_visualize_stan_model.R'
library(cmdstanr)
library(rstan)
library(StanHeaders)
library(shinystan)
library(posterior)
library(bayesplot)

## Libraries for 'analyze_visualize_stan_model.R'
library(egg)
library(Cairo) # requires installation of Xcode

## Libraries for 'simulate.R' and 'plot_simulation'
library(deSolve)
library(reshape2)
library(ggplot2)
library(egg)
library(tidyselect)
library(scales)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## check wd is appropriate
getwd()

## relative file paths
code <- "../code"
data <- "../data"
results <- "../results"
figs <- "../figs"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source('format_data.R')
source('fit_stan_model.R')
source('analyze_visualize.R')
source('simulate.R')
source('simulate_process.R')
source('plot_simulation.R')


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~