## Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

## Libraries for parallelism + progress reporting (used by fit_all_models.R,
## RunMe.R, simulate.R)
load_pkg("future")
load_pkg("future.apply")
load_pkg("progressr")

## Libraries for 'empirical.R'
load_pkg("tidyverse")
load_pkg("scales")
load_pkg("patchwork")

## Libraries for 'format_data.R'
load_pkg("stats")
load_pkg("gridExtra")
load_pkg("gtable")
load_pkg("grid")
load_pkg("deSolve")
load_pkg("rstudioapi")
load_pkg("reshape2")
load_pkg("hexbin")
load_pkg("janitor")
load_pkg("diffdf")

## Libraries for 'fit_all_models.R' & 'analyze_visualize.R'
load_pkg("cmdstanr")
load_pkg("rstan")
load_pkg("StanHeaders")
load_pkg("shinystan")
load_pkg("posterior")
load_pkg("bayesplot")

## Libraries for 'compare_models.R'
load_pkg("loo")
load_pkg("bayesplot")

## Libraries for 'analyze_visualize.R'
load_pkg("patchwork")
load_pkg("Cairo") # requires installation of XQuartz on machine

## Libraries for 'simulate.R' and 'plot_simulation'
load_pkg("deSolve")
load_pkg("reshape2")
load_pkg("ggplot2")
load_pkg("tidyselect")
load_pkg("scales")
load_pkg("pbapply")
load_pkg("parallel")
load_pkg("grid")
