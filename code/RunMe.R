## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
closeAllConnections()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('load_packages.R')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Relative file paths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
code    <- "../code"
models  <- "../code/models"
data    <- "../data"
results <- "../results"
tmp     <- "../results/tmp"
tables  <- "../tables"
figs    <- "../figs"

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Choices ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Models to fit and compare. Each entry must have a 'stan_model_<name>.stan'
# file in code/models/.
# Set to a single-element vector to fit just one model.
model_names <- c("Logistic", 
                 "LogisticClearance",
                 "vanLeeuwen_q",
                 "vanLeeuwen_wq")

# Fit models concurrently? FALSE = sequential (safer for RAM).
# TRUE multiplies memory by length(model_names) since each fit also runs 4 parallel chains.
parallel_models <- TRUE

# Skip re-fitting a model whose .stan file is unchanged since the previous fit.
# Detected by md5sum of the .stan source, stored alongside model_output_<name>.RDS.
# Set FALSE to force a fresh fit regardless of the cache.
reuse_existing_fits <- TRUE

# Specify number of MCMC iterations
warmup_iter <- 200
sampling_iter <- 500

# Specify number of cores to use for parallel computing
n_cores <- max(1L, detectCores() - 1L)

# Number of posterior draws to use in simulate_model(). NULL = keep all draws.
# Lower this (e.g. 200) to speed up the ODE simulation at the cost of CI estimation.
n_sim_draws <- 200

# Specify kelp abundances (grams) at which to simulate
A.level = c("low" = 50,
            "high" = 250)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Source files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source('plot_juvenile_kelp.R')
source('empirical.R')
source('format_data.R')

source('fit_all_models.R')
source('compare_models.R')

source('visualize_stan_model.R')          # defines visualize_model()
source('simulate.R')                      # defines simulate_model()
source('process_simulation.R')            # defines process_model_sim()
source('plot_simulation.R')               # defines plot_model_sim()

if (parallel_models) {
  mc_n <- min(length(model_names), n_cores)
  cat("Running visualize / simulate / process / plot stages for",
      length(model_names), "models in parallel (",
      paste(model_names, collapse = ", "), ") on", mc_n, "core(s).\n")
  invisible(parallel::mclapply(model_names, visualize_model,   mc.cores = mc_n))
  invisible(parallel::mclapply(model_names, simulate_model,
                               n_draws = n_sim_draws, internal_cores = 1L,
                               mc.cores = mc_n))
  invisible(parallel::mclapply(model_names, process_model_sim, mc.cores = mc_n))
  invisible(parallel::mclapply(model_names, plot_model_sim,    mc.cores = mc_n))
} else {
  invisible(lapply(model_names, visualize_model))
  invisible(lapply(model_names, simulate_model, n_draws = n_sim_draws))
  invisible(lapply(model_names, process_model_sim))
  invisible(lapply(model_names, plot_model_sim))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~