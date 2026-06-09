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
                 "Logistic_z",
                 "vanLeeuwen_q",
                 "vanLeeuwen_wq")

# Fit models concurrently? FALSE = sequential (safer for RAM).
# TRUE multiplies memory by length(model_names) since each fit also runs parallel chains.
parallel_models <- FALSE

# Skip re-fitting a model whose .stan file is unchanged since the previous fit.
# Detected by md5sum of the .stan source, stored at tmp/model_output_<name>.hash
# alongside the cached fit at tmp/model_output_<name>.RDS.
# Set FALSE to force a fresh fit regardless of the cache.
reuse_existing_fits <- TRUE

# Skip re-simulating a model whose ODE inputs are unchanged since the previous
# simulation. Detected by a signature of the .stan source, posterior_draws_<name>.RDA,
# simulate.R, A.level, and n_sim_draws, stored at tmp/ODE_kelp_<name>.hash.
# Set FALSE to force fresh simulations regardless of the cache.
reuse_existing_sims <- TRUE

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

source('preference_helpers.R')            # defines preference() / log_switch_point()
source('visualize_stan_model.R')          # defines visualize_model()
source('simulate.R')                      # defines simulate_model()
source('process_simulation.R')            # defines process_model_sim()
source('plot_simulation.R')               # defines plot_model_sim()

## compare_models.R runs at source time and reuses parse_param_block() defined
## by visualize_stan_model.R, so it must be sourced after the definitions above.
source('compare_models.R')

## Parallel mode: each model gets its own multisession worker and
## simulate_model's inner ODE loop runs sequentially inside it
## (internal_cores = 1L) to avoid nested cluster oversubscription. A progressr
## bar ticks once per model per stage.
## Serial mode: run_stage() loops over models with a plain for so each stage's
## console output (cmdstanr per-chain progress, simulate_model's per-kelp-level
## prints) prints live. simulate_model uses all n_cores for its inner ODE loop.
if (parallel_models) {
  mc_n <- min(length(model_names), n_cores)
  cat("Running visualize / simulate / process / plot stages for",
      length(model_names), "models in parallel (",
      paste(model_names, collapse = ", "), ") on", mc_n, "core(s).\n")
  plan(multisession, workers = mc_n)
  on.exit(plan(sequential), add = TRUE)
  handlers(global = TRUE)
  handlers("progress")
  internal_cores <- 1L
} else {
  internal_cores <- n_cores
}

run_stage <- function(label, fn, ...) {
  if (parallel_models) {
    with_progress({
      p <- progressor(steps = length(model_names))
      future_lapply(model_names, function(m) {
        source("load_packages.R")
        options(device = function(...) pdf(NULL))
        out <- fn(m, ...)
        p(sprintf("[%s %s] done", label, m))
        out
      }, future.globals = TRUE, future.seed = TRUE, future.stdout = NA)
    })
  } else {
    out <- vector("list", length(model_names))
    for (i in seq_along(model_names)) {
      m <- model_names[i]
      cat(sprintf("\n[%s %d/%d] %s...\n",
                  label, i, length(model_names), m))
      out[[i]] <- fn(m, ...)
    }
    out
  }
}

invisible(run_stage("visualize", visualize_model))
invisible(run_stage("simulate",  simulate_model,
                    n_draws = n_sim_draws, internal_cores = internal_cores))
invisible(run_stage("process",   process_model_sim))
invisible(run_stage("plot",      plot_model_sim))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~