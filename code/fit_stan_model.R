## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit system of ODEs to empirical data in STAN ~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Defines `fit_model(name)` which compiles 'stan_model_<name>.stan', samples
## using the shared `loss_dat` object, and writes the resulting CmdStanMCMC fit
## to results/model_output_<name>.RDS.
## Callers (fit_all_models.R) supply `code`, `results`, `warmup_iter`,
## `sampling_iter`, and `n_cores` from the global environment.

fit_model <- function(name,
                      data_path     = paste0(results, "/loss_dat.Rdata"),
                      warmup        = warmup_iter,
                      sampling      = sampling_iter,
                      chains        = 4L,
                      parallel      = n_cores,
                      adapt_delta   = 0.80,
                      force_recompile = FALSE) {

  loss_dat <- readRDS(data_path)

  stan_file <- paste0(code, "/stan_model_", name, ".stan")
  if (!file.exists(stan_file)) {
    stop(sprintf("Stan file not found: %s", stan_file))
  }

  model <- cmdstan_model(stan_file, force_recompile = force_recompile)

  fit <- model$sample(data            = loss_dat,
                      chains          = chains,
                      iter_warmup     = warmup,
                      iter_sampling   = sampling,
                      adapt_delta     = adapt_delta,
                      parallel_chains = parallel)

  out_file <- paste0(results, "/model_output_", name, ".RDS")
  fit$save_object(file = out_file)
  invisible(fit)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
