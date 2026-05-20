## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit every model in the registry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Reads from the caller's environment (set in RunMe.R):
##   code, results, warmup_iter, sampling_iter, n_cores
##   model_names      character vector of models to fit
##   parallel_models  TRUE = fit models concurrently via mclapply (memory-heavy);
##                    FALSE = sequential. Each model fit already uses 4 parallel
##                    chains internally, so concurrent models multiply RAM/CPU
##                    usage by length(model_names).

source('fit_stan_model.R')             # defines fit_model()

if (parallel_models) {
  ## Concurrent fits: give each fit a single core for chain parallelism so total
  ## load = length(model_names) * 4 chains. Adjust if memory is tight.
  n_per_fit <- max(1L, floor(n_cores / length(model_names)))
  fits <- parallel::mclapply(model_names,
                             FUN = function(m) fit_model(m, parallel = n_per_fit),
                             mc.cores = length(model_names))
} else {
  fits <- lapply(model_names, fit_model)
}

names(fits) <- model_names
cat("Fitted models:", paste(model_names, collapse = ", "), "\n")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
