## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit every model in the registry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Reads from the caller's environment (set in RunMe.R):
##   code, results, warmup_iter, sampling_iter, n_cores
##   model_names      character vector of models to fit
##   parallel_models  TRUE = fit models concurrently via a PSOCK cluster
##                    (cross-platform: macOS, Windows, Linux; memory-heavy);
##                    FALSE = sequential. Each model fit already uses 4 parallel
##                    chains internally, so concurrent models multiply RAM/CPU
##                    usage by length(model_names).

source('fit_stan_model.R')             # defines fit_model()

## Validate that every name in `model_names` has a matching .stan file.
stan_files <- paste0(models, "/stan_model_", model_names, ".stan")
missing    <- model_names[!file.exists(stan_files)]
if (length(missing) > 0) {
  stop("Missing Stan model file(s) in ", normalizePath(models, mustWork = FALSE), ":\n  ",
       paste(paste0("stan_model_", missing, ".stan"), collapse = ", "),
       "\nCheck spelling of model_names in RunMe.R or add the .stan file(s).",
       call. = FALSE)
}

if (parallel_models) {
  ## Concurrent fits: give each fit a single core for chain parallelism so total
  ## load = length(model_names) * 4 chains. Adjust if memory is tight.
  n_per_fit <- max(1L, floor(n_cores / length(model_names)))
  cat("Fitting", length(model_names), "models in parallel (",
      paste(model_names, collapse = ", "), ") with",
      n_per_fit, "chain core(s) per model.\n")

  ## PSOCK cluster (cross-platform: macOS, Windows, Linux). Workers start with
  ## empty environments, so paths, fit_model(), and the iteration settings are
  ## exported and load_packages.R is sourced on each one. First parallel block
  ## in the pipeline pays a one-time startup cost (~seconds) for worker spawn
  ## and package loading; subsequent stages in RunMe.R rebuild their own
  ## cluster but the package install/compile cache is warm by then.
  cl <- parallel::makeCluster(length(model_names))
  wd <- getwd()
  parallel::clusterExport(
    cl,
    varlist = c("fit_model", "models", "results", "tmp",
                "warmup_iter", "sampling_iter", "n_cores",
                "reuse_existing_fits", "n_per_fit"),
    envir = .GlobalEnv)
  parallel::clusterCall(cl, function(w) {
    setwd(w); source("load_packages.R")
  }, wd)

  fits <- tryCatch(
    parallel::parLapply(cl, model_names,
                        fun = function(m) fit_model(m, parallel = n_per_fit)),
    finally = parallel::stopCluster(cl))
} else {
  fits <- lapply(model_names, fit_model)
}

names(fits) <- model_names
cat("Fitted models:", paste(model_names, collapse = ", "), "\n")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
