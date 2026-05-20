## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit system of ODEs to empirical data in STAN ~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Defines `fit_model(name)` which compiles 'models/stan_model_<name>.stan',
## samples using the shared `loss_dat` object, and writes the resulting
## CmdStanMCMC fit to results/model_output_<name>.RDS.
## Callers (fit_all_models.R) supply `models`, `results`, `warmup_iter`,
## `sampling_iter`, and `n_cores` from the global environment.
##
## Skip-if-unchanged: when `reuse_existing` is TRUE (default driven by the
## `reuse_existing_fits` global), the function md5sums the .stan source and
## compares against results/model_output_<name>.hash. If the hash matches the
## previous fit's, it loads and returns the saved RDS instead of resampling.

fit_model <- function(name,
                      data_path     = paste0(results, "/loss_dat.Rdata"),
                      warmup        = warmup_iter,
                      sampling      = sampling_iter,
                      chains        = 4L,
                      parallel      = n_cores,
                      adapt_delta   = 0.80,
                      force_recompile = FALSE,
                      reuse_existing  = if (exists("reuse_existing_fits", inherits = TRUE))
                                          reuse_existing_fits else FALSE) {

  stan_file <- paste0(models, "/stan_model_", name, ".stan")
  if (!file.exists(stan_file)) {
    stop(sprintf("Stan file not found: %s", stan_file))
  }

  out_file     <- paste0(results, "/model_output_", name, ".RDS")
  hash_file    <- paste0(results, "/model_output_", name, ".hash")
  current_hash <- unname(tools::md5sum(stan_file))

  ## Reuse cached fit if .stan unchanged and a prior fit exists
  if (reuse_existing && file.exists(out_file) && file.exists(hash_file)) {
    saved_hash <- readLines(hash_file, warn = FALSE)
    if (length(saved_hash) > 0 && identical(saved_hash[1], current_hash)) {
      message("[", name, "] reusing cached fit (.stan unchanged).")
      return(invisible(readRDS(out_file)))
    }
  }

  loss_dat <- readRDS(data_path)
  model    <- cmdstan_model(stan_file, force_recompile = force_recompile)

  fit <- model$sample(data            = loss_dat,
                      chains          = chains,
                      iter_warmup     = warmup,
                      iter_sampling   = sampling,
                      adapt_delta     = adapt_delta,
                      parallel_chains = parallel)

  fit$save_object(file = out_file)
  writeLines(current_hash, hash_file)
  invisible(fit)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
