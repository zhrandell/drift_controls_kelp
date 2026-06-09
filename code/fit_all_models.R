## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit every model in the registry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Reads from the caller's environment (set in RunMe.R):
##   code, results, tmp, models, warmup_iter, sampling_iter, n_cores
##   model_names      character vector of models to fit
##   parallel_models  TRUE = fit models concurrently via future::multisession
##                    (cross-platform: macOS, Windows, Linux; memory-heavy);
##                    FALSE = sequential. Each model fit already uses 4 parallel
##                    chains internally, so concurrent models multiply RAM/CPU
##                    usage by length(model_names).

## fit_model(name): compile 'models/stan_model_<name>.stan', sample using the
## shared `loss_dat` object, and write the CmdStanMCMC fit to
## tmp/model_output_<name>.RDS.
##
## Skip-if-unchanged: when `reuse_existing` is TRUE (default driven by the
## `reuse_existing_fits` global), the function md5sums the .stan source and
## compares against tmp/model_output_<name>.hash. If the hash matches the
## previous fit's, it loads and returns the saved RDS instead of resampling.
fit_model <- function(name,
                      data_path     = paste0(tmp, "/loss_dat.RData"),
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

  out_file     <- paste0(tmp, "/model_output_", name, ".RDS")
  hash_file    <- paste0(tmp, "/model_output_", name, ".hash")
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

## Validate that every name in `model_names` has a matching .stan file.
stan_files <- paste0(models, "/stan_model_", model_names, ".stan")
missing    <- model_names[!file.exists(stan_files)]
if (length(missing) > 0) {
  stop("Missing Stan model file(s) in ", normalizePath(models, mustWork = FALSE), ":\n  ",
       paste(paste0("stan_model_", missing, ".stan"), collapse = ", "),
       "\nCheck spelling of model_names in RunMe.R or add the .stan file(s).",
       call. = FALSE)
}

## Wrap fit_model() so any error -- e.g. cmdstanr's "No chains finished
## successfully" -- is captured as an error object instead of unwinding the
## worker. Errors thrown inside future_lapply otherwise propagate through
## the master and abort the whole batch with no per-model attribution.
safe_fit <- function(m, parallel) {
  tryCatch(fit_model(m, parallel = parallel),
           error = function(e) e)
}

## Parallel mode uses future::multisession + progressr for a per-model bar.
## Serial mode uses a plain for loop so each fit's cmdstanr per-chain progress
## prints directly to the console without any future_lapply stdout capture.
if (parallel_models) {
  n_per_fit <- max(1L, floor(n_cores / length(model_names)))
  cat("Fitting", length(model_names), "models in parallel (",
      paste(model_names, collapse = ", "), ") with",
      n_per_fit, "chain core(s) per model.\n")
  plan(multisession, workers = length(model_names))
  on.exit(plan(sequential), add = TRUE)
  handlers(global = TRUE)
  handlers("progress")

  fits <- with_progress({
    p <- progressor(steps = length(model_names))
    future_lapply(
      model_names,
      function(m) {
        source("load_packages.R")
        out <- safe_fit(m, parallel = n_per_fit)
        p(sprintf("[%s] done", m))
        out
      },
      future.globals = TRUE,
      future.seed    = TRUE,
      future.stdout  = NA
    )
  })
} else {
  n_per_fit <- n_cores
  fits <- vector("list", length(model_names))
  for (i in seq_along(model_names)) {
    m <- model_names[i]
    cat(sprintf("\n[%d/%d] fitting %s...\n", i, length(model_names), m))
    fits[[i]] <- safe_fit(m, parallel = n_per_fit)
  }
}

names(fits) <- model_names

## Surface failures with model name + underlying error message. Abort only
## after reporting every failed model so the user sees the full picture.
failed <- vapply(fits, inherits, logical(1), what = "error")
if (any(failed)) {
  cat("\n", sum(failed), " of ", length(model_names),
      " model fit(s) failed:\n", sep = "")
  for (m in names(fits)[failed]) {
    cat("  [", m, "] ", conditionMessage(fits[[m]]), "\n", sep = "")
  }
  stop(sum(failed), " model fit(s) failed (see messages above).",
       call. = FALSE)
}

cat("Fitted models:", paste(model_names, collapse = ", "), "\n")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
