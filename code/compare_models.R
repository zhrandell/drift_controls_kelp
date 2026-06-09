## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Compare fitted models via PSIS-LOO + posterior preds ~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Reads from the caller's environment (RunMe.R):
##   results, figs, n_cores, model_names

## Set the session-wide core count so loo and any nested parallel routines pick
## it up (loo_model_weights() has no explicit `cores` argument).
options(mc.cores = n_cores)

## ---- 1. Load fits and compute loo per model ------------------------------- ##

fits <- setNames(lapply(model_names, function(m) {
  readRDS(paste0(tmp, "/model_output_", m, ".RDS"))
}), model_names)

loo_list <- setNames(lapply(model_names, function(m) {
  fits[[m]]$loo(cores = n_cores)
}), model_names)

## ---- 2. Pairwise comparison (LaTeX export) -------------------------------- ##
## loo_compare() / loo_model_weights() require >= 2 models, so the comparison
## table is only meaningful then. With a single model we just print its ELPD
## summary to screen and skip the table.
##
## Plain ASCII column headers are required: stargazer's `out=` path runs an
## internal text-width pass that returns NA from nchar() when headers contain
## LaTeX math markers ($...$, backslashes), erroring in .text.column.width.

if (length(model_names) >= 2) {

  loo_cmp   <- loo::loo_compare(loo_list)
  stack_wts <- loo::loo_model_weights(loo_list, method = "stacking")
  pbma_wts  <- loo::loo_model_weights(loo_list, method = "pseudobma")
  print(loo_cmp)

  cmp_tex <- data.frame(
    Model            = rownames(loo_cmp),
    `ELPD diff`      = round(loo_cmp[, "elpd_diff"],    2),
    `SE diff`        = round(loo_cmp[, "se_diff"],      2),
    `ELPD loo`       = round(loo_cmp[, "elpd_loo"],     2),
    `SE ELPD loo`    = round(loo_cmp[, "se_elpd_loo"],  2),
    `p loo`          = round(loo_cmp[, "p_loo"],        2),
    `Stacking wt.`   = round(as.numeric(stack_wts[rownames(loo_cmp)]), 2),
    `Pseudo-BMA wt.` = round(as.numeric(pbma_wts[rownames(loo_cmp)]),  2),
    check.names      = FALSE,
    row.names        = NULL
  )

  invisible(capture.output(
    stargazer::stargazer(
      cmp_tex,
      summary  = FALSE,
      rownames = FALSE,
      out      = paste0(tables, "/Summary_model_comparison.tex")
    )
  ))

} else {
  print(loo_list[[model_names[1]]])
}

## ---- 3. Pareto-k diagnostics (flag unreliable obs per model) -------------- ##

pareto_k_summary <- data.frame(
  model         = model_names,
  pct_k_above_0p7 = vapply(model_names, function(m) {
    k <- loo_list[[m]]$diagnostics$pareto_k
    100 * mean(k > 0.7)
  }, numeric(1))
)
print(pareto_k_summary)

## ---- 4. Posterior predictive checks --------------------------------------- ##
## Build the observed-y vector that y_rep was generated against, in the same
## order as the Stan generated-quantities loops (drift_1, kelp_1, drift_2,
## kelp_2 - subjects x periods, row-major).

loss_dat <- readRDS(paste0(tmp, "/loss_dat.RData"))
y_obs <- c(
  as.numeric(t(loss_dat$S_obs_1)),
  as.numeric(t(loss_dat$A_obs_1)),
  as.numeric(t(loss_dat$S_obs_2)),
  as.numeric(t(loss_dat$A_obs_2))
)

for (m in model_names) {
  y_rep <- posterior::as_draws_matrix(fits[[m]]$draws("y_rep"))
  ## sample at most 200 draws for plotting speed
  keep <- sample.int(nrow(y_rep), size = min(200L, nrow(y_rep)))
  p <- bayesplot::ppc_dens_overlay(y = y_obs, yrep = y_rep[keep, , drop = FALSE]) +
    ggplot2::ggtitle(paste("PPC density overlay:", m))
  ggplot2::ggsave(paste0(figs, "/ppc_dens_", m, ".pdf"),
                  plot = p, width = 7, height = 5)
}

## ---- 5. Combined preference summary (LaTeX, all models) ------------------- ##
## For each model, compute the equal-preference switch point (g kelp per 1 g
## drift) and the equal-abundance baseline preference (probability, log-odds,
## odds).

## preference() and log_switch_point() (defined in preference_helpers.R) dispatch
## on the model name via model_family().

fmt_med_ci <- function(v, digits) {
  qs <- stats::quantile(v, c(0.025, 0.5, 0.975), na.rm = TRUE)
  sprintf("%.*f (%.*f, %.*f)",
          digits, qs[2], digits, qs[1], digits, qs[3])
}

summary_rows <- lapply(model_names, function(m) {
  parms <- parse_param_block(paste0(models, "/stan_model_", m, ".stan"))
  if (!all(c("w", "q") %in% parms)) {
    message("[", m, "] preference summary skipped: 'w' and 'q' not both sampled.")
    return(NULL)
  }
  draws <- posterior::as_draws_df(fits[[m]]$draws(c("w", "q")))
  w <- draws$w; q <- draws$q

  pref     <- preference(0, w, q, m)
  switch_g <- 1 / exp(log_switch_point(0.5, w, q, m))
  logodds  <- log(pref / (1 - pref))
  odds     <- exp(logodds)

  c(
    fmt_med_ci(switch_g, digits = 2),
    fmt_med_ci(pref,     digits = 3),
    fmt_med_ci(logodds,  digits = 3),
    fmt_med_ci(odds,     digits = 2)
  )
})

## Drop models that were skipped (e.g. missing w/q) so the row count matches.
keep       <- !vapply(summary_rows, is.null, logical(1))
summary_df <- data.frame(
  Model = model_names[keep],
  do.call(rbind, summary_rows[keep]),
  stringsAsFactors = FALSE,
  check.names      = FALSE,
  row.names        = NULL
)
colnames(summary_df) <- c(
  "Model",
  "Switch point (g kelp per 1 g drift)",
  "Baseline preference (drift)",
  "Baseline log-odds (drift to kelp)",
  "Baseline odds (drift to kelp)"
)

invisible(capture.output(
  stargazer::stargazer(
    summary_df,
    summary  = FALSE,
    rownames = FALSE,
    out      = paste0(tables, "/Summary_preference.tex")
  )
))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
