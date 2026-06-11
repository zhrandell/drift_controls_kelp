## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Compare fitted models via PSIS-LOO + posterior preds ~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Reads from the caller's environment (RunMe.R):
##   results, figs, n_cores, model_names, exclude_from_comparison

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
## `compare_names` is the subset of `model_names` kept in the LaTeX summary
## tables (section 2 below and section 5). LOO is still computed for every
## model in `model_names` above so sections 3-4 (Pareto-k, PPC) see all of
## them; only the table rows are filtered.
##
## The pairwise-comparison table is written as raw LaTeX (not via stargazer)
## so that math-mode markup ($...$, backslashes) can appear in the column
## labels; stargazer's `out=` path returns NA from nchar() on such headers
## and errors in .text.column.width (same workaround as
## visualize_stan_model.R).

compare_names <- setdiff(model_names, exclude_from_comparison)

if (length(compare_names) >= 2) {

  loo_cmp   <- loo::loo_compare(loo_list[compare_names])
  pbma_wts  <- loo::loo_model_weights(loo_list[compare_names],
                                      method = "pseudobma")
  print(loo_cmp)

  mods <- rownames(loo_cmp)
  ## Adopt loo_compare's best-to-worst ELPD order for every downstream summary
  ## so Summary_preference.tex rows match Summary_model_comparison.tex rows.
  compare_names <- mods
  vals <- data.frame(
    "ELPD"                     = sprintf("%.0f", loo_cmp[, "elpd_loo"]),
    "$SE_{ELPD}$"              = sprintf("%.0f", loo_cmp[, "se_elpd_loo"]),
    "$p$"                      = sprintf("%.1f", loo_cmp[, "p_loo"]),
    "$\\Delta \\text{ELPD}$"     = sprintf("%.1f", loo_cmp[, "elpd_diff"]),
    "$SE(\\Delta \\text{ELPD})$" = sprintf("%.1f", loo_cmp[, "se_diff"]),
    "Weight"                   = sprintf("%.2f", as.numeric(pbma_wts[mods])),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )

  cmp_tex <- data.frame(
    Model       = model_label(mods),
    vals,
    check.names = FALSE,
    row.names   = NULL
  )

  tab_caption <- paste0(
    "Relative model performance as assessed by the Bayesian LOO ",
    "estimate of the expected log pointwise predictive density.
    $p$ is the effective number of parameters.
    Model weights estimated using the pseudo-BMA method."
  )
  tab_label <- "tab:modelcomp_LOO"
  tab_path  <- paste0(tables, "/Summary_model_comparison.tex")
  tex_rows  <- apply(cmp_tex, 1, function(r) {
    paste0(paste(r, collapse = " & "), " \\\\")
  })
  writeLines(c(
    "\\begin{table}[!htbp] \\centering",
    paste0("  \\caption{", tab_caption, "}"),
    paste0("  \\label{",   tab_label,   "}"),
    "\\begin{tabular}{@{\\extracolsep{5pt}} lcccccc}",
    "\\\\[-1.8ex]\\hline",
    "\\hline \\\\[-1.8ex]",
    paste0(paste(colnames(cmp_tex), collapse = " & "), " \\\\"),
    "\\hline \\\\[-1.8ex]",
    tex_rows,
    "\\hline \\\\[-1.8ex]",
    "\\end{tabular}",
    "\\end{table}"
  ), tab_path)

} else if (length(compare_names) == 1) {
  print(loo_list[[compare_names[1]]])
} else {
  message("Summary_model_comparison.tex skipped: ",
          "all models were excluded via exclude_from_comparison.")
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

summary_rows <- lapply(compare_names, function(m) {
  parms  <- parse_param_block(paste0(models, "/stan_model_", m, ".stan"))
  family <- model_family(m)
  has_w  <- "w" %in% parms
  has_q  <- "q" %in% parms

  ## vanLeeuwen variants without a sampled `w` (e.g. vanLeeuwen_q) fix
  ## w = q - 4 in the Stan model; derive it here so the summary uses the
  ## same value the model was fit with.
  if (!has_q || (family == "Logistic" && !has_w)) {
    message("[", m, "] preference summary skipped: required parameters not sampled.")
    return(NULL)
  }
  cols  <- if (has_w) c("w", "q") else "q"
  draws <- posterior::as_draws_df(fits[[m]]$draws(cols))
  q <- draws$q
  w <- if (has_w) draws$w else q - 4

  pref     <- preference(0, w, q, m)
  switch_g <- 1 / exp(log_switch_point(0.5, w, q, m))
  logodds  <- log(pref / (1 - pref))

  c(
    fmt_med_ci(switch_g, digits = 2),
    fmt_med_ci(pref,     digits = 3),
    fmt_med_ci(logodds,  digits = 3)
  )
})

## Drop models that were skipped (e.g. missing w/q) so the row count matches.
keep       <- !vapply(summary_rows, is.null, logical(1))
summary_df <- data.frame(
  Model = model_label(compare_names[keep]),
  do.call(rbind, summary_rows[keep]),
  stringsAsFactors = FALSE,
  check.names      = FALSE,
  row.names        = NULL
)
colnames(summary_df) <- c(
  "Model",
  "Switch point",
  "Baseline preference",
  "Baseline log-odds"
)

if (nrow(summary_df) > 0) {
  ## Written as raw LaTeX (not via stargazer) so model_labels containing math
  ## markup ($...$) survive — same workaround as the comparison table above.
  pref_caption <- paste0(
    "Model-specific estimates (and 95\\% credible intervals) of the relative ",
    "abundance of drift versus kelp at which urchins preference for the two ",
    "resources is equal (expressed as $g$ of kelp per 1 $g$ of drift), ",
    "and of their baseline preference (expressed as proportional preference ",
    "and log-odds of consumption) for drift over kelp when the abundance of ",
    "the two resources is equal."
  )
  pref_label <- "tab:modelcomp_prefs"
  pref_path  <- paste0(tables, "/Summary_preference.tex")
  pref_rows  <- apply(summary_df, 1, function(r) {
    paste0(paste(r, collapse = " & "), " \\\\")
  })
  writeLines(c(
    "\\begin{table}[!htbp] \\centering",
    paste0("  \\caption{", pref_caption, "}"),
    paste0("  \\label{",   pref_label,   "}"),
    "\\begin{tabular}{@{\\extracolsep{5pt}} lccc}",
    "\\\\[-1.8ex]\\hline",
    "\\hline \\\\[-1.8ex]",
    paste0(paste(colnames(summary_df), collapse = " & "), " \\\\"),
    "\\hline \\\\[-1.8ex]",
    pref_rows,
    "\\hline \\\\[-1.8ex]",
    "\\end{tabular}",
    "\\end{table}"
  ), pref_path)
} else {
  message("Summary_preference.tex skipped: no eligible models after applying ",
          "exclude_from_comparison.")
}
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
