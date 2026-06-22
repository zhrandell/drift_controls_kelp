## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Compare fitted models via PSIS-LOO + posterior preds ~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Reads from the caller's environment (RunMe.R):
##   results, figs, n_cores, model_names, exclude_from_comparison, weight_method

## Set the session-wide core count so loo and any nested parallel routines pick
## it up (loo_model_weights() has no explicit `cores` argument).
options(mc.cores = n_cores)

## ---- Helpers for prior/bounds extraction (section 6) ---------------------- ##
## Reuse parse_param_block() (defined in visualize_stan_model.R) for parameter
## names; these helpers extend it to also extract bounds and prior specs.

## Read a stan file with comments stripped, ready for regex parsing.
.read_stan_stripped <- function(stan_file) {
  src <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  src <- gsub("//[^\n]*",   "", src)               # strip line comments
  src <- gsub("/\\*.*?\\*/", "", src, perl = TRUE) # strip block comments
  src
}

## Parse the `parameters { ... }` block and return per-parameter bounds.
## Returns a data.frame with columns name, lower, upper (character).
parse_param_bounds <- function(stan_file) {
  src <- .read_stan_stripped(stan_file)
  m   <- regmatches(src, regexpr("\\bparameters\\s*\\{[^}]*\\}", src, perl = TRUE))
  if (length(m) == 0) return(data.frame(name = character(0),
                                        lower = character(0),
                                        upper = character(0),
                                        stringsAsFactors = FALSE))
  body  <- sub("^parameters\\s*\\{", "", m)
  body  <- sub("\\}$", "", body)
  decls <- strsplit(body, ";", fixed = TRUE)[[1]]
  out   <- data.frame(name = character(0), lower = character(0),
                      upper = character(0), stringsAsFactors = FALSE)
  for (d in decls) {
    d  <- trimws(d)
    if (nchar(d) == 0) next
    nm <- regmatches(d, regexpr("[A-Za-z_][A-Za-z0-9_]*\\s*$", d))
    if (length(nm) == 0) next
    lo <- regmatches(d, regexpr("lower\\s*=\\s*[^,>\\s]+", d, perl = TRUE))
    up <- regmatches(d, regexpr("upper\\s*=\\s*[^,>\\s]+", d, perl = TRUE))
    lo <- if (length(lo)) trimws(sub("^lower\\s*=\\s*", "", lo)) else NA_character_
    up <- if (length(up)) trimws(sub("^upper\\s*=\\s*", "", up)) else NA_character_
    out <- rbind(out, data.frame(name = trimws(nm), lower = lo, upper = up,
                                 stringsAsFactors = FALSE))
  }
  out
}

## Parse the `model { ... }` block and return per-parameter prior expressions.
## Returns a data.frame with columns name, prior (e.g. "exponential(10)").
## Only `<name> ~ <dist>(...)` statements on a single line are matched, where
## <name> is one of the declared parameters -- so likelihood lines like
## `target += normal_lpdf(...)` and indexed `~` statements are excluded.
parse_priors <- function(stan_file, param_names) {
  src <- .read_stan_stripped(stan_file)
  m   <- regmatches(src, regexpr("\\bmodel\\s*\\{.*?\\n\\}", src, perl = TRUE))
  if (length(m) == 0) {
    ## Fallback: model block may not end on its own line; grab up to the next
    ## top-level `}` by greedy match.
    m <- regmatches(src, regexpr("\\bmodel\\s*\\{[\\s\\S]*", src, perl = TRUE))
    if (length(m) == 0) return(data.frame(name = character(0),
                                          prior = character(0),
                                          stringsAsFactors = FALSE))
  }
  out <- data.frame(name = character(0), prior = character(0),
                    stringsAsFactors = FALSE)
  for (nm in param_names) {
    pat <- paste0("(?m)^\\s*\\Q", nm,
                  "\\E\\s*~\\s*([A-Za-z_][A-Za-z0-9_]*\\s*\\([^;]*\\))\\s*;")
    hit <- regmatches(m, regexpr(pat, m, perl = TRUE))
    if (length(hit) == 0) next
    pri <- sub(pat, "\\1", hit, perl = TRUE)
    out <- rbind(out, data.frame(name = nm, prior = trimws(pri),
                                 stringsAsFactors = FALSE))
  }
  out
}

## Bounds → "[lower, upper]" (plain text).
format_bounds <- function(lower, upper) {
  if (is.na(lower) && is.na(upper)) return("")
  paste0("[", ifelse(is.na(lower), "", lower), ", ",
              ifelse(is.na(upper), "", upper), "]")
}

## Stan distribution syntax → LaTeX math. Unknown distributions fall through
## wrapped in \text{...} so future families work without code changes.
format_prior <- function(stan_expr) {
  if (is.na(stan_expr) || !nzchar(stan_expr)) return("")
  m <- regmatches(stan_expr,
                  regexec("^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\((.*)\\)\\s*$",
                          stan_expr))[[1]]
  if (length(m) < 3) return(paste0("$\\text{", stan_expr, "}$"))
  fam  <- m[2]
  args <- trimws(m[3])
  tex <- switch(fam,
    exponential = paste0("$\\text{Exp}(", args, ")$"),
    normal      = paste0("$\\mathcal{N}(", args, ")$"),
    uniform     = paste0("$\\mathcal{U}(", args, ")$"),
    lognormal   = paste0("$\\text{LogNormal}(", args, ")$"),
    cauchy      = paste0("$\\text{Cauchy}(", args, ")$"),
    student_t   = paste0("$t(", args, ")$"),
    gamma       = paste0("$\\text{Gamma}(", args, ")$"),
    beta        = paste0("$\\text{Beta}(", args, ")$"),
    paste0("$\\text{", fam, "}(", args, ")$")
  )
  tex
}

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

  ## loo_model_weights() (via validate_psis_loo_list) requires every psis_loo
  ## object to have identical dim() = c(n_draws, n_obs). loo_compare() only
  ## checks n_obs, so models fit with different post-warmup iteration counts
  ## pass the comparison above but error here with "Each object in the list
  ## must have the same dimensions." Detect that case and emit the table
  ## without a Weight column instead of aborting the whole script.
  loo_dims  <- sapply(loo_list[compare_names], dim)
  same_dims <- length(unique(loo_dims[1L, ])) == 1L &&
               length(unique(loo_dims[2L, ])) == 1L
  if (same_dims) {
    model_wts <- loo::loo_model_weights(loo_list[compare_names],
                                        method = weight_method)
  } else {
    message("loo_model_weights skipped: models have mismatched dim() ",
            "(rows = n_draws, cols = n_obs); refit with matching iter ",
            "counts to recover ", weight_method, " weights:")
    print(loo_dims)
    model_wts <- setNames(rep(NA_real_, length(compare_names)),
                          compare_names)
  }
  print(loo_cmp)

  ## Sort models by stacking weight (highest first) for every downstream
  ## summary so Summary_preference.tex rows match Summary_model_comparison.tex
  ## rows. Falls back to loo_compare's best-to-worst ELPD order when weights
  ## are unavailable (mismatched dim() above).
  mods <- rownames(loo_cmp)
  if (all(is.finite(model_wts))) {
    mods <- names(sort(model_wts, decreasing = TRUE))
    loo_cmp <- loo_cmp[mods, , drop = FALSE]
  }
  compare_names <- mods
  vals <- data.frame(
    "ELPD"                     = sprintf("%.0f", loo_cmp[, "elpd_loo"]),
    "$SE_{ELPD}$"              = sprintf("%.0f", loo_cmp[, "se_elpd_loo"]),
    "$p$"                      = sprintf("%.1f", loo_cmp[, "p_loo"]),
    "$\\Delta \\text{ELPD}$"     = sprintf("%.1f", loo_cmp[, "elpd_diff"]),
    "$SE(\\Delta \\text{ELPD})$" = sprintf("%.1f", loo_cmp[, "se_diff"]),
    "Weight"                   = sprintf("%.2f", as.numeric(model_wts[mods])),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )

  cmp_tex <- data.frame(
    Model       = model_label(mods),
    vals,
    check.names = FALSE,
    row.names   = NULL
  )

  wt_method_label <- if (weight_method == "pseudobma") "pseudo-BMA+" else weight_method
  tab_caption <- paste0(
    "Relative model performance as assessed by the Bayesian LOO
    estimate of the expected log pointwise predictive density.
    $p$ is the effective number of parameters.
    Model weights estimated using the ", wt_method_label, " method.
    Models ordered by weight.
    A parenthetical $z$ indicates a model including a gut
    evacuation rate, assumed absent in other models.
    A parenthetical $m=1$ indicates a model with $m$ fixed to 1
    (equivalent to $b=0$), meaning no movement suppression at high resource
    abundances."
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
    fmt_med_ci(switch_g, digits = 1, ci_level = ci_level),
    fmt_med_ci(pref,     digits = 2, ci_level = ci_level),
    fmt_med_ci(logodds,  digits = 1, ci_level = ci_level)
  )
})

## Drop models that were skipped so the row count matches.
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
    "Model-specific median posterior estimates (and ", round(ci_level * 100), "\\% credible intervals) of the relative ",
    "abundance of drift versus kelp at which urchins preference for the two ",
    "resources is equal (expressed as $g$ of kelp per 1 $g$ of drift), ",
    "and of their baseline preference (expressed as proportional preference ",
    "and log-odds of consumption) for drift over kelp when the abundance of ",
    "the two resources is equal.  Models ordered and averaged by ", wt_method_label, " weight.
    See Table \\ref{tab:modelcomp_LOO} for model name interpretations."
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

## ---- 6. Combined priors, bounds, and posterior summary (LaTeX, all models)  ##
## For each model, scan the .stan source for declared parameter bounds and
## prior specifications, then summarize the corresponding posterior draws as
## the median and 95% credible interval. All models are emitted to a single
## summary_priors_posteriors.tex. Iterates in the order specified by
## `model_names` in RunMe.R (not the LOO best-to-worst order used by the
## comparison and preference tables), with `exclude_from_comparison` removed.

priors_posteriors_names <- setdiff(model_names, exclude_from_comparison)

if (length(priors_posteriors_names) >= 1) {

  ## Posterior median / CI formatter matches the per-model posterior tables
  ## that visualize_stan_model.R used to write: signif to 4 then 4 decimals.
  .fmt_post <- function(x) formatC(signif(x, 4), 4, format = "f")

  prior_post_rows <- do.call(rbind, lapply(priors_posteriors_names, function(m) {
    stan_file <- paste0(models, "/stan_model_", m, ".stan")
    bnds      <- parse_param_bounds(stan_file)
    prs       <- parse_priors(stan_file, bnds$name)
    df <- merge(bnds, prs, by = "name", all.x = TRUE, sort = FALSE)
    df <- df[match(bnds$name, df$name), , drop = FALSE]   # preserve decl order

    draws <- posterior::as_draws_df(fits[[m]]$draws(df$name))
    posts <- as.data.frame(draws)[, df$name, drop = FALSE]
    tail_p <- (1 - ci_level) / 2
    qs    <- apply(posts, 2, stats::quantile,
                   probs = c(tail_p, 0.5, 1 - tail_p), na.rm = TRUE)
    post_med <- vapply(qs[2, ], .fmt_post, character(1))
    post_ci  <- paste0("(", vapply(qs[1, ], .fmt_post, character(1)),
                       "---", vapply(qs[3, ], .fmt_post, character(1)), ")")

    par_lbls <- vapply(df$name, function(p) {
      lbl <- param_labels[[p]]
      if (is.null(lbl)) p else label_to_tex(lbl)
    }, character(1))
    data.frame(
      Model              = c(model_label(m), rep("", nrow(df) - 1L)),
      Parameter          = par_lbls,
      Prior              = vapply(df$prior, format_prior, character(1)),
      Bounds             = mapply(format_bounds, df$lower, df$upper),
      `Posterior median` = post_med,
      CI                 = post_ci,
      stringsAsFactors = FALSE,
      check.names      = FALSE,
      row.names        = NULL
    )
  }))

  ## Build tex rows, inserting an \hline before every model header row
  ## except the first so model groups read as distinct blocks.
  group_starts <- which(nzchar(prior_post_rows$Model))
  tex_rows <- character(0)
  for (i in seq_len(nrow(prior_post_rows))) {
    if (i %in% group_starts && i != group_starts[1]) {
      tex_rows <- c(tex_rows, "\\hline \\\\[-1.8ex]")
    }
    tex_rows <- c(tex_rows,
                  paste0(paste(unlist(prior_post_rows[i, ]),
                               collapse = " & "),
                         " \\\\"))
  }

  combined_caption <- paste0(
    "Parameter bounds, prior specifications, and posterior median estimates ",
    "with ", round(ci_level * 100), "\\% credible intervals.
    See Table \\ref{tab:modelcomp_LOO} for model name interpretations."
  )
  combined_label <- "tab:priors_posteriors"
  combined_path  <- paste0(tables, "/Summary_priors_posteriors.tex")
  ncol_tab       <- ncol(prior_post_rows)
  header_row     <- paste0(paste(colnames(prior_post_rows), collapse = " & "),
                           " \\\\")
  writeLines(c(
    "\\begingroup\\footnotesize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\begin{longtable}{@{} llllcc @{}}",
    paste0("\\caption{", combined_caption, "}\\label{",
           combined_label, "}\\\\"),
    "\\hline\\hline \\\\[-1.8ex]",
    header_row,
    "\\hline \\\\[-1.8ex]",
    "\\endfirsthead",
    paste0("\\multicolumn{", ncol_tab,
           "}{l}{\\textit{Table \\thetable{} continued from previous page}}\\\\"),
    "\\hline\\hline \\\\[-1.8ex]",
    header_row,
    "\\hline \\\\[-1.8ex]",
    "\\endhead",
    "\\hline \\\\[-1.8ex]",
    paste0("\\multicolumn{", ncol_tab,
           "}{r}{\\textit{Continued on next page}} \\\\"),
    "\\endfoot",
    "\\hline \\\\[-1.8ex]",
    "\\endlastfoot",
    tex_rows,
    "\\end{longtable}",
    "\\endgroup"
  ), combined_path)

} else {
  message("Summary_priors_posteriors.tex skipped: no eligible models after ",
          "applying exclude_from_comparison.")
}
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
