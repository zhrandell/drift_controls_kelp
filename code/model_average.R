## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~ Model averaging via stacking weights ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Defines compute_model_average(model_wts, compare_names, fits, n_avg):
## pools ODE trajectories and (w, q) posterior draws proportionally to
## stacking weights and produces:
##   results/tmp/ODE_kelp_{low,high}_model_average.RDA
##   results/tmp/ODE_toPlot_kelp_{low,high}_model_average.RDA
##   figs/ODE_simulation_model_average.pdf
##   figs/preference_model_average.pdf
##   tables/Summary_preference.tex  (Model average row prepended)
##
## Reads from the caller's environment:
##   tmp, figs, tables, models, A.level, n_sim_draws
##   fmt_med_ci, preference, log_switch_point, model_family (preference_helpers.R)
##   parse_param_block (visualize_stan_model.R)
##   process_model_sim (process_simulation.R)
##   plot_model_sim, my.theme (plot_simulation.R)

compute_model_average <- function(model_wts, compare_names, fits,
                                  n_avg = n_sim_draws,
                                  plot_curves = TRUE) {

  ## Distribute n_total draws among models by weight, applying a rounding
  ## correction to the highest-weight model so the total is exact.
  alloc_draws <- function(wts, n_total) {
    n           <- round(wts * n_total)
    top         <- which.max(wts)
    n[top]      <- n_total - sum(n[-top])
    n
  }

  wts        <- as.numeric(model_wts[compare_names])
  names(wts) <- compare_names

  ## ---- A. Pool ODE outputs ------------------------------------------------ ##

  n_per_model <- alloc_draws(wts, n_avg)
  ref_model   <- compare_names[which.max(wts)]

  for (level_name in names(A.level)) {

    ## load shared simulation parameters from the top-weighted model
    load(paste0(tmp, "/ODE_kelp_", level_name, "_", ref_model, ".RDA"))
    P1_ref <- P1; P2_ref <- P2; P3_ref <- P3
    S0_ref <- S0; A0_ref <- A0; F0_ref <- F0

    pooled_outs <- list()
    for (k in seq_along(compare_names)) {
      n_k <- n_per_model[k]
      if (n_k <= 0) next
      load(paste0(tmp, "/ODE_kelp_", level_name, "_", compare_names[k], ".RDA"))
      n_avail     <- length(outs_parms)
      idx         <- sample.int(n_avail, size = n_k, replace = n_k > n_avail)
      pooled_outs <- c(pooled_outs, outs_parms[idx])
    }

    outs_parms <- pooled_outs
    P1 <- P1_ref; P2 <- P2_ref; P3 <- P3_ref
    S0 <- S0_ref; A0 <- A0_ref; F0 <- F0_ref
    save(outs_parms, P1, P2, P3, S0, A0, F0,
         file = paste0(tmp, "/ODE_kelp_", level_name, "_model_average.RDA"))

    cat(sprintf("[model_average] pooled %d ODE draws for kelp %s\n",
                length(pooled_outs), level_name))
  }

  process_model_sim("model_average", ci_level = ci_level)
  plot_model_sim("model_average")

  ## ---- B. Pool (w, q) draws ----------------------------------------------- ##

  n_pref_per_model <- alloc_draws(wts, n_avg)

  w_pool   <- numeric(0)
  q_pool   <- numeric(0)
  fam_pool <- character(0)

  for (k in seq_along(compare_names)) {
    m   <- compare_names[k]
    n_k <- n_pref_per_model[k]
    if (n_k <= 0) next

    parms  <- parse_param_block(paste0(models, "/stan_model_", m, ".stan"))
    family <- model_family(m)
    has_w  <- "w" %in% parms
    has_q  <- "q" %in% parms
    if (!has_q || (family == "Logistic" && !has_w)) {
      message("[model_average] preference skipped for ", m,
              ": required parameters not sampled.")
      next
    }

    cols  <- if (has_w) c("w", "q") else "q"
    draws <- posterior::as_draws_df(fits[[m]]$draws(cols))
    q_m   <- draws$q
    w_m   <- if (has_w) draws$w else q_m - 4

    n_avail  <- length(q_m)
    idx      <- sample.int(n_avail, size = n_k, replace = n_k > n_avail)
    w_pool   <- c(w_pool,   w_m[idx])
    q_pool   <- c(q_pool,   q_m[idx])
    fam_pool <- c(fam_pool, rep(m, n_k))
  }

  if (length(w_pool) == 0) {
    message("compute_model_average: no eligible models for preference averaging.")
    return(invisible(NULL))
  }

  ## ---- C. Model-averaged preference statistics ---------------------------- ##

  ## dispatch preference() and log_switch_point() by model family per draw
  pref_draws <- numeric(length(w_pool))
  lsp_draws  <- numeric(length(w_pool))
  for (fam in unique(fam_pool)) {
    idx              <- which(fam_pool == fam)
    pref_draws[idx]  <- preference(0, w_pool[idx], q_pool[idx], fam)
    lsp_draws[idx]   <- log_switch_point(0.5, w_pool[idx], q_pool[idx], fam)
  }
  switch_g <- 1 / exp(lsp_draws)
  logodds  <- log(pref_draws / (1 - pref_draws))

  avg_row <- c(
    "Model average",
    fmt_med_ci(switch_g,   digits = 1, ci_level = ci_level),
    fmt_med_ci(pref_draws, digits = 2, ci_level = ci_level),
    fmt_med_ci(logodds,    digits = 1, ci_level = ci_level)
  )

  ## prepend "Model average" row to Summary_preference.tex
  pref_path <- paste0(tables, "/Summary_preference.tex")
  if (file.exists(pref_path)) {
    lines   <- readLines(pref_path)
    hdr_idx <- grep("Switch point", lines, fixed = TRUE)
    if (length(hdr_idx) == 1L) {
      ins_after <- hdr_idx + 1L   # line after the \hline that follows column headers
      avg_tex   <- paste0(paste(avg_row, collapse = " & "), " \\\\")
      lines <- c(
        lines[seq_len(ins_after)],
        avg_tex,
        "\\hline \\\\[-1.8ex]",
        lines[seq(ins_after + 1L, length(lines))]
      )
      writeLines(lines, pref_path)
      cat("[model_average] Model average row prepended to Summary_preference.tex\n")
    }
  }

  ## ---- D. Preference curve plot ------------------------------------------- ##

  xlims  <- c(-6, 6)
  x_grid <- seq(xlims[1], xlims[2], length.out = 300)
  pref_mat <- matrix(NA_real_, nrow = length(w_pool), ncol = length(x_grid))
  for (fam in unique(fam_pool)) {
    idx             <- which(fam_pool == fam)
    pref_mat[idx, ] <- sapply(x_grid, function(x)
      preference(x, w_pool[idx], q_pool[idx], fam))
  }

  ## summary quantities reused from section C
  Po2o <- median(pref_draws)
  Lsp  <- median(lsp_draws)

  pref_med <- apply(pref_mat, 2, median, na.rm = TRUE)

  ## initial drift:kelp ratios from experimental data (for the range polygon)
  loss_dat       <- readRDS(paste0(tmp, "/loss_dat.RData"))
  rdat_avg       <- data.frame(rbind(loss_dat$y1_init_s_a, loss_dat$y2_init_s_a,
                                     loss_dat$y3_init_s_a, loss_dat$y4_init_s_a,
                                     loss_dat$y5_init_s_a))
  colnames(rdat_avg) <- c('Drift_Initial', 'Kelp_Initial', 'Urchins')
  rdat_avg       <- rdat_avg[rdat_avg$Drift_Initial > 0 & rdat_avg$Kelp_Initial > 0, ]
  init_ratios    <- log2(rdat_avg$Drift_Initial / rdat_avg$Kelp_Initial) / log2(exp(1))
  init_ratios    <- init_ratios[is.finite(init_ratios)]

  pdf(paste0(figs, '/preference_model_average.pdf'),
      height = 4,
      width  = 8)
  par(mar  = c(3, 3, 1, 1),
      mgp  = c(2, 0.2, 0),
      tcl  = -0.1,
      xaxs = 'i',
      yaxs = 'i')

    plot(1, 1,
         xlim = xlims,
         ylim = c(0, 1),
         xlab = 'Relative abundance [Drift:Kelp]',
         ylab = 'Preference for Drift',
         type = 'n',
         axes = FALSE)

    rng <- range(init_ratios)
    polygon(c(rng, rev(rng)),
            c(0, 0, 1, 1),
            border = NA,
            col    = 'gray95')

    if (plot_curves) {
      n_draws <- nrow(pref_mat)
      for (i in seq_len(n_draws)) {
        lines(x_grid, pref_mat[i, ],
              col = alpha('black', 100 / n_draws))
      }
    } else {
      tail_p  <- (1 - ci_level) / 2
      pref_lo <- apply(pref_mat, 2, quantile, tail_p,       na.rm = TRUE)
      pref_hi <- apply(pref_mat, 2, quantile, 1 - tail_p,   na.rm = TRUE)
      polygon(c(x_grid, rev(x_grid)),
              c(pref_lo, rev(pref_hi)),
              col    = adjustcolor('steelblue', alpha.f = 0.3),
              border = NA)
    }
    lines(x_grid, pref_med, col = 'grey90', lwd = 5)
    lines(x_grid, pref_med, col = 'black',  lwd = 2)

    segments(c(0, 0), c(0, Po2o),
             c(0, -10), c(Po2o, Po2o),
             lty = 1, lwd = 2, col = 'grey80')
    segments(c(0, 0), c(0, Po2o),
             c(0, -10), c(Po2o, Po2o),
             lty = 3, lwd = 2, col = 'grey40')

    segments(c(Lsp, Lsp), c(0, 0.5),
             c(Lsp, -10), c(0.5, 0.5),
             lty = 1, lwd = 2, col = 'grey80')
    segments(c(Lsp, Lsp), c(0, 0.5),
             c(Lsp, -10), c(0.5, 0.5),
             lty = 3, lwd = 2, col = 'grey40')

    x2.lim  <- 10
    x2.step <- 2
    x2.vals <- 2^seq(-x2.lim, x2.lim, x2.step)
    x2.ats  <- log2(x2.vals) / log2(exp(1))
    x2.labs <- c(rev(paste0('1:', 2^seq(0, x2.lim, x2.step))),
                 paste0(2^seq(0, x2.lim, x2.step)[-1], ':1'))
    axis(1, at = x2.ats, labels = x2.labs)
    axis(2, las = 1)
    box(lwd = 1)

  dev.off()

  cat("[model_average] Done.\n")
  invisible(NULL)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
