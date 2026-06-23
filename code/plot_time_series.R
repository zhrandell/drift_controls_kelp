## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Time-series visualization of model-averaged remaining prey ~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Defines plot_timeseries_model_average(): reads model-average ODE trajectories
## from results/tmp/ODE_kelp_{low,high}_model_average.RDA (produced by
## compute_model_average() in model_average.R) and writes:
##   figs/timeseries_remaining_model_average.pdf     (2 rows x 2 cols)
##
## Reads from caller environment: tmp, figs, A.level
## ci_level is read from caller if present, otherwise defaults to 0.95.

plot_timeseries_model_average <- function(
  S0_targets = c(20, 100, 200, 300)   ## representative initial drift values (g)
) {

  ci_level_ts <- if (exists("ci_level", inherits = TRUE)) ci_level else 0.95
  tail_p      <- (1 - ci_level_ts) / 2
  ci_probs    <- c(tail_p, 1 - tail_p)

  ## blue gradient scaled to however many S0_targets are requested
  S0_cols <- colorRampPalette(c("#9ECAE1", "#08306B"))(length(S0_targets))

  ## graphical constants (mirrors plot_simulation.R)
  alph  <- 0.35
  linew <- 1

  ## ---- 1. Load period lengths from first available level ------------------- ##
  load(paste0(tmp, "/ODE_kelp_", names(A.level)[1], "_model_average.RDA"))
  ## P1, P2, P3, S0, A0, F0, outs_parms now in scope

  ## ---- 2. Build trajectory data frame for all kelp levels and S0 targets --- ##

  traj_list <- list()
  A0_vals   <- setNames(numeric(length(A.level)), names(A.level))

  for (level_name in names(A.level)) {

    load(paste0(tmp, "/ODE_kelp_", level_name, "_model_average.RDA"))
    ## reloads outs_parms, P1, P2, P3, S0, A0, F0 for this kelp level

    A0_val              <- A0[1]   ## actual initial kelp used in the simulation
    A0_vals[level_name] <- A0_val

    S0_rep_idx <- sapply(S0_targets, function(s) which.min(abs(S0 - s)))

    for (ii in seq_along(S0_rep_idx)) {

      init_idx <- S0_rep_idx[ii]
      S0_val   <- S0[init_idx]

      ## stack all posterior-draw matrices for this init condition
      ## outs_parms[[parm_idx]][[init_idx]]: n_time x 4 matrix [time, S, A, F]
      all_mats <- lapply(outs_parms, function(x) x[[init_idx]])
      arr      <- simplify2array(all_mats)   ## n_time x 4 x n_parms

      times <- arr[, 1, 1]   ## time vector (identical across draws)

      ## helper: median + CI across draws for a matrix (time x n_parms)
      summ <- function(mat) {
        list(
          med = apply(mat, 1, median),
          lo  = apply(mat, 1, quantile, ci_probs[1]),
          hi  = apply(mat, 1, quantile, ci_probs[2])
        )
      }

      sS <- summ(arr[, 2, ])
      sA <- summ(arr[, 3, ])

      df <- data.frame(
        time       = times,
        S0_val     = S0_val,
        kelp_level = level_name,
        S0_rank    = ii,
        S_med = sS$med, S_lo = sS$lo, S_hi = sS$hi,
        A_med = sA$med, A_lo = sA$lo, A_hi = sA$hi
      )

      ## prepend t = 0 row: prey at initial values
      t0       <- df[1, ]
      t0$time  <- 0
      t0$S_med <- S0_val; t0$S_lo <- S0_val; t0$S_hi <- S0_val
      t0$A_med <- A0_val; t0$A_lo <- A0_val; t0$A_hi <- A0_val

      traj_list[[length(traj_list) + 1]] <- rbind(t0, df)
    }
  }

  traj_df <- do.call(rbind, traj_list)

  ## ordered factor for legend (low -> mid -> high S0)
  S0_labs <- paste0(
    round(S0[sapply(S0_targets, function(s) which.min(abs(S0 - s)))]),
    " g"
  )
  traj_df$S0_label <- factor(
    paste0(round(traj_df$S0_val), " g"),
    levels = S0_labs
  )

  ## ---- 3. Shared graphical elements ---------------------------------------- ##

  t_total <- P1 + P2 + P3

  kelp_titles <- setNames(
    paste0(
      paste0(toupper(substring(names(A0_vals), 1, 1)), substring(names(A0_vals), 2)),
      " kelp (A0 = ", A0_vals, " g)"
    ),
    names(A0_vals)
  )

  restock_x <- c(P1, P1 + P2)

  x_breaks <- c(0, restock_x, t_total)

  ## axis-suppression helpers (mirrors plot_simulation.R style)
  both_blank <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                      axis.title.x = element_blank(), axis.text.x = element_blank())
  x_blank    <- theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  y_blank    <- theme(axis.title.y = element_blank(), axis.text.y = element_blank())

  ## ---- 4. Panel-building helper -------------------------------------------- ##

  make_ts_panel <- function(df_sub, y_med, y_lo, y_hi,
                            y_lab, title) {

    title_layer <- if (is.null(title)) {
      theme(plot.title = element_blank())
    } else {
      ggplot2::labs(title = title)
    }

    ggplot(df_sub, aes(x = time, group = S0_label)) +
      my.theme +
      geom_vline(xintercept = restock_x,
                 linetype = "dashed", color = "grey55", linewidth = 0.5) +
      geom_ribbon(aes(ymin = {{ y_lo }}, ymax = {{ y_hi }},
                      fill = S0_label),
                  alpha = alph) +
      geom_line(aes(y = {{ y_med }}, color = S0_label),
                linewidth = linew) +
      scale_fill_manual(
        values = setNames(S0_cols, S0_labs),
        name   = expression(S[0]~"(g)")
      ) +
      scale_color_manual(
        values = setNames(S0_cols, S0_labs),
        name   = expression(S[0]~"(g)")
      ) +
      scale_x_continuous(limits = c(0, t_total), breaks = x_breaks) +
      coord_cartesian(ylim = c(0, NA)) +
      xlab("Time (hours)") +
      ylab(y_lab) +
      title_layer
  }

  ## ============================================================
  ## Figure: remaining prey over time (2 rows x 2 cols)
  ## ============================================================

  panels_f2 <- list()
  for (lvl in names(A.level)) {
    sub <- traj_df[traj_df$kelp_level == lvl, ]
    panels_f2[[paste0("S_", lvl)]] <- make_ts_panel(
      sub, S_med, S_lo, S_hi,
      "Remaining drift (g)", kelp_titles[lvl])

    panels_f2[[paste0("A_", lvl)]] <- make_ts_panel(
      sub, A_med, A_lo, A_hi,
      "Remaining kelp (g)", NULL)
  }

  panels_f2$S_low  <- panels_f2$S_low  + x_blank
  panels_f2$S_high <- panels_f2$S_high + both_blank
  panels_f2$A_high <- panels_f2$A_high + y_blank

  fig2 <- wrap_plots(
    panels_f2$S_low,  panels_f2$S_high,
    panels_f2$A_low,  panels_f2$A_high,
    nrow = 2, ncol = 2
  ) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = rel(1.2), face = "plain"),
          legend.position = "right")

  ggplot2::ggsave(
    filename = paste0(figs, "/timeseries_remaining_model_average.pdf"),
    plot = fig2, dpi = 1200, width = 11, height = 6, units = "in"
  )
  invisible(NULL)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
