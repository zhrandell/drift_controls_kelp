## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Process ODE simulation output for plotting ~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Defines process_model_sim(model_name): consumes
## tmp/ODE_kelp_<level>_<model_name>.RDA produced by simulate_model() and
## writes tmp/ODE_toPlot_kelp_<level>_<model_name>.RDA holding the per-
## period median, 95% CI, and initial-condition columns expected by
## plot_simulation.R. Reads `tmp`, `A.level` from the caller's environment.

process_model_sim <- function(model_name, ci_level = 0.95) {
  tail_p <- (1 - ci_level) / 2
  ci_probs <- c(tail_p, 1 - tail_p)

for (AL in 1:length(A.level)) { # initial kelp abundance (low and high)

  ## calculate and plot 95% CI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## loads outs_parms plus P1, P2, P3, S0, A0, F0 stashed by simulate_model().
  load(paste0(tmp, "/ODE_kelp_", names(A.level[AL]), '_', model_name, ".RDA"))

  ## functions to extract state values from list of lists
  extract_Sloss_P1 = function(y) {
    y[[1, 1]] <- y[[1, 2]] - y[[P1, 2]]
  }
  extract_Aloss_P1 = function(y) {
    y[[1, 1]] <- y[[1, 3]] - y[[P1, 3]]
  }
  extract_Ffill_P1 = function(y) {
    y[[1, 1]] <- (-1 * (y[[1, 4]] - y[[P1, 4]]))
  }

  extract_Sloss_P2 = function(y) {
    y[[1, 1]] <- y[[1, 2]] - y[[(P1 + P2), 2]]
  }
  extract_Aloss_P2 = function(y) {
    y[[1, 1]] <- y[[1, 3]] - y[[(P1 + P2), 3]]
  }
  extract_Ffill_P2 = function(y) {
    y[[1, 1]] <- (-1 * (y[[1, 4]] - y[[(P1 + P2), 4]]))
  }

  extract_Sloss_P3 = function(y) {
    y[[1, 1]] <- y[[1, 2]] - y[[(P1 + P2 + P3), 2]]
  }
  extract_Aloss_P3 = function(y) {
    y[[1, 1]] <- y[[1, 3]] - y[[(P1 + P2 + P3), 3]]
  }
  extract_Ffill_P3 = function(y) {
    y[[1, 1]] <- (-1 * (y[[1, 4]] - y[[(P1 + P2 + P3), 4]]))
  }


  ## extract all 9 state values in a single pass (avoids 9 separate list traversals)
  extracted_raw <- lapply(outs_parms, function(x) {
    vapply(x, function(y) {
      c(extract_Sloss_P1(y), extract_Aloss_P1(y), extract_Ffill_P1(y),
        extract_Sloss_P2(y), extract_Aloss_P2(y), extract_Ffill_P2(y),
        extract_Sloss_P3(y), extract_Aloss_P3(y), extract_Ffill_P3(y))
    }, numeric(9))  # returns 9 x len_init matrix per parameter set
  })

  ## stack into 9 x len_init x n_parms array, slice by variable into data.frames
  all_data <- simplify2array(extracted_raw)  # dims: 9 x len_init x n_parms
  make_df  <- function(k) na.omit(as.data.frame(t(all_data[k, , ])))

  S_loss_P1 <- make_df(1); A_loss_P1 <- make_df(2); F_fill_P1 <- make_df(3)
  S_loss_P2 <- make_df(4); A_loss_P2 <- make_df(5); F_fill_P2 <- make_df(6)
  S_loss_P3 <- make_df(7); A_loss_P3 <- make_df(8); F_fill_P3 <- make_df(9)


  ## rarify using a single consistent set of indices across all variables and periods
  rarify <- min(nrow(S_loss_P1), nrow(A_loss_P1), nrow(F_fill_P1),
                nrow(S_loss_P2), nrow(A_loss_P2), nrow(F_fill_P2),
                nrow(S_loss_P3), nrow(A_loss_P3), nrow(F_fill_P3), 990)
  idx <- sample(nrow(S_loss_P1), rarify)

  S_loss_P1 <- S_loss_P1[idx, ]; A_loss_P1 <- A_loss_P1[idx, ]; F_fill_P1 <- F_fill_P1[idx, ]
  S_loss_P2 <- S_loss_P2[idx, ]; A_loss_P2 <- A_loss_P2[idx, ]; F_fill_P2 <- F_fill_P2[idx, ]
  S_loss_P3 <- S_loss_P3[idx, ]; A_loss_P3 <- A_loss_P3[idx, ]; F_fill_P3 <- F_fill_P3[idx, ]



  ## calculate median from simulations
  S_P1 <- as.data.frame(matrix(apply(S_loss_P1, 2, median), ncol = 1))
  A_P1 <- as.data.frame(matrix(apply(A_loss_P1, 2, median), ncol = 1))
  F_P1 <- as.data.frame(matrix(apply(F_fill_P1, 2, median), ncol = 1))

  S_P2 <- as.data.frame(matrix(apply(S_loss_P2, 2, median), ncol = 1))
  A_P2 <- as.data.frame(matrix(apply(A_loss_P2, 2, median), ncol = 1))
  F_P2 <- as.data.frame(matrix(apply(F_fill_P2, 2, median), ncol = 1))

  S_P3 <- as.data.frame(matrix(apply(S_loss_P3, 2, median), ncol = 1))
  A_P3 <- as.data.frame(matrix(apply(A_loss_P3, 2, median), ncol = 1))
  F_P3 <- as.data.frame(matrix(apply(F_fill_P3, 2, median), ncol = 1))


  ## calculate CI at the requested ci_level
  S_loss_P1 <- apply(S_loss_P1, 2, quantile, ci_probs)
  A_loss_P1 <- apply(A_loss_P1, 2, quantile, ci_probs)
  F_fill_P1 <- apply(F_fill_P1, 2, quantile, ci_probs)

  S_loss_P2 <- apply(S_loss_P2, 2, quantile, ci_probs)
  A_loss_P2 <- apply(A_loss_P2, 2, quantile, ci_probs)
  F_fill_P2 <- apply(F_fill_P2, 2, quantile, ci_probs)

  S_loss_P3 <- apply(S_loss_P3, 2, quantile, ci_probs)
  A_loss_P3 <- apply(A_loss_P3, 2, quantile, ci_probs)
  F_fill_P3 <- apply(F_fill_P3, 2, quantile, ci_probs)


  ## bind df
  df <- as.data.frame(t(
    rbind(
      S_loss_P1,
      A_loss_P1,
      F_fill_P1,
      S_loss_P2,
      A_loss_P2,
      F_fill_P2,
      S_loss_P3,
      A_loss_P3,
      F_fill_P3,
      S0,
      A0,
      F0
    )
  ))

  names(df)[1] = "Smin_1"
  names(df)[2] = "Smax_1"
  names(df)[3] = "Amin_1"
  names(df)[4] = "Amax_1"
  names(df)[5] = "Fmin_1"
  names(df)[6] = "Fmax_1"
  names(df)[7] = "Smin_2"
  names(df)[8] = "Smax_2"
  names(df)[9] = "Amin_2"
  names(df)[10] = "Amax_2"
  names(df)[11] = "Fmin_2"
  names(df)[12] = "Fmax_2"
  names(df)[13] = "Smin_3"
  names(df)[14] = "Smax_3"
  names(df)[15] = "Amin_3"
  names(df)[16] = "Amax_3"
  names(df)[17] = "Fmin_3"
  names(df)[18] = "Fmax_3"


  ## bind data frames with median
  means <- as.data.frame(cbind(S_P1, A_P1, F_P1, S_P2, A_P2, F_P2, S_P3, A_P3, F_P3))

  names(means)[1] = "S_P1"
  names(means)[2] = "A_P1"
  names(means)[3] = "F_P1"
  names(means)[4] = "S_P2"
  names(means)[5] = "A_P2"
  names(means)[6] = "F_P2"
  names(means)[7] = "S_P3"
  names(means)[8] = "A_P3"
  names(means)[9] = "F_P3"

  ## END data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  varName <- paste0('combined_', names(A.level[AL]))
  assign(varName, cbind(df, means))
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} # end AL (kelp low or high) loop

save(combined_low,
     file = paste0(tmp, "/ODE_toPlot_kelp_low_",
                   model_name, ".RDA"))
save(combined_high,
     file = paste0(tmp, "/ODE_toPlot_kelp_high_",
                   model_name, ".RDA"))

invisible(NULL)
}  # end process_model_sim()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
