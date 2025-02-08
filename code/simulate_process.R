
for(AL in 1:length(A.level)){ # initial kelp abundance (low and high)
  
  ## calculate and plot 95% CI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  load(paste0(results,"/ODE_kelp_", names(A.level[AL]), '_', sel.model, ".RDA"))
  
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
  
  
  ## Drift loss
  S_loss_P1 <- lapply(outs_parms, function(x)
    lapply(x, function(y) extract_Sloss_P1(y)))
  A_loss_P1 <- lapply(outs_parms, function(x)
    lapply(x, function(y) extract_Aloss_P1(y)))
  F_fill_P1 <- lapply(outs_parms, function(x)
    lapply(x, function(y) extract_Ffill_P1(y)))
  
  
  S_loss_P2 <- lapply(outs_parms, function(x)
    lapply(x, function(y) extract_Sloss_P2(y)))
  A_loss_P2 <- lapply(outs_parms, function(x)
    lapply(x, function(y) extract_Aloss_P2(y)))
  F_fill_P2 <- lapply(outs_parms, function(x)
    lapply(x, function(y) extract_Ffill_P2(y)))
  
  
  S_loss_P3 <- lapply(outs_parms, function(x)
    lapply(x, function(y) extract_Sloss_P3(y)))
  A_loss_P3 <- lapply(outs_parms, function(x)
    lapply(x, function(y) extract_Aloss_P3(y)))
  F_fill_P3 <- lapply(outs_parms, function(x)
    lapply(x, function(y) extract_Ffill_P3(y)))
  
  
  ## convert to df 
  S_loss_P1 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, S_loss_P1)), 
                                            ncol=len_init)))
  A_loss_P1 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, A_loss_P1)), 
                                            ncol=len_init)))
  F_fill_P1 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, F_fill_P1)), 
                                            ncol=len_init)))
  
  S_loss_P2 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, S_loss_P2)), 
                                            ncol=len_init)))
  A_loss_P2 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, A_loss_P2)), 
                                            ncol=len_init)))
  F_fill_P2 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, F_fill_P2)), 
                                            ncol=len_init)))
  
  S_loss_P3 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, S_loss_P3)), 
                                            ncol=len_init)))
  A_loss_P3 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, A_loss_P3)), 
                                            ncol=len_init)))
  F_fill_P3 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, F_fill_P3)), 
                                            ncol=len_init)))
  
  
  ## rarify to compensate for NA's produced during simulation 
  rarify <- 990
  S_loss_P1 <- S_loss_P1[sample(1:nrow(S_loss_P1), rarify), ]
  A_loss_P1 <- A_loss_P1[sample(1:nrow(A_loss_P1), rarify), ]
  F_fill_P1 <- F_fill_P1[sample(1:nrow(F_fill_P1), rarify), ]
  
  S_loss_P2 <- S_loss_P2[sample(1:nrow(S_loss_P1), rarify), ]
  A_loss_P2 <- A_loss_P2[sample(1:nrow(A_loss_P1), rarify), ]
  F_fill_P2 <- F_fill_P2[sample(1:nrow(F_fill_P1), rarify), ]
  
  S_loss_P3 <- S_loss_P3[sample(1:nrow(S_loss_P1), rarify), ]
  A_loss_P3 <- A_loss_P3[sample(1:nrow(A_loss_P1), rarify), ]
  F_fill_P3 <- F_fill_P3[sample(1:nrow(F_fill_P1), rarify), ]
  
  
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
  
  
  ## calculate 95 CI
  S_loss_P1 <- apply(S_loss_P1, 2, quantile, c(0.025, 0.975))
  A_loss_P1 <- apply(A_loss_P1, 2, quantile, c(0.025, 0.975))
  F_fill_P1 <- apply(F_fill_P1, 2, quantile, c(0.025, 0.975))
  
  S_loss_P2 <- apply(S_loss_P2, 2, quantile, c(0.025, 0.975))
  A_loss_P2 <- apply(A_loss_P2, 2, quantile, c(0.025, 0.975))
  F_fill_P2 <- apply(F_fill_P2, 2, quantile, c(0.025, 0.975))
  
  S_loss_P3 <- apply(S_loss_P3, 2, quantile, c(0.025, 0.975))
  A_loss_P3 <- apply(A_loss_P3, 2, quantile, c(0.025, 0.975))
  F_fill_P3 <- apply(F_fill_P3, 2, quantile, c(0.025, 0.975))
  
  
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
  
  ## use to plot both in final figure: high kelp simulation
  varName <- paste0('combined_', names(A.level[AL]))
  assign(varName, cbind(df, means))
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} # end AL (kelp low or high) loop

save(combined_low, 
     file = paste0(results, "/ODE_toPlot_kelp_low_", 
                   sel.model, ".RDA"))
save(combined_high, 
     file = paste0(results, "/ODE_toPlot_kelp_high_", 
                   sel.model, ".RDA"))



