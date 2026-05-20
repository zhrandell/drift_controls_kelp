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
  readRDS(paste0(results, "/model_output_", m, ".RDS"))
}), model_names)

loo_list <- setNames(lapply(model_names, function(m) {
  fits[[m]]$loo(cores = n_cores)
}), model_names)

## ---- 2. Pairwise comparison ----------------------------------------------- ##
## loo_compare() / loo_model_weights() require >= 2 models. With a single model
## we just emit its ELPD summary; the diff / weight columns are not defined.

if (length(model_names) >= 2) {

  loo_cmp <- loo::loo_compare(loo_list)
  print(loo_cmp)

  ## Stacking / pseudo-BMA weights for model averaging
  stack_wts <- loo::loo_model_weights(loo_list, method = "stacking")
  pbma_wts  <- loo::loo_model_weights(loo_list, method = "pseudobma")

  cmp_df <- data.frame(
    model     = rownames(loo_cmp),
    elpd_diff = loo_cmp[, "elpd_diff"],
    se_diff   = loo_cmp[, "se_diff"],
    elpd_loo  = loo_cmp[, "elpd_loo"],
    se_elpd   = loo_cmp[, "se_elpd_loo"],
    p_loo     = loo_cmp[, "p_loo"],
    stacking_weight  = as.numeric(stack_wts[rownames(loo_cmp)]),
    pseudobma_weight = as.numeric(pbma_wts[rownames(loo_cmp)])
  )

} else {

  m <- model_names[1]
  est <- loo_list[[m]]$estimates    # rows: elpd_loo, p_loo, looic; cols: Estimate, SE
  print(loo_list[[m]])
  cmp_df <- data.frame(
    model     = m,
    elpd_diff = NA_real_,
    se_diff   = NA_real_,
    elpd_loo  = est["elpd_loo", "Estimate"],
    se_elpd   = est["elpd_loo", "SE"],
    p_loo     = est["p_loo",    "Estimate"],
    stacking_weight  = NA_real_,
    pseudobma_weight = NA_real_
  )

}

write.csv(cmp_df,
          file = paste0(results, "/model_comparison.csv"),
          row.names = FALSE)

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
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
