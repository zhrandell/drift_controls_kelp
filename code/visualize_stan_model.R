## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Inspect and visualize STAN output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Defines visualize_model(model_name): produces per-model trace / pairs /
## posterior / preference plots for the fit previously written by fit_model().
## Reads `results`, `figs`, `models` from the caller's environment. The
## preference form is selected by model name via model_family() in
## preference_helpers.R (names starting with "Logistic" or "vanLeeuwen").

## label map for known parameters; unknown ones get their raw name. Defined at
## script scope so compare_models.R (sourced after this file) can reuse it for
## the summary_priors_posteriors.tex table.
param_labels <- list(
  a     = expression("Search rate "          (italic(a))),
  b     = expression("Supression rate "      (italic(b))),
  w     = expression("Baseline preference "  (tilde(italic("ω")))),
  q     = expression("Switching sensitivity " (italic("φ"))),
  s     = expression("Stomach sensitivity "  (italic(v))),
  z     = expression("Stomach clearance "    (italic(z))),
  sigma = expression("Variance "             (italic(sigma)))
)

## Render a param_labels plot expression as a LaTeX-friendly string for
## table cells: "Search rate " (italic(a)) -> "Search rate ($a$)".
math_to_tex <- function(x) {
  if (is.call(x)) {
    fn    <- as.character(x[[1]])
    inner <- math_to_tex(x[[2]])
    if (fn == "tilde") return(paste0("\\tilde{", inner, "}"))
    return(inner)  # italic() etc. -- math mode handles italics
  }
  greek <- c("ω" = "\\omega",
             "φ" = "\\varphi",
             "sigma"  = "\\sigma")
  s <- as.character(x)
  if (s %in% names(greek)) return(unname(greek[s]))
  s
}
label_to_tex <- function(e) {
  if (is.null(e)) return(NA_character_)
  ex <- e[[1]]
  paste0(trimws(as.character(ex[[1]])), " ($", math_to_tex(ex[[2]]), "$)")
}

## Parse the names declared in a Stan file's `parameters { ... }` block.
## Returns the base identifiers only (no indices), e.g. c("a", "b", ...).
parse_param_block <- function(stan_file) {
  src <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  src <- gsub("//[^\n]*",   "", src)               # strip line comments
  src <- gsub("/\\*.*?\\*/", "", src, perl = TRUE) # strip block comments
  m   <- regmatches(src, regexpr("\\bparameters\\s*\\{[^}]*\\}", src, perl = TRUE))
  if (length(m) == 0) return(character(0))
  body  <- sub("^parameters\\s*\\{", "", m)
  body  <- sub("\\}$", "", body)
  decls <- strsplit(body, ";", fixed = TRUE)[[1]]
  out   <- character(0)
  for (d in decls) {
    d  <- trimws(d)
    if (nchar(d) == 0) next
    nm <- regmatches(d, regexpr("[A-Za-z_][A-Za-z0-9_]*\\s*$", d))
    if (length(nm) > 0) out <- c(out, trimws(nm))
  }
  out
}

visualize_model <- function(model_name) {

## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## specify directory and open saved RDS files
fit <- readRDS(paste0(tmp, "/model_output_", model_name, ".RDS"))

# rdat <- read.csv(paste0(data, "/drift_kelp_loss.csv"))
rdat <- readRDS(paste0(tmp, "/loss_dat.RData"))
rdat <- data.frame(rbind(rdat$y1_init_s_a,
                         rdat$y2_init_s_a,
                         rdat$y3_init_s_a,
                         rdat$y4_init_s_a,
                         rdat$y5_init_s_a
                         ))
colnames(rdat) <- c('Drift_Initial','Kelp_Initial','Urchins')
sel <- rdat$Drift_Initial > 0 & rdat$Kelp_Initial > 0
rdat <- rdat[sel ,]

## END start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## basic posterior check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## auto-detect sampled parameters by parsing the Stan source's `parameters`
## block. fit$metadata()$model_params is unreliable here -- it returns every
## variable in the CSV (transformed params, generated quantities, lp__).
parms  <- parse_param_block(paste0(models, "/stan_model_", model_name, ".stan"))
family <- model_family(model_name)  # "Logistic" or "vanLeeuwen"


## print param list output
print(fit, parms)


## extract draw information (sampled parameters only; pulling fit$draws() with
## no args includes log_lik, y_rep, theta, drift_*, kelp_*, etc., which inflates
## the array ~100x and triggers spurious bayesplot "dropped" warnings).
draws_array <- fit$draws(parms)
draws_df <- posterior::as_draws_df(draws_array)


## trace plot
t1 <- mcmc_trace(draws_array, pars = parms)

ggplot2::ggsave(filename = paste0(figs, "/trace_", model_name, ".pdf"), 
                plot = t1, 
                dpi = 1200, 
                width = 11,
                height = 6, 
                units = "in")


## pairs plot
pairsplot <- mcmc_pairs(draws_array,
                        pars = parms,
                        off_diag_args = list(size = 0.75))

ggplot2::ggsave(filename = paste0(figs, "/pairs_", model_name, ".pdf"), 
                plot = pairsplot, 
                dpi = 1200, 
                width = 11,
                height = 8, 
                units = "in")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## extract posteriors for subsequent simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
posts_df_raw <- as.data.frame(draws_df)[, parms, drop = FALSE]

## vanLeeuwen variants without a sampled `w` (e.g. vanLeeuwen_q) fix w = q - 4
## in the Stan model; derive it here so downstream preference math uses the
## same value the model was fit with.
if (family == "vanLeeuwen" && !("w" %in% parms) && "q" %in% parms) {
  posts_df_raw$w <- posts_df_raw$q - 4
}

## save RDA file with posteriors from a single chain
save(posts_df_raw, file = paste0(tmp, "/posterior_draws_", model_name, ".RDA"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Custom posterior plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function to visualize posteriors w/ CI and median values 
plot.posterior <- function(x, param='a', label='a', tag = '(a)'){
  ## graphical parameters
  col <- "#2FAA96"
  alp <- 1
  CI_col <- "black"
  med_col <- "black"
  lty1 <- 1
  lty2 <- 2
  sz1 <- 0.5
  sz2 <- 0.5

  CI <- quantile(data.frame(x)[param], c(0.0250, 0.5, 0.975), na.rm = TRUE)
  
  fig <- ggplot(x, aes(.data[[param]])) + 
    geom_histogram(fill = col, bins = 100) +
    # geom_density(fill = col, alpha = alp) +
    xlab(label) + 
    ylab('Density') + 
    geom_vline(xintercept = CI[1], 
               color = CI_col, linewidth = sz2, linetype = lty2) +
    geom_vline(xintercept = CI[2], 
               color = med_col, linewidth = sz1, linetype= lty1) + 
    geom_vline(xintercept = CI[3], 
               color = CI_col, linewidth = sz2, linetype = lty2) + 
    labs(tag = tag) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size=16),
          axis.text = element_text(size=14),
          plot.title = element_text(size=16))
  return(fig)
}


## one panel per sampled parameter
posterior_panels <- lapply(parms, function(p) {
  lbl <- if (!is.null(param_labels[[p]])) param_labels[[p]] else p
  plot.posterior(posts_df_raw, p, lbl, NULL)
})

n_panels  <- length(posterior_panels)
ncol_grid <- ceiling(sqrt(n_panels))
nrow_grid <- ceiling(n_panels / ncol_grid)

allparms <- wrap_plots(posterior_panels,
                       nrow = nrow_grid,
                       ncol = ncol_grid) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(size = rel(1.2), face = "plain"))

ggplot2::ggsave(filename = paste0(figs, "/posteriors_", model_name,".pdf"), 
                plot = allparms, 
                device = cairo_pdf,
                dpi = 1200, 
                width = 16,
                height = 8, 
                units = "in")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Plot preference-for-drift function
# For logistic, x = log(S / A).
# For vanLeeuwen, divide numerator and denominator by A and substitute S/A = exp(x).
# For van Leeuwen we use parameters 'w' and 'q' for convenience though in the notes we use \nu for w and \psi for q

## Preference plot requires `w` and `q` to be available in posts_df_raw.
## For vanLeeuwen variants without a sampled `w`, w was derived above as q - 4.
## Skip cleanly if the inputs still aren't there.
if (!all(c("w", "q") %in% names(posts_df_raw))) {
  message("[", model_name, "] preference plot skipped: parameters 'w' and 'q' not both available.")
  return(invisible(NULL))
}

## One-arg closures over the loop-local `w`, `q` (and `model_name`) so the
## existing `curve(Preference, ...)` and `LogSwitchPoint(sp)` call sites work
## unchanged. Math lives in preference_helpers.R, shared with compare_models.R.
Preference     <- function(x) preference(x, w, q, model_name)
LogSwitchPoint <- function(y) log_switch_point(y, w, q, model_name)

initial.ratios <- rdat$Drift_Initial / rdat$Kelp_Initial
initial.ratios <- log2(initial.ratios[is.finite(initial.ratios)])/log2(exp(1))

ratio.lim <- 10
ratio.step <- 0.1
ratio.vals <- log10(10^seq(-ratio.lim, ratio.lim, ratio.step))/log2(exp(1))
Pref.Predicts <- array(NA,
                     dim = c(nrow(posts_df_raw), length(ratio.vals)))
Pref.One2Ones <- vector()

for(i in 1:nrow(posts_df_raw)){
  w <- posts_df_raw$w[i]
  q <- posts_df_raw$q[i]
  Pref.Predicts[i,] <- Preference(ratio.vals)
  Pref.One2Ones[i] <- Preference(0)
}

Pref.Predictions <- apply(Pref.Predicts, 2, median)
Pref.One2One <- round(quantile(Pref.One2Ones, c(0.025, 0.5, 0.975)), 3)

pdf(paste0(figs, '/preference_', model_name,'.pdf'),
    height = 4,
    width = 8)
par(mar = c(3, 3, 1, 1),
    mgp = c(2, 0.2, 0),
    tcl = -0.1,
    xaxs = 'i',
    yaxs = 'i')

  xlims <- c(-6, 6)
  plot(1,1,
       xlim = xlims,
       ylim = c(0, 1),
       xlab = 'Relative abundance [Drift:Kelp]',
       ylab = 'Preference for Drift',
       type = 'n',
       axes = FALSE
  )
  
  rng <- range(initial.ratios)
  polygon(c(rng, rev(rng)), 
          c(0,0,1,1),
          border = NA,
          col = 'gray95')
  
  for(i in 1:nrow(posts_df_raw)){
    w <- posts_df_raw$w[i]
    q <- posts_df_raw$q[i]
    curve(Preference,
          min(xlims), max(xlims),
          add = TRUE,
          col = alpha('black', 0.1))
  }
  w <- median(posts_df_raw$w, na.rm = TRUE)
  q <- median(posts_df_raw$q, na.rm = TRUE)
  curve(Preference, 
        min(xlims), max(xlims), 
        add = TRUE,
        col = 'grey90',
        lwd = 5)
  curve(Preference, 
        min(xlims), max(xlims), 
        add = TRUE,
        col = 'black',
        lwd = 2)
  
  Po2o <- Pref.One2One[2]
  segments(c(0, 0), c(0, Po2o), 
           c(0, -10), c(Po2o, Po2o),
           lty = 1,
           lwd = 2,
           col = 'grey80')
  segments(c(0, 0), c(0, Po2o), 
           c(0, -10), c(Po2o, Po2o),
           lty = 3,
           lwd = 2,
           col = 'grey40')
  
  sp <- 0.5
  Lsp <- LogSwitchPoint(sp)
  segments(c(Lsp, Lsp), c(0, sp), 
           c(Lsp, -10), c(sp, sp),
           lty = 1,
           lwd = 2,
           col = 'grey80')
  segments(c(Lsp, Lsp), c(0, sp), 
           c(Lsp, -10), c(sp, sp),
           lty = 3,
           lwd = 2,
           col = 'grey40')
  
  x2.lim <- 10
  x2.step <- 2
  x2.vals <- 2^seq(-x2.lim, x2.lim, x2.step)
  x2.ats <- log2(x2.vals)/log2(exp(1))
  x2.labs <- c(rev(
    paste0('1:',2^seq(0, x2.lim, x2.step))), 
    paste0(2^seq(0, x2.lim, x2.step)[-1],':1'))
  axis(1, at = x2.ats, labels = x2.labs)
  axis(2, las = 1)
  box(lwd = 1)

dev.off()

invisible(NULL)
}  # end visualize_model()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


