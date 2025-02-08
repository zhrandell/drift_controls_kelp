## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Inspect and visualize STAN output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## specify directory and open saved RDS file
fit <- readRDS(paste0(results, "/model_output_", sel.model, ".RDS"))

rdat <- read.csv(paste0(data, "/drift_kelp_loss.csv"))

## END start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## basic posterior check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## list of params
parms <- c("a", "v", "w", "q", "sigma")


## print param list output
print(fit, parms)


## extract draw information
draws_array <- fit$draws()
draws_df <- posterior::as_draws_df(draws_array)


## plot posteriors
posts <- mcmc_hist(fit$draws(parms))
# print(posts)


## sampling diagnostics 
diagnostic_df <- as_draws_df(fit$sampler_diagnostics()) 


## trace plot 
t1 <- mcmc_trace(draws_array, pars = parms) 
print(t1)  

ggplot2::ggsave(filename = paste0(figs, "/trace_", sel.model, ".pdf"), 
                plot = t1, 
                dpi = 1200, 
                width = 11,
                height = 6, 
                units = "in")


## pairs plot 
pairsplot <- mcmc_pairs(draws_array, 
                        pars = parms,
                        off_diag_args = list(size = 0.75))
print(pairsplot)

ggplot2::ggsave(filename = paste0(figs, "/pairs_", sel.model, ".pdf"), 
                plot = pairsplot, 
                dpi = 1200, 
                width = 11,
                height = 8, 
                units = "in")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## extract posteriors for subsequent simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draws2pull <- 1000
posts_df_raw <- suppressWarnings(
  draws_df[1:min(nrow(draws_df), draws2pull), 
                                          1+(1:length(parms))]
  )

## save RDA file with posteriors from a single chain
save(posts_df_raw, file = paste0(results, "/posterior_draws_", sel.model, ".RDA"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Transform baseline preference parameter to Yodzis w form
Logit <- function(x){
  log(x / (1-x) )
}
InvLogit <- function(x){
  1 / (1 + exp( -x ))
}

if(sel.model=='Logistic'){
  # Yodzis preference (for Drift)
  posts_df_raw$w_y <- InvLogit(posts_df_raw$w)
}else{
  # vanLeewen preference (for Drift)
  posts_df_raw$w_y <- NA
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Median and 95% CI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CI <- data.frame(apply(posts_df_raw, 2, 
                       quantile, c(0.0250, 0.5, 0.975), 
                       na.rm = TRUE))



## Custom posterior plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function to visualize posteriors w/ CI and median values 
plot.posterior <- function(x, param='a', label='a'){
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
    geom_density(fill = col, alpha = alp) +
    xlab(label) + 
    ylab('Density') + 
    geom_vline(xintercept = CI[1], 
               color = CI_col, linewidth = sz2, linetype = lty2) +
    geom_vline(xintercept = CI[2], 
               color = med_col, linewidth = sz1, linetype= lty1) + 
    geom_vline(xintercept = CI[3], 
               color = CI_col, linewidth = sz2, linetype = lty2) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size=16),
          axis.text = element_text(size=14),
          plot.title = element_text(size=16))
  return(fig)
}

## create posterior plots
a_post <- plot.posterior(posts_df_raw, 'a', "Encounter rate (\u03B1)")
v_post <- plot.posterior(posts_df_raw, 'v', "Satiation sensitivity (\u03B7)")
w_post <- plot.posterior(posts_df_raw, 'w', "Baseline preference (\u03c9)")
wy_post <- plot.posterior(posts_df_raw, 'w_y', "Baseline preference (w)")
q_post <- plot.posterior(posts_df_raw, 'q', "Switching rate (\u03C6)")
sigma_post <- plot.posterior(posts_df_raw, 'sigma', "Variance (\u03C3)")


allparms <- ggarrange(tag_facet(a_post + facet_wrap(~"time"), tag_pool = "a"),
                  tag_facet(v_post + facet_wrap(~"time"), tag_pool = "b"),
                  tag_facet(q_post + facet_wrap(~"time"), tag_pool = "c"),
                  tag_facet(w_post + facet_wrap(~"time"), tag_pool = "d"),
                  tag_facet(wy_post + facet_wrap(~"time"), tag_pool = "e"),
                  tag_facet(sigma_post + facet_wrap(~"time"), tag_pool = "f"),
                  nrow = 2, ncol = 3)

ggplot2::ggsave(filename = paste0(figs, "/posteriors_", sel.model,".pdf"), 
                plot = allparms, 
                device = cairo_pdf,
                dpi = 1200, 
                width = 16,
                height = 4, 
                units = "in")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Plot preference-for-drift function 
# For logistic, x = log(S / A), 
# thus
# for vanLeeuwen, divide numerator and denominator by A and substitute S/A = exp(x)
if(sel.model == 'Logistic'){
  
  Preference <- function(x){
    1 - ( 1 / ( 1 + exp( w + q * x ) ) )
  }

  LogSwitchPoint <- function(y){
    ( log( -y / ( y - 1 ) ) - w ) / q
  }
  
}else{
  
  Preference <- function(x){
    1 - ( 1 + exp( w + x )) / 
      ( 1 + exp( log(2) + w + x ) + exp( q + 2 * x )) 
  }
  
  LogSwitchPoint <- function(y){
    log(- (2 * y) / 
     ( exp(w) * (2 * y - 1) - exp(q) * sqrt(exp(-2 * q) * (exp(2 * w) * (1 - 2 * y)^2 - 4 * exp(q) * (y - 1) * y)) ))
  }
  
}

initial.ratios <- rdat$Drift_Initial / rdat$Kelp_Initial
initial.ratios <- log2(initial.ratios[is.finite(initial.ratios)])/log2(exp(1))

ratio.lim <- 10
ratio.step <- 0.1
ratio.vals <- log10(10^seq(-ratio.lim, ratio.lim, ratio.step))/log2(exp(1))
Pref.Predicts <- array(NA,
                     dim = c(nrow(posts_df_raw), length(ratio.vals)))
Pref.One2Ones <- vector()
Switch.Predictions <- vector()

for(i in 1:nrow(posts_df_raw)){
  w <- posts_df_raw$w[i]
  q <- posts_df_raw$q[i]
  Pref.Predicts[i,] <- Preference(ratio.vals)
  Pref.One2Ones[i] <- Preference(0)
  Switch.Predictions[i] <- 1/exp(LogSwitchPoint(0.5))
}

Pref.Predictions <- apply(Pref.Predicts, 2, median)
Pref.One2One <- round(quantile(Pref.One2Ones, c(0.025, 0.5, 0.975)), 3)
Switch.Prediction <- round(quantile(Switch.Predictions, c(0.025, 0.5, 0.975)), 3)

pdf(paste0(figs, '/preference_', sel.model,'.pdf'),
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
  w <- CI$w[2]
  q <- CI$q[2]
  curve(Preference, 
        min(xlims), max(xlims), 
        add = TRUE,
        col = 'black',
        lwd = 5)
  curve(Preference, 
        min(xlims), max(xlims), 
        add = TRUE,
        col = 'grey',
        lwd = 3)
  
  # points(ratio.vals, Pref.Predictions,
  #        type = 'l',
  #        col = 'black',
  #        lwd = 5)
  # points(ratio.vals, Pref.Predictions,
  #        type = 'l',
  #        col = 'blue',
  #        lwd = 3)
  
  Po2o <- Pref.One2One[2]
  segments(c(0, 0), c(0, Po2o), 
           c(0, -10), c(Po2o, Po2o),
           lty = 3,
           lwd = 2,
           col = 'grey40')
  
  sp <- 0.5
  Lsp <- LogSwitchPoint(sp)
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

sink(paste0(results, "/Summary_", sel.model, ".txt"))
cat('Abundance switch point (equal preference):\n 1g drift to ', 
             round(1/exp(Lsp), 2), 
             'g (',
             Switch.Prediction[1],
             'g-',
             Switch.Prediction[3],
             'g) kelp.\n')
cat('Baseline preference (equal abundance):\n', 
             Pref.One2One[2], 
             ' (',
             Pref.One2One[1],
             '-',
             Pref.One2One[3],
             ') drift to ',
             1-Pref.One2One[2], 
             ' (',
             1-Pref.One2One[3],
             '-',
             1-Pref.One2One[1],
             ') kelp.\n')
cat('Baseline preference (equal abundance) log-odds:\n', 
    Logit(Pref.One2One[2]), 
    ' (',
    Logit(Pref.One2One[1]),
    '-',
    Logit(Pref.One2One[3]),
    ') drift to kelp.\n')
sink()
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


