## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Inspect and visualize STAN output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##


## specify directory and open saved RDS file
fit <- readRDS(paste0(data_output, "/model_output.RDS"))


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

ggplot2::ggsave(filename = paste0(figs, "/trace.pdf"), 
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

ggplot2::ggsave(filename = paste0(figs, "/pairs.pdf"), 
                plot = pairsplot, 
                dpi = 1200, 
                width = 11,
                height = 8, 
                units = "in")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## extract posteriors for subsequent simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## extract a, v, p, q
posts_df_raw <- suppressWarnings(draws_df[1:min(nrow(draws_df), 1000), 
                                          1+(1:length(parms))])

## save RDA file with posteriors from a single chain
save(posts_df_raw, file = paste0(data_output, "/posts_new_All.RDA"))


## load RDA file  
#load("posts_df_raw.Rda")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Transform baseline preference parameter to Yodzis w form
InvLogit <- function(x){
  1 / (1 + exp(-x))
}
# Yodzis preference (for Drift)
posts_df_raw$w_y <- InvLogit(posts_df_raw$w)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Median and 95% CI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CI <- data.frame(apply(posts_df_raw, 2, quantile, c(0.0250, 0.5, 0.975), na.rm = TRUE))



## Custom posterior plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function to visualize posteriors w/ CI and median values 
plot.posterior <- function(data, param, label){
  ## graphical parameters
  col <- "#2FAA96"
  alp <- 1
  CI_col <- "black"
  med_col <- "black"
  lty1 <- 1
  lty2 <- 2
  sz1 <- 0.5
  sz2 <- 0.5

  CI <- quantile(data.frame(data)[param], c(0.0250, 0.5, 0.975), na.rm = TRUE)
  
  fig <- ggplot(data, aes(.data[[param]])) + 
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
wy_post <- plot.posterior(posts_df_raw, 'w_y', "Baseline preference (Yodzis w)")
q_post <- plot.posterior(posts_df_raw, 'q', "Switching rate (\u03C6)")
sigma_post <- plot.posterior(posts_df_raw, 'sigma', "Variance (\u03C3)")


allparms <- ggarrange(tag_facet(a_post + facet_wrap(~"time"), tag_pool = "a"),
                  tag_facet(v_post + facet_wrap(~"time"), tag_pool = "b"),
                  tag_facet(q_post + facet_wrap(~"time"), tag_pool = "c"),
                  tag_facet(w_post + facet_wrap(~"time"), tag_pool = "d"),
                  tag_facet(wy_post + facet_wrap(~"time"), tag_pool = "e"),
                  tag_facet(sigma_post + facet_wrap(~"time"), tag_pool = "f"),
                  nrow = 2, ncol = 3)

ggplot2::ggsave(filename = paste0(figs, "/posteriors.pdf"), 
                plot = allparms, 
                device = cairo_pdf,
                dpi = 1200, 
                width = 16,
                height = 4, 
                units = "in")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Plot preference function
LogisticPreference <- function(x){
    1 - (1 / (1 + exp(w + q * x) ))
}

ratio.lim <- 10
ratio.step <- 0.5
ratio.vals <- log2(2^seq(-ratio.lim, ratio.lim, ratio.step))/log2(exp(1))
Predictions <- array(NA, 
                     dim = c(nrow(posts_df_raw), length(ratio.vals)))

for(i in 1:nrow(posts_df_raw)){
  w <- posts_df_raw$w[i]
  q <- posts_df_raw$q[i]
  Predictions[i,] <- LogisticPreference(ratio.vals)
}
Prediction <- apply(Predictions, 2, mean)

# pdf(paste0(figs, '/preference.pdf'),
#     height = 4,
#     width = 8)
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
  
  segments(0, 0, 0, CI$w_y[2],
           lty = 3,
           col = 'grey')

  
  for(i in 1:nrow(posts_df_raw)){
    w <- posts_df_raw$w[i]
    q <- posts_df_raw$q[i]
    curve(LogisticPreference, 
          min(xlims), max(xlims), 
          add = TRUE,
          col = alpha('black', 0.1))
  }
  points(ratio.vals, Prediction,
         type = 'l',
         col = 'red',
         lwd = 2)

# dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

