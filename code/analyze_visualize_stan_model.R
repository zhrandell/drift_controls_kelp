## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit system of ODEs to empirical data in STAN ~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
rm(list = ls())


## read in libraries
library(tidyverse)
library(cmdstanr)
library(rstan)
library(StanHeaders)
library(shinystan)
library(posterior)
library(bayesplot)


## check wd is appropriate
getwd()


## hardcode relative file paths
code <- "../code"
data_input <- "../data_input"
data_output <- "../data_output"
figs <- "../figs"


## specify directory and open saved RDS file
setwd(data_output)
fit <- readRDS("model_output_V1.RDS")


## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=16),
                 axis.text=element_text(size=14),
                 plot.title = element_text(size=16))
## END start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## basic posterior check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## list of params
parms <- c("a","v","q","p","sigma")


## print param list output
print(fit, parms)


## extract draw information
draws_array <- fit$draws()
str(draws_array)
draws_df <- posterior::as_draws_df(draws_array)


## plot posteriors
posts <- mcmc_hist(fit$draws(parms))
print(posts)


## sampling diagnostics 
diagnostic_df <- as_draws_df(fit$sampler_diagnostics()) 


## trace plot 
t1 <- mcmc_trace(draws_array, pars = parms) 
print(t1)  

setwd(msFigs)
ggplot2::ggsave(filename = "trace.eps", 
                plot = t1, 
                device = cairo_ps, 
                dpi = 1200, 
                width = 11,
                height = 8, 
                units = "in")


## pairs plot 
pairsplot <- mcmc_pairs(draws_array, pars = parms,
                        off_diag_args = list(size = 0.75))

setwd(msFigs)
ggplot2::ggsave(filename = "pairs.eps", 
                plot = pairsplot, 
                device = cairo_ps, 
                dpi = 1200, 
                width = 11,
                height = 8, 
                units = "in")


print(pairsplot)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## extract posterior for subsequent simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## extract a, v, p, q, from Chain 1 
posts_df_raw <- draws_df[c(1:10000),c(2:6)]


## redefine for posterior plot below
dat <- posts_df_raw

## save RDA file with posteriors from a single chain
setwd(modOutput)
save(posts_df_raw, file="posts_new_All.Rda")


## load RDA file  
#load("posts_df_raw.Rda")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom posteior plot with median and 95% CI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(egg)
CI <- apply(posts_df_raw, 2, quantile, c(0.0250, 0.975), na.rm=TRUE)
med <- apply(posts_df_raw, 2, median, na.rm=T)


## extract CI values 
a_lower <- CI[1,1]
a_upper <- CI[2,1]
v_lower <- CI[1,2]
v_upper <- CI[2,2]
q_lower <- CI[1,3]
q_upper <- CI[2,3]
p_lower <- CI[1,4]
p_upper <- CI[2,4]
sig_lower <- CI[1,5]
sig_upper <- CI[2,5]


## extract median values
a_med <- med[1]
v_med <- med[2]
q_med <- med[3]
p_med <- med[4]
sig_med <- med[5]


## graphical parameters
col <- "#2FAA96"
alp <- 1
CI_col <- "black"
med_col <- "black"
lty1 <- 1
lty2 <- 2
sz1 <- 0.5
sz2 <- 0.5


## plot
a_post <- ggplot(data = dat, aes(a)) + geom_density(fill=col, alpha=alp) +
  ggtitle("encounter rate \u03B1") + xlab("\u03B1") + my.theme + 
  geom_vline(xintercept = a_lower, color = CI_col, size = sz2, linetype=lty2) +
  geom_vline(xintercept = a_upper, color = CI_col, size = sz2, linetype=lty2) + 
  geom_vline(xintercept = a_med, color = med_col, size = sz1, linetype=lty1) 

q_post <- ggplot(data = dat, aes(q)) + geom_density(fill=col, alpha=alp) +
  ggtitle("resource preference \u03C6") + xlab("\u03C6") + my.theme + 
  geom_vline(xintercept = q_lower, color = CI_col, size = sz2, linetype=lty2) +
  geom_vline(xintercept = q_upper, color = CI_col, size = sz2, linetype=lty2) + 
  geom_vline(xintercept = q_med, color = med_col, size = sz1, linetype=lty1) +
  theme(axis.title.y = element_blank())

v_post <- ggplot(data = dat, aes(v)) + geom_density(fill=col, alpha=alp) +
  ggtitle("max gut fullness \u03B7") + xlab("\u03B7") + my.theme + 
  geom_vline(xintercept = v_lower, color = CI_col, size = sz2, linetype=lty2) +
  geom_vline(xintercept = v_upper, color = CI_col, size = sz2, linetype=lty2) + 
  geom_vline(xintercept = v_med, color = med_col, size = sz1, linetype=lty1) +
  theme(axis.title.y = element_blank())

p_post <- ggplot(data = dat, aes(p)) + geom_density(fill=col, alpha=alp) +
  ggtitle("gut clearance \u03B5") + xlab("\u03B5") + my.theme + 
  geom_vline(xintercept = p_lower, color = CI_col, size = sz2, linetype=lty2) +
  geom_vline(xintercept = p_upper, color = CI_col, size = sz2, linetype=lty2) + 
  geom_vline(xintercept = p_med, color = med_col, size = sz1, linetype=lty1) +
  theme(axis.title.y = element_blank())

sigma_post <- ggplot(data = dat, aes(sigma)) + geom_density(fill=col, alpha=alp) +
  ggtitle("variance \u03C3") + xlab("\u03C3") + my.theme + 
  geom_vline(xintercept = sig_lower, color = CI_col, size = sz2, linetype=lty2) +
  geom_vline(xintercept = sig_upper, color = CI_col, size = sz2, linetype=lty2) + 
  geom_vline(xintercept = sig_med, color = med_col, size = sz1, linetype=lty1) +
  theme(axis.title.y = element_blank())


windows(20,4,record=T)

all6 <- ggarrange(tag_facet(a_post + facet_wrap(~"time"), tag_pool="a"),
                  tag_facet(q_post + facet_wrap(~"time"), tag_pool="b"),
                  tag_facet(v_post + facet_wrap(~"time"), tag_pool="c"),
                  tag_facet(p_post + facet_wrap(~"time"), tag_pool="d"),
                  tag_facet(sigma_post + facet_wrap(~"time"), tag_pool="e"),
                  nrow=1, ncol=5)

print(all6)

setwd(msFigs)
ggplot2::ggsave(filename = "posts.eps", 
                plot = all6, 
                device = cairo_ps, 
                dpi = 1200, 
                width = 20,
                height = 4, 
                units = "in")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
