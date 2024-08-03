## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ inspect and visualize stan output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## specify directory and open saved RDS file
fit <- readRDS(paste0(data_output, "/model_output.RDS"))
## END start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## basic posterior check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## list of params
parms <- c("a", "v", "p", "w", "q", "sigma")


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


## save plot
ggplot2::ggsave(filename = paste0(figs, "/trace.eps"), 
                plot = t1, 
                device = cairo_ps, 
                dpi = 1200, 
                width = 11,
                height = 8, 
                units = "in")


## pairs plot 
pairsplot <- mcmc_pairs(draws_array, pars = parms,
                        off_diag_args = list(size = 0.75))

ggplot2::ggsave(filename = paste0(figs, "/pairs.eps"), 
                plot = pairsplot, 
                device = cairo_ps, 
                dpi = 1200, 
                width = 11,
                height = 8, 
                units = "in")


print(pairsplot)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## extract posteriors for subsequent simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## extract a, v, p, q
posts_df_raw <- suppressWarnings(draws_df[1:min(nrow(draws_df), 1000), 
                                          1+(1:length(parms))])


## redefine for posterior plot below
dat <- posts_df_raw

## save RDA file with posteriors from a single chain
save(posts_df_raw, file = paste0(data_output, "/posts_new_All.RDA"))


## load RDA file  
#load("posts_df_raw.Rda")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom posterior plot with median and 95% CI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CI <- data.frame(apply(posts_df_raw, 2, quantile, c(0.0250, 0.975), na.rm = TRUE))
med <- apply(posts_df_raw, 2, median, na.rm = T)


## extract CI values 
a_lower <- CI$a[1]
a_upper <- CI$a[2]
v_lower <- CI$v[1]
v_upper <- CI$v[2]
p_lower <- CI$p[1]
p_upper <- CI$p[2]
w_lower <- CI$w[1]
w_upper <- CI$w[2]
q_lower <- CI$q[1]
q_upper <- CI$q[2]
sig_lower <- CI$sigma[1]
sig_upper <- CI$sigma[2]


## extract median values
a_med <- med['a']
v_med <- med['v']
p_med <- med['p']
w_med <- med['w']
q_med <- med['q']
sig_med <- med['sigma']


## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title = element_text(size=16),
                 axis.text = element_text(size=14),
                 plot.title = element_text(size=16))


## graphical parameters
col <- "#2FAA96"
alp <- 1
CI_col <- "black"
med_col <- "black"
lty1 <- 1
lty2 <- 2
sz1 <- 0.5
sz2 <- 0.5


## function to visualize posteriors w/ CI and median values 
plot.posts <- function(data, param, text, label, lower_CI, upper_CI, median){
  
  fig <- ggplot(data, aes(param)) + geom_density(fill=col, alpha=alp) +
    ggtitle(text) + xlab(label) + my.theme + 
    
    geom_vline(xintercept = lower_CI, color = CI_col, size = sz2, linetype=lty2) +
    geom_vline(xintercept = upper_CI, color = CI_col, size = sz2, linetype=lty2) + 
    geom_vline(xintercept = median, color = med_col, size = sz1, linetype=lty1)
  
  return(fig)
}


## create posterior plots
a_post <- plot.posts(dat, dat$a, "encounter rate \u03B1", "\u03B1", a_lower, a_upper, a_med)
v_post <- plot.posts(dat, dat$v, "satiation sensitivity \u03B7", "\u03B7", v_lower, v_upper, v_med)
p_post <- plot.posts(dat, dat$p, "gut clearance \u03B5", "\u03B5", p_lower, p_upper, p_med)
w_post <- plot.posts(dat, dat$w, "baseline preference \u03c9", "\u03c9", w_lower, w_upper, w_med)
q_post <- plot.posts(dat, dat$q, "switching rate \u03C6", "\u03C6", q_lower, q_upper, q_med)
sigma_post <- plot.posts(dat, dat$sigma, "variance \u03C3", "\u03C3", sig_lower, sig_upper, sig_med)


## arrange all 6 posteriors in single ms figure
all6 <- ggarrange(tag_facet(a_post + facet_wrap(~"time"), tag_pool = "a"),
                  tag_facet(v_post + facet_wrap(~"time"), tag_pool = "b"),
                  tag_facet(p_post + facet_wrap(~"time"), tag_pool = "c"),
                  tag_facet(w_post + facet_wrap(~"time"), tag_pool = "d"),
                  tag_facet(q_post + facet_wrap(~"time"), tag_pool = "e"),
                  tag_facet(sigma_post + facet_wrap(~"time"), tag_pool = "f"),
                  nrow = 2, ncol = 3)


## save ms fig 
ggplot2::ggsave(filename = paste0(figs, "/posts.eps"), 
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
