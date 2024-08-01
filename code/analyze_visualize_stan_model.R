## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit system of ODEs to empirical data in STAN ~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##


## specify directory and open saved RDS file
fit <- readRDS(paste0(data_output, "/model_output.RDS"))


## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title = element_text(size=16),
                 axis.text = element_text(size=14),
                 plot.title = element_text(size=16))
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





## custom posteior plot with median and 95% CI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
a_post <- ggplot(data = dat, aes(a)) + geom_density(fill = col, alpha = alp) +
  ggtitle("encounter rate \u03B1") + xlab("\u03B1") + my.theme + 
  geom_vline(xintercept = a_lower, color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = a_upper, color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = a_med, color = med_col, size = sz1, linetype = lty1) 

v_post <- ggplot(data = dat, aes(v)) + geom_density(fill = col, alpha = alp) +
  ggtitle("satiation sensitivity \u03B7") + xlab("\u03B7") + my.theme + 
  geom_vline(xintercept = v_lower, color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = v_upper, color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = v_med, color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())

p_post <- ggplot(data = dat, aes(p)) + geom_density(fill = col, alpha = alp) +
  ggtitle("gut clearance \u03B5") + xlab("\u03B5") + my.theme + 
  geom_vline(xintercept = p_lower, color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = p_upper, color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = p_med, color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())

w_post <- ggplot(data = dat, aes(w)) + geom_density(fill = col, alpha = alp) +
  ggtitle("baseline preference \u03c9") + xlab("\u03c9") + my.theme + 
  geom_vline(xintercept = w_lower, color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = w_upper, color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = w_med, color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())

q_post <- ggplot(data = dat, aes(q)) + geom_density(fill = col, alpha = alp) +
  ggtitle("switching rate \u03C6") + xlab("\u03C6") + my.theme + 
  geom_vline(xintercept = q_lower, color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = q_upper, color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = q_med, color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())

sigma_post <- ggplot(data = dat, aes(sigma)) + geom_density(fill = col, alpha = alp) +
  ggtitle("variance \u03C3") + xlab("\u03C3") + my.theme + 
  geom_vline(xintercept = sig_lower, color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = sig_upper, color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = sig_med, color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())


all6 <- ggarrange(tag_facet(a_post + facet_wrap(~"time"), tag_pool = "a"),
                  tag_facet(v_post + facet_wrap(~"time"), tag_pool = "b"),
                  tag_facet(p_post + facet_wrap(~"time"), tag_pool = "c"),
                  tag_facet(w_post + facet_wrap(~"time"), tag_pool = "d"),
                  tag_facet(q_post + facet_wrap(~"time"), tag_pool = "e"),
                  tag_facet(sigma_post + facet_wrap(~"time"), tag_pool = "f"),
                  nrow = 2, ncol = 3)


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

