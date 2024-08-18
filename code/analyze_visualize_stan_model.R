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

# Transform baseline preference parameter to Yodzis w form
InvLogit <- function(x){
  1 / (1 + exp(-x))
}
# Yodzis preference (for Drift)
dat$w_y <- InvLogit(dat$w)
posts_df_raw$w_y <- InvLogit(posts_df_raw$w)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## custom posteior plot with median and 95% CI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CI <- data.frame(apply(posts_df_raw, 2, quantile, c(0.0250, 0.5, 0.975), na.rm = TRUE))

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
  xlab("encounter rate \u03B1") + my.theme + theme(axis.text.y = element_blank()) + 
  geom_vline(xintercept = CI$a[1], color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = CI$a[3], color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = CI$a[2], color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())

v_post <- ggplot(data = dat, aes(v)) + geom_density(fill = col, alpha = alp) +
  xlab("satiation sensitivity \u03B7") + my.theme + theme(axis.text.y = element_blank()) + 
  geom_vline(xintercept = CI$v[1], color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = CI$v[3], color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = CI$v[2], color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())

w_post <- ggplot(data = dat, aes(w)) + geom_density(fill = col, alpha = alp) +
  xlab("baseline preference \u03c9") + my.theme + theme(axis.text.y = element_blank()) + 
  geom_vline(xintercept = CI$w[1], color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = CI$w[3], color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = CI$w[2], color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())

wy_post <- ggplot(data = dat, aes(w_y)) + geom_density(fill = col, alpha = alp) +
  xlab("baseline preference (Yodzis) w") + my.theme + theme(axis.text.y = element_blank()) + 
  geom_vline(xintercept = CI$w_y[1], color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = CI$w_y[3], color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = CI$w_y[2], color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())

q_post <- ggplot(data = dat, aes(q)) + geom_density(fill = col, alpha = alp) +
  xlab("switching rate \u03C6") + my.theme + theme(axis.text.y = element_blank()) + 
  geom_vline(xintercept = CI$q[1], color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = CI$q[3], color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = CI$q[2], color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())

sigma_post <- ggplot(data = dat, aes(sigma)) + geom_density(fill = col, alpha = alp) +
  xlab("variance \u03C3") + my.theme + theme(axis.text.y = element_blank()) + 
  geom_vline(xintercept = CI$sigma[1], color = CI_col, size = sz2, linetype = lty2) +
  geom_vline(xintercept = CI$sigma[3], color = CI_col, size = sz2, linetype = lty2) + 
  geom_vline(xintercept = CI$sigma[2], color = med_col, size = sz1, linetype = lty1) +
  theme(axis.title.y = element_blank())


allparms <- ggarrange(
  tag_facet(a_post, tag_pool = "a"),
  tag_facet(v_post, tag_pool = "b"),
  tag_facet(w_post, tag_pool = "c"),
  tag_facet(q_post, tag_pool = "d"),
  tag_facet(sigma_post, tag_pool = "e"),
  tag_facet(wy_post, tag_pool = "f"),
  nrow = 2, ncol = 3)


ggplot2::ggsave(filename = paste0(figs, "/posts.eps"), 
                plot = allparms, 
                device = cairo_ps, 
                dpi = 1200, 
                width = 20,
                height = 4, 
                units = "in")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Plot preference function
LogisticPreference <- function(x){
    1 - (1 / (1 + exp(w + q * x) ))
}
xlims <- c(-6, 6)
plot(1,1,
     xlim = xlims,
     ylim = c(0, 1),
     xlab = 'log(Drift/Kelp)',
     ylab = 'Preference for drift',
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


abline(v = 0, 
       lty = 3,
       col = 'grey')
for(i in 1:nrow(dat)){
  w <- dat$w[i]
  q <- dat$q[i]
  curve(LogisticPreference, min(xlims), max(xlims), 
        add = TRUE,
        col = alpha('black', 0.1))
}
abline(h = CI$w_y[2],
       col = 'grey20',
       lwd = 2,
       lty = 2)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

