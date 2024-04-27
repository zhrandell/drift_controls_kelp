## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit system of ODEs to empirical data in STAN ~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
  rm(list = ls())
  
  
## read in libraries
library(tidyverse)
library(stats)
library(rstan)
library(StanHeaders)
library(gridExtra)
library(gtable)
library(grid)
library(deSolve)
library(rstudioapi)
library(reshape2)
library(bayesplot)
library(hexbin)
library(shinystan)
library(janitor)
library(diffdf)
library(cmdstanr)
library(posterior)


## check wd is appropriate
getwd()


## hardcode relative file paths
code <- "../code"
data_input <- "../data_input"
data_output <- "../data_output"
figs <- "../figs"


## specify directory and open data
setwd(data_input)


## read in data
dat <- read.csv("drift_kelp_loss.csv", header = TRUE)
## END start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Data Configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subtract wet weight (20g) from drift bags
dat$Drift_Initial <- dat$Drift_Initial - 20
dat$Drift_Remaining <- dat$Drift_Remaining - 20
dat$Drift_Consumed <- dat$Drift_Initial - dat$Drift_Remaining

## substract wet weight (80g) from kelp bags 
dat$Kelp_Initial <- dat$Kelp_Initial - 80
dat$Kelp_Remaining <- dat$Kelp_Remaining - 80
dat$Kelp_Consumed <- dat$Kelp_Initial - dat$Kelp_Remaining

## set as factor
dat$Trial <- as.factor(dat$Trial)

## filter down to desired data: 
dat <- na.omit(dat)

## 24 hour data 
dat1 <- filter(dat, Trial %in% c("5","6")) 
dat1 <- filter(dat1, Treatment %in% c("Low", "High"))
dat3 <- filter(dat, Trial %in% c("7","8")) 
dat3 <- filter(dat3, Treatment %in% c("Low", "High"))
dat1 <- rbind(dat1, dat3)
remove(dat3)

## 48 hour data 
dat2 <- filter(dat, Trial %in% c("1","2","3","4")) 
dat2 <- filter(dat2, Period %in% c("1","2","3"))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






## assign unique key to urchin cohorts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## paste _1_ for "Low kelp" and _2_ for "High Kelp" 
dat1$key <- paste(dat1$Level, c("_2_", "_1_")[(dat1$Treatment=="Low")+1])
dat2$key <- paste(dat2$Level, c("_2_", "_1_")[(dat2$Treatment=="Low")+1])


## paste0 Trial # -- 1, 2, 3, 4, -- onto key column 
dat1$key <- with(dat1, paste0(key, Trial))
dat2$key <- with(dat2, paste0(key, Trial))


## remove white space from a character strings within a data frame  
dat1 <- as.data.frame(apply(dat1,2,function(x)gsub('\\s+', '',x)))
dat2 <- as.data.frame(apply(dat2,2,function(x)gsub('\\s+', '',x)))


## remove missing data; verify missing data with: 
dat1 <- dat1[!grepl("3_1_5", dat1$key), ]
dat1 <- dat1[!grepl("7_2_7", dat1$key), ]
dat2 <- dat2[!grepl("6_2_1", dat2$key), ]
dat2 <- dat2[!grepl("6_2_4", dat2$key), ]
dat2 <- dat2[!grepl("7_2_1", dat2$key), ]
dat2 <- dat2[!grepl("4_2_1", dat2$key), ]
dat2 <- dat2[!grepl("8_1_1", dat2$key), ]


## select Drift and Kelp initial conditions for each Period 
p1_1 <- filter(dat1, Period %in% c("1"))
p2_1 <- filter(dat1, Period %in% c("2"))
p3_2 <- filter(dat2, Period %in% c("1"))
p4_2 <- filter(dat2, Period %in% c("2"))
p5_2 <- filter(dat2, Period %in% c("3"))


## select focal columns 
p1_init <- p1_1[, c("key", "Drift_Initial","Kelp_Initial", "Drift_Consumed", 
                    "Drift_Remaining", "Kelp_Remaining", "Urchins")]
p2_init <- p2_1[, c("key", "Drift_Initial", "Kelp_Initial", "Drift_Consumed", 
                    "Drift_Remaining", "Kelp_Remaining", "Urchins")]
p3_init <- p3_2[, c("key", "Drift_Initial","Kelp_Initial", "Drift_Consumed", 
                    "Drift_Remaining", "Kelp_Remaining", "Urchins")]
p4_init <- p4_2[, c("key", "Drift_Initial", "Kelp_Initial", "Drift_Consumed", 
                    "Drift_Remaining", "Kelp_Remaining", "Urchins")]
p5_init <- p5_2[, c("key", "Drift_Initial", "Kelp_Initial", "Drift_Consumed", 
                    "Drift_Remaining", "Kelp_Remaining", "Urchins")]


## order data by key to for Stan to line up urchin cohorts
p1_init <- p1_init[order(p1_init$key), ]
p2_init <- p2_init[order(p2_init$key), ]
p3_init <- p3_init[order(p3_init$key), ]
p4_init <- p4_init[order(p4_init$key), ]
p5_init <- p5_init[order(p5_init$key), ]


## remove intermediary data frames 
remove(p1_1, p2_1, p3_2, p4_2, p5_2)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## initial conditions and time periods for Stan ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## number of subjects (replicates) we following through the three periods
n_subject_1 <- length(p1_init[,2]) 
n_subject_2 <- length(p3_init[,2])


## initial conditions for Drift S and Kelp A 
y1_init_SA <- as.matrix(p1_init[, c(2:3,7)])    ## Period 1
y2_init_SA <- as.matrix(p2_init[, c(2:3,7)])    ## Period 2
y3_init_SA <- as.matrix(p3_init[, c(2:3,7)])    ## Period 3
y4_init_SA <- as.matrix(p4_init[, c(2:3,7)])    ## Period 4
y5_init_SA <- as.matrix(p5_init[, c(2:3,7)])    ## Period 5 


## data from all three periods 
p1_Sremain <- p1_init[,5]
p2_Sremain <- p2_init[,5]
p1_Aremain <- p1_init[,6]
p2_Aremain <- p2_init[,6] 

p3_Sremain <- p3_init[,5]
p4_Sremain <- p4_init[,5]
p5_Sremain <- p5_init[,5]
p3_Aremain <- p3_init[,6]
p4_Aremain <- p4_init[,6]
p5_Aremain <- p5_init[,6]


## set dimensions 
col_1 <- 1
col_3 <- ncol(y1_init_SA)
nrow_1 <- n_subject_1
nrow_2 <- n_subject_2


## convert y1_init_SA from character matrix, to numeric, to numeric with proper dimensions
y1_init_SA <- mapply(y1_init_SA, FUN=as.numeric)
y1_init_SA <- matrix(data=y1_init_SA, ncol=col_3, nrow=nrow_1)

y2_init_SA <- mapply(y2_init_SA, FUN=as.numeric)
y2_init_SA <- matrix(data=y2_init_SA, ncol=col_3, nrow=nrow_1)

y3_init_SA <- mapply(y3_init_SA, FUN=as.numeric)
y3_init_SA <- matrix(data=y3_init_SA, ncol=col_3, nrow=nrow_2)

y4_init_SA <- mapply(y4_init_SA, FUN=as.numeric)
y4_init_SA <- matrix(data=y4_init_SA, ncol=col_3, nrow=nrow_2)

y5_init_SA <- mapply(y5_init_SA, FUN=as.numeric)
y5_init_SA <- matrix(data=y5_init_SA, ncol=col_3, nrow=nrow_2)


## convert drift consumed from character matrix, to numeric, to numeric with proper dimensions
## observations of drift remaining  
p1_Sremain <- mapply(p1_Sremain, FUN=as.numeric)
p1_Sremain <- matrix(data=p1_Sremain, ncol=col_1, nrow=nrow_1)
p2_Sremain <- mapply(p2_Sremain, FUN=as.numeric)
p2_Sremain <- matrix(data=p2_Sremain, ncol=col_1, nrow=nrow_1)

p3_Sremain <- mapply(p3_Sremain, FUN=as.numeric)
p3_Sremain <- matrix(data=p3_Sremain, ncol=col_1, nrow=nrow_2)
p4_Sremain <- mapply(p4_Sremain, FUN=as.numeric)
p4_Sremain <- matrix(data=p4_Sremain, ncol=col_1, nrow=nrow_2)
p5_Sremain <- mapply(p5_Sremain, FUN=as.numeric)
p5_Sremain <- matrix(data=p5_Sremain, ncol=col_1, nrow=nrow_2)

## observations of kelp remaining 
p1_Aremain <- mapply(p1_Aremain, FUN=as.numeric)
p1_Aremain <- matrix(data=p1_Aremain, ncol=col_1, nrow=nrow_1)
p2_Aremain <- mapply(p2_Aremain, FUN=as.numeric)
p2_Aremain <- matrix(data=p2_Aremain, ncol=col_1, nrow=nrow_1)

p3_Aremain <- mapply(p3_Aremain, FUN=as.numeric)
p3_Aremain <- matrix(data=p3_Aremain, ncol=col_1, nrow=nrow_2)
p4_Aremain <- mapply(p4_Aremain, FUN=as.numeric)
p4_Aremain <- matrix(data=p4_Aremain, ncol=col_1, nrow=nrow_2)
p5_Aremain <- mapply(p5_Aremain, FUN=as.numeric)
p5_Aremain <- matrix(data=p5_Aremain, ncol=col_1, nrow=nrow_2)


## bind resources remaining for Stan 
s0_1 <- cbind(p1_Sremain, p2_Sremain)
a0_1 <- cbind(p1_Aremain, p2_Aremain)
s0_2 <- cbind(p3_Sremain, p4_Sremain, p5_Sremain)
a0_2 <- cbind(p3_Aremain, p4_Aremain, p5_Aremain)


## row length per temporal sequence
dim.n_subject_1 <- length(y1_init_SA[,1])
dim.n_subject_2 <- length(y3_init_SA[,1])


## number of columns i.e. number of Periods per temporal sequence
dim.n_total_1 <- ncol(s0_1[,])
dim.n_total_2 <- ncol(s0_2[,])


## number of observations per each cohort per each of the three Periods.
nts1 <- 1
nts2 <- 1
nts3 <- 1
nts4 <- 1
nts5 <- 1


## sampling times for first initialization, sample Period 1, Period 2, and Period 3
time_1 <- c(1, 18, 36)
time_2 <- c(1, 44, 89, 134)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data list for Stan ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Loss.Dat <- list(
  n_subject_1 = n_subject_1,
  n_subject_2 = n_subject_2,
  nts1 = nts1, 
  nts2 = nts2,
  nts3 = nts3,
  nts4 = nts4,
  nts5 = nts5,
  t0_1 = time_1[1],
  t0_2 = time_2[1],
  ts1 = array(time_1[2], dim=c(1)),
  ts2 = array(time_1[3], dim=c(1)),
  ts3 = array(time_2[2], dim=c(1)), 
  ts4 = array(time_2[3], dim=c(1)), 
  ts5 = array(time_2[4], dim=c(1)), 
  y1_init_s_a = array(y1_init_SA, dim=c(dim.n_subject_1, col_3)), 
  y2_init_s_a = array(y2_init_SA, dim=c(dim.n_subject_1, col_3)),
  y3_init_s_a = array(y3_init_SA, dim=c(dim.n_subject_2, col_3)), 
  y4_init_s_a = array(y4_init_SA, dim=c(dim.n_subject_2, col_3)),
  y5_init_s_a = array(y5_init_SA, dim=c(dim.n_subject_2, col_3)),
  S_obs_1 = array(s0_1, dim=c(dim.n_subject_1, dim.n_total_1)),
  A_obs_1 = array(a0_1, dim=c(dim.n_subject_1, dim.n_total_1)),
  S_obs_2 = array(s0_2, dim=c(dim.n_subject_2, dim.n_total_2)),
  A_obs_2 = array(a0_2, dim=c(dim.n_subject_2, dim.n_total_2)))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## run Stan ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## load directory with Stan script
setwd(stan_code)


## compile model
model <- cmdstan_model(file <- "Restocking_2024.stan")


## initiate sampling 
fit <- model$sample(data = Loss.Dat, 
                  chains = 1,
                  iter_warmup = 1000,
                  iter_sampling = 2000,
                  adapt_delta=0.80, 
                  parallel_chains = 4) 





# set path to Model output
setwd(stanFit)

## save via cmdstan's preferred method
fit$save_object(file="fit_new_All_logNormalQ.RDS")

## open saved RDS file
fit <- readRDS("fit_new2_All.RDS")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## basic posterior check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
windows(h=8,w=11, record=TRUE)


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
