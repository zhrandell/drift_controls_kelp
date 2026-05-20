## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit system of ODEs to empirical data in STAN ~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##


## read in data
dat <- read.csv(paste0(data, "/drift_kelp_loss.csv"), header = TRUE)
## END start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Data Configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subtract drift bag wet weight (20g)
dat$Drift_Initial <- pmax(0, dat$Drift_Initial - 20)
dat$Drift_Remaining <- pmax(0, dat$Drift_Remaining - 20)
dat$Drift_Consumed <- pmax(0, dat$Drift_Initial - dat$Drift_Remaining)

## subtract kelp bag wet weight (80g)
dat$Kelp_Initial <- pmax(0, dat$Kelp_Initial - 80)
dat$Kelp_Remaining <- pmax(0, dat$Kelp_Remaining - 80)
dat$Kelp_Consumed <- pmax(0, dat$Kelp_Initial - dat$Kelp_Remaining)

## set as factor
dat$Trial <- as.factor(dat$Trial)

## filter down to desired data:
dat <- na.omit(dat)


dat$Drift_Initial[dat$Drift_Initial==0] <- 0
dat$Kelp_Initial[dat$Kelp_Initial==0] <- 0

keepTrtmts <- c("Low", "High") #, "Low_Control", "High_Control", "Drift_Control")

## 24 hour data 
dat1 <- filter(dat, Trial %in% c("5","6")) 
dat1 <- filter(dat1, Treatment %in% keepTrtmts)
dat3 <- filter(dat, Trial %in% c("7","8")) 
dat3 <- filter(dat3, Treatment %in% keepTrtmts)
dat1 <- rbind(dat1, dat3)
remove(dat3)

## 48 hour data 
dat2 <- filter(dat, Trial %in% c("1","2","3","4")) 
dat2 <- filter(dat2, Period %in% c("1","2","3"))
dat2 <- filter(dat2, Treatment %in% keepTrtmts)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END data filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Drift to Kelp ratio
DKratio <- c(dat1$Drift_Initial/dat1$Kelp_Initial,
             dat2$Drift_Initial/dat2$Kelp_Initial)
DKratio.rng <- range(DKratio[is.finite(DKratio) & DKratio!=0])


## assign unique key to urchin cohorts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to rename the treatment level to a numeric value to create proper key
relabel.levels <- function(data){
  data <- data %>%
    mutate(Level = case_when(
      Level == "Drift_Control" ~ "9",
      Level == "Low_Control" ~ "10",
      Level == "High_Control" ~ "11",
      TRUE ~ Level
    ))
  return(data)
}


## function to create key 
create.key <- function(data){
  data <- data %>%
    mutate(key = paste(Level, case_when(
      Treatment == "Low" ~ "_1_",
      Treatment == "Drift_Control" ~ "_3_",
      Treatment == "Low_Control" ~ "_4_",
      Treatment == "High_Control" ~ "_5_",
      TRUE ~ "_2_"  # Default case, if none of the conditions are met
    )))
  return(data)
}


## invoke function to rename levels to accommodate controls
dat1 <- relabel.levels(dat1)
dat2 <- relabel.levels(dat2)


## create key
dat1 <- create.key(dat1)
dat2 <- create.key(dat2)


## paste0 Trial onto key column
dat1$key <- with(dat1, paste0(key, Trial))
dat2$key <- with(dat2, paste0(key, Trial))


## remove white space from a character strings within a data frame
dat1 <- as.data.frame(apply(dat1, 2, function(x)
  gsub('\\s+', '', x)))
dat2 <- as.data.frame(apply(dat2, 2, function(x)
  gsub('\\s+', '', x)))
## END key creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 





## remove data via key ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## remove keys associated with missing data in dat1, the 24hr data
dat1 <- dat1[!grepl("3_1_5", dat1$key),] 
dat1 <- dat1[!grepl("9_3_8", dat1$key),] 


## remove keys associated with missing data in dat2, the 48hr data
dat2 <- dat2[!grepl("6_2_1", dat2$key),] 
dat2 <- dat2[!grepl("6_2_4", dat2$key),] 
dat2 <- dat2[!grepl("7_2_1", dat2$key),] 
dat2 <- dat2[!grepl("4_2_1", dat2$key),] 
dat2 <- dat2[!grepl("8_1_1", dat2$key),] 


## select Drift and Kelp initial conditions for each sampling Period
## 24 hour data
p1_1 <- filter(dat1, Period %in% c("1")) ## 24 hour data, period 1
p2_1 <- filter(dat1, Period %in% c("2")) ## 24 hour data, period 2


## 48 hour data  
p3_2 <- filter(dat2, Period %in% c("1")) ## 48 hour data, period 1
p4_2 <- filter(dat2, Period %in% c("2")) ## 48 hour data, period 2
p5_2 <- filter(dat2, Period %in% c("3")) ## 48 hour data, period 3


## select focal columns
p1_init <-
  p1_1[, c(
    "key",
    "Drift_Initial",
    "Kelp_Initial",
    "Drift_Consumed",
    "Drift_Remaining",
    "Kelp_Remaining",
    "Urchins"
  )]
p2_init <-
  p2_1[, c(
    "key",
    "Drift_Initial",
    "Kelp_Initial",
    "Drift_Consumed",
    "Drift_Remaining",
    "Kelp_Remaining",
    "Urchins"
  )]
p3_init <-
  p3_2[, c(
    "key",
    "Drift_Initial",
    "Kelp_Initial",
    "Drift_Consumed",
    "Drift_Remaining",
    "Kelp_Remaining",
    "Urchins"
  )]
p4_init <-
  p4_2[, c(
    "key",
    "Drift_Initial",
    "Kelp_Initial",
    "Drift_Consumed",
    "Drift_Remaining",
    "Kelp_Remaining",
    "Urchins"
  )]
p5_init <-
  p5_2[, c(
    "key",
    "Drift_Initial",
    "Kelp_Initial",
    "Drift_Consumed",
    "Drift_Remaining",
    "Kelp_Remaining",
    "Urchins"
  )]


## order data by key to for Stan to line up urchin cohorts
p1_init <- p1_init[order(p1_init$key),]
p2_init <- p2_init[order(p2_init$key),]
p3_init <- p3_init[order(p3_init$key),]
p4_init <- p4_init[order(p4_init$key),]
p5_init <- p5_init[order(p5_init$key),]


## remove intermediary data frames
remove(p1_1, p2_1, p3_2, p4_2, p5_2)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## initial conditions and time periods for Stan ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## number of subjects (replicates) we following through the three periods
n_subject_1 <- length(p1_init[, 2])
n_subject_2 <- length(p3_init[, 2])


## initial conditions for Drift S and Kelp A
y1_init_SA <- as.matrix(p1_init[, c(2:3, 7)])    ## Period 1
y2_init_SA <- as.matrix(p2_init[, c(2:3, 7)])    ## Period 2
y3_init_SA <- as.matrix(p3_init[, c(2:3, 7)])    ## Period 3
y4_init_SA <- as.matrix(p4_init[, c(2:3, 7)])    ## Period 4
y5_init_SA <- as.matrix(p5_init[, c(2:3, 7)])    ## Period 5


## data from all three periods
p1_Sremain <- p1_init[, 5]
p2_Sremain <- p2_init[, 5]

p1_Aremain <- p1_init[, 6]
p2_Aremain <- p2_init[, 6]

p3_Sremain <- p3_init[, 5]
p4_Sremain <- p4_init[, 5]
p5_Sremain <- p5_init[, 5]

p3_Aremain <- p3_init[, 6]
p4_Aremain <- p4_init[, 6]
p5_Aremain <- p5_init[, 6]


## set dimensions
col_1 <- 1
col_3 <- ncol(y1_init_SA)
nrow_1 <- n_subject_1
nrow_2 <- n_subject_2


## convert to numeric matrices

## initial drift and kelp
class(y1_init_SA) <- "numeric"
class(y2_init_SA) <- "numeric"
class(y3_init_SA) <- "numeric"
class(y4_init_SA) <- "numeric"
class(y5_init_SA) <- "numeric"

## observations of drift remaining
class(p1_Sremain) <- "numeric"
class(p2_Sremain) <- "numeric"
class(p3_Sremain) <- "numeric"
class(p4_Sremain) <- "numeric"
class(p5_Sremain) <- "numeric"

p1_Sremain <- matrix(data = p1_Sremain, ncol = col_1, nrow = nrow_1)
p2_Sremain <- matrix(data = p2_Sremain, ncol = col_1, nrow = nrow_1)
p3_Sremain <- matrix(data = p3_Sremain, ncol = col_1, nrow = nrow_2)
p4_Sremain <- matrix(data = p4_Sremain, ncol = col_1, nrow = nrow_2)
p5_Sremain <- matrix(data = p5_Sremain, ncol = col_1, nrow = nrow_2)

## observations of kelp remaining
class(p1_Aremain) <- "numeric"
class(p2_Aremain) <- "numeric"
class(p3_Aremain) <- "numeric"
class(p4_Aremain) <- "numeric"
class(p5_Aremain) <- "numeric"

p1_Aremain <- matrix(data = p1_Aremain, ncol = col_1, nrow = nrow_1)
p2_Aremain <- matrix(data = p2_Aremain, ncol = col_1, nrow = nrow_1)
p3_Aremain <- matrix(data = p3_Aremain, ncol = col_1, nrow = nrow_2)
p4_Aremain <- matrix(data = p4_Aremain, ncol = col_1, nrow = nrow_2)
p5_Aremain <- matrix(data = p5_Aremain, ncol = col_1, nrow = nrow_2)


## bind resources remaining for Stan
s0_1 <- cbind(p1_Sremain, p2_Sremain)
a0_1 <- cbind(p1_Aremain, p2_Aremain)
s0_2 <- cbind(p3_Sremain, p4_Sremain, p5_Sremain)
a0_2 <- cbind(p3_Aremain, p4_Aremain, p5_Aremain)


## row length per temporal sequence
dim.n_subject_1 <- length(y1_init_SA[, 1])
dim.n_subject_2 <- length(y3_init_SA[, 1])


## number of columns i.e. number of Periods per temporal sequence
dim.n_total_1 <- ncol(s0_1[, ])
dim.n_total_2 <- ncol(s0_2[, ])


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
loss_dat <- list(
  n_subject_1 = n_subject_1,
  n_subject_2 = n_subject_2,
  nts1 = nts1,
  nts2 = nts2,
  nts3 = nts3,
  nts4 = nts4,
  nts5 = nts5,
  t0_1 = time_1[1],
  t0_2 = time_2[1],
  ts1 = array(time_1[2], dim = c(1)),
  ts2 = array(time_1[3], dim = c(1)),
  ts3 = array(time_2[2], dim = c(1)),
  ts4 = array(time_2[3], dim = c(1)),
  ts5 = array(time_2[4], dim = c(1)),
  y1_init_s_a = array(y1_init_SA, dim = c(dim.n_subject_1, col_3)),
  y2_init_s_a = array(y2_init_SA, dim = c(dim.n_subject_1, col_3)),
  y3_init_s_a = array(y3_init_SA, dim = c(dim.n_subject_2, col_3)),
  y4_init_s_a = array(y4_init_SA, dim = c(dim.n_subject_2, col_3)),
  y5_init_s_a = array(y5_init_SA, dim = c(dim.n_subject_2, col_3)),
  S_obs_1 = array(s0_1, dim = c(dim.n_subject_1, dim.n_total_1)),
  A_obs_1 = array(a0_1, dim = c(dim.n_subject_1, dim.n_total_1)),
  S_obs_2 = array(s0_2, dim = c(dim.n_subject_2, dim.n_total_2)),
  A_obs_2 = array(a0_2, dim = c(dim.n_subject_2, dim.n_total_2))
)


saveRDS(loss_dat,
        file = paste0(tmp, "/loss_dat.RData"))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
