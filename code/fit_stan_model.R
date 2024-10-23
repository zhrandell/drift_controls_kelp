## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit system of ODEs to empirical data in STAN ~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##


## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

loss_dat <- readRDS(paste0(results,"/loss_dat.Rdata"))

## END start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## run Stan ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## load directory with Stan script

## compile model
model <- cmdstan_model(file <- paste0(code, "/restocking_stan_model.stan"))


## initiate sampling 
fit <- model$sample(data = loss_dat,
                  chains = 4,
                  iter_warmup = 200,
                  iter_sampling = 500,
                  adapt_delta = 0.80, 
                  parallel_chains = 4)


## save via cmdstan's preferred method
fit$save_object(file = paste0(results, "/model_output.RDS"))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## run Stan ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## load directory with Stan script

## compile model
model <- cmdstan_model(file <- paste0(code, "/restocking_stan_model_vL.stan"))


## initiate sampling 
fit <- model$sample(data = loss_dat,
                    chains = 4,
                    iter_warmup = 200,
                    iter_sampling = 500,
                    adapt_delta = 0.80, 
                    parallel_chains = 4)


## save via cmdstan's preferred method
fit$save_object(file = paste0(results, "/model_output_vL.RDS"))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
