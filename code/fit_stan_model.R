## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~ Fit system of ODEs to empirical data in STAN ~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

loss_dat <- readRDS(paste0(results,"/loss_dat.Rdata"))

## run Stan ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## compile model
model <- cmdstan_model(file <- paste0(code,
                                      "/stan_model_",
                                      sel.model,
                                      ".stan"))


## initiate sampling 
fit <- model$sample(data = loss_dat,
                  chains = 4,
                  iter_warmup = 1000,
                  iter_sampling = 3000,
                  adapt_delta = 0.80, 
                  parallel_chains = 4)


## save via cmdstan's preferred method
fit$save_object(file = paste0(results, "/model_output_", 
                              sel.model, ".RDS"))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
