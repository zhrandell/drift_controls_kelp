## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Simulate ODE system with posteriors from STAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# Initiate clusters for parallel implementation
cat('Initiating parallel simulation of model.\n')
cl <- makeCluster(n_cores)
invisible(clusterEvalQ(cl, library(deSolve)))


## ode() with multiple params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## load posts_Df

load(paste0(results,"/posterior_draws_", sel.model,".RDA"))

## generate resourceLoss_*.R from the bracketed ODE body in the Stan file
make_resourceLoss <- function(sel.model, code_dir) {
  stan_lines <- readLines(paste0(code_dir, "/stan_model_", sel.model, ".stan"))

  start <- which(grepl("## ODE_BODY_START ##", stan_lines)) + 1
  end   <- which(grepl("## ODE_BODY_END ##",   stan_lines)) - 1

  body <- stan_lines[start:end]
  body <- gsub("//", "#",  body)                    # Stan comments -> R comments
  body <- gsub(";(\\s*#|\\s*$)", "\\1", body)       # remove semicolons before # or at end of line

  writeLines(c(
    paste0("resourceLoss_", sel.model, " <- function(S0, A0, F0, params) {"),
    "  with(as.list(c(S0, A0, F0)), {",
    body,
    "    return(list(c(dS_dt, dA_dt, dF_dt)))",
    "  })",
    "}"
  ), paste0(code_dir, "/resourceLoss_", sel.model, ".R"))
}

make_resourceLoss(sel.model, code)
source(paste0(code, "/resourceLoss_", sel.model, ".R"))
resourceLoss <- get(paste0("resourceLoss_", sel.model))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(AL in 1:length(A.level)){ # initial kelp abundance (low and high)

  ## create sequences of initial conditions
  len_init <- 80 # 100
  A0 <- seq(A.level[AL], A.level[AL], length.out = len_init)    # Kelp
  S0 <- seq(  1, 300, length.out = len_init)    # Drift
  F0 <- seq(  0,   0, length.out = len_init)    # Stomach Fullness
  U <- 20                                       # Urchins
  
  
  ## set time points (in hrs) for the three Periods
  P1 <- 44 # 44 hrs
  P2 <- 45 # 89 hrs
  P3 <- 45 # 134 hrs
  
  
  ## time sequences to pass to ode()
  t.list_P1 <- seq(1, P1, by = 1)
  t.list_P2 <- seq(1, P2, by = 1)
  t.list_P3 <- seq(1, P3, by = 1)
  
  
  ## set initial conditions 
  init_P1 = c(S = S0[1], 
              A = A0[1], 
              F = F0[1])
  
  
  ## concatenate into list of lists (mapply() doesn't like 'F', so we rename after)
  inits_P1 <- mapply(c, 
                     S = S0, A = A0, ff = F0, 
                     SIMPLIFY = FALSE)
  inits_P1 <- lapply(inits_P1, function(x){ names(x) <- sub("ff", "F", names(x)); x})
  

  ## derive ODE parameter columns automatically from posts_df_raw
  non_ode_parms <- c("sigma")   # Stan params not used in the ODE
  ode_parm_cols <- setdiff(
    names(posts_df_raw)[!startsWith(names(posts_df_raw), ".")],
    non_ode_parms
  )

  ## first set of params (single draw, for test ODE run)
  parm_list <- setNames(as.numeric(posts_df_raw[1, ode_parm_cols]), ode_parm_cols)


  ## run single ODE 
  out_P1 <- ode(init_P1, 
                times = t.list_P1,
                func = resourceLoss,
                parms = parm_list)
  

  ## full list of params (all draws, for parallel run)
  full_parm_list <- lapply(seq_len(nrow(posts_df_raw)), function(i)
    setNames(as.numeric(posts_df_raw[i, ode_parm_cols]), ode_parm_cols)
  )


  ## t.lists for restocking model
  t.list_P1_restock <- seq(1, P1, by = 1)
  t.list_P2_restock <- seq(P1 + 1, P1 + P2, by = 1)
  t.list_P3_restock <- seq(P1 + P2 + 1, P1 + P2 + P3, by = 1)
  
  
  ## flatten parameter x initial-condition grid for better parallel utilisation
  print(paste('Kelp level', AL, 'of', length(A.level)))
  n_parms <- length(full_parm_list)
  n_inits <- length(inits_P1)
  combos  <- expand.grid(init_idx = seq_len(n_inits), parm_idx = seq_len(n_parms))

  clusterExport(cl,
    c("full_parm_list", "inits_P1", "combos",
      "t.list_P1_restock", "t.list_P2_restock", "t.list_P3_restock",
      "resourceLoss", "P1", "P2", "U"),
    envir = environment())

  flat_results <- pblapply(seq_len(nrow(combos)), function(i) {
    x <- full_parm_list[[combos$parm_idx[i]]]
    y <- inits_P1[[combos$init_idx[i]]]

    p1 <- ode(y,
              times = t.list_P1_restock,
              func  = resourceLoss,
              parms = x)

    p2 <- ode(c(y[1], y[2], p1[P1, 4]),
              times = t.list_P2_restock,
              func  = resourceLoss,
              parms = x)

    p3 <- ode(c(y[1], y[2], p2[P2, 4]),
              times = t.list_P3_restock,
              func  = resourceLoss,
              parms = x)

    rbind(p1, p2, p3)
  }, cl = cl)

  ## re-nest into original outs_parms[[parm_idx]][[init_idx]] structure
  outs_parms <- split(flat_results, combos$parm_idx)

save(outs_parms, 
     file = paste0(results,"/ODE_kelp_", names(A.level[AL]), '_', sel.model, ".RDA"))

} # end AL (kelp low or high) loop

## END ODE simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stopCluster(cl)
cat('Clusters closed.\n')