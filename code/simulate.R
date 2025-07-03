## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Simulate ODE system with posteriors from STAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## ode() with multiple params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## load posts_Df

load(paste0(results,"/posterior_draws_", sel.model,".RDA"))

A.level = c("low" = 30,
            "high" = 300)

for(AL in 1:length(A.level)){ # initial kelp abundance (low and high)

  ## create sequences of initial conditions
  len_init <- 40 # 100
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
  
  
  ## list of params 
  aL <- posts_df_raw$a 
  wL <- posts_df_raw$w
  qL <- posts_df_raw$q
  vL <- posts_df_raw$v
  zL <- posts_df_raw$z
  
  
  ## first set of params 
  parm_list = c(a = aL[1],
                w = wL[1],
                q = qL[1],
                v = vL[1],
                z = zL[1]
                )
  
  
  ## specify ODE functions
  
  resourceLoss_Logistic <- function(S0, A0, F0, params) {
    with(as.list(c(S0, A0, F0)),{
      
      # H = exp(- v * F)
      # H = 1 / (1 + exp( v * (F - z)))
      H = 2 / (1 + exp( ( F / z )^v ))
   
      ## logistic
      p = (1 - (1 / (1 + exp(w + q * log(S / A)))))
      
      f_S = S * H * a * p
      f_A = A * H * a * (1 - p)
      
      dS_dt = - U * f_S
      dA_dt = - U * f_A
      dF_dt = f_S + f_A
      
      return(list(c(dS_dt, dA_dt, dF_dt)))
    })
  }
  
  resourceLoss_vanLeeuwen <- function(S0, A0, F0, params) {
    with(as.list(c(S0, A0, F0)),{
      
      # H = exp(- v * F)
      H = 1 / (1 + exp(v * (F - z)))
      
      
      ## vanLeeuwen
      # We use parameters 'w' and 'q' for convenience though in the notes we use \nu for w and \psi for q
      p = (1 - (1 + exp(w + log(S / A))) / (1 + exp(log(2) + w + log(S / A)) + exp(q + 2 * log(S / A))))
      
      f_S = S * H * a * p
      f_A = A * H * a * (1 - p)
    
      dS_dt = - U * f_S
      dA_dt = - U * f_A
      dF_dt = f_S + f_A
      
      return(list(c(dS_dt, dA_dt, dF_dt)))
    })
  }
  
  if(sel.model == 'Logistic'){
    resourceLoss <- resourceLoss_Logistic
  }
  if(sel.model == 'vanLeeuwen'){
    resourceLoss <- resourceLoss_vanLeeuwen
  }
  
  
  ## run single ODE 
  out_P1 <- ode(init_P1, 
                times = t.list_P1,
                func = resourceLoss,
                parms = parm_list)
  
  
  ## concatenate params into list of lists
  full_parm_list <- mapply(c, 
                           a = aL, w = wL, q = qL, v = vL, z = zL,
                           SIMPLIFY = FALSE)
  
  
  ## t.lists for restocking model
  t.list_P1_restock <- seq(1, P1, by = 1)
  t.list_P2_restock <- seq(P1 + 1, P1 + P2, by = 1)
  t.list_P3_restock <- seq(P1 + P2 + 1, P1 + P2 + P3, by = 1)
  
  
  ## nested lapply 
  outs_parms <- lapply(full_parm_list, function(x){
    lapply(inits_P1, function(y){
      
      p1 <- ode(y,
                times = t.list_P1_restock,
                func = resourceLoss,
                parms = x)
  
      p2 <- ode(c(y[1], y[2], p1[P1, 4]), 
                times = t.list_P2_restock,
                func = resourceLoss,
                parms = x)
  
      p3 <- ode(c(y[1], y[2], p2[P2, 4]), 
                times = t.list_P3_restock,
                func = resourceLoss,
                parms = x)
      
      return(rbind(p1, p2, p3))
  
    })
  })

save(outs_parms, 
     file = paste0(results,"/ODE_kelp_", names(A.level[AL]), '_', sel.model, ".RDA"))

} # end AL (kelp low or high) loop

## END ODE simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


