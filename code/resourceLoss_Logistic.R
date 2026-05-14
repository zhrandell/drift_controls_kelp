resourceLoss_Logistic <- function(S0, A0, F0, params) {
  with(as.list(c(S0, A0, F0)), {
  
    # Hunger level
    H = exp(- F)
    # H = 2 / (1 + exp( ( F / z )^s ))
    
    # Movement slowdown
    M = 1 / (1 + b * S^2)
    
    # Logistic preference - additive log formulation [permitting Normal priors on w and q]
    p = ( 1 - ( 1 / ( 1 + exp( w + q * log(S / A) ))))
  
    # Consumption rates
    f_S = S * H * a * p * M        #   / ( 1 + a * b * ( p * S + (1-p) * A )^2 )
    f_A = A * H * a * (1-p) * M    #   / ( 1 + a * b * ( p * S + (1-p) * A )^2 )
  
    # Drift, Kelp, Stomach
  	dS_dt = - U * f_S
  	dA_dt = - U * f_A
  	dF_dt =  (f_S + f_A) - z * F
  	
    return(list(c(dS_dt, dA_dt, dF_dt)))
  })
}
