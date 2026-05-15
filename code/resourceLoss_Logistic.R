resourceLoss_Logistic <- function(S0, A0, F0, params) {
  with(as.list(c(S0, A0, F0)), {

    # Hunger level
    H = exp(- s * F)
    
    # Logistic preference - additive log formulation [permitting Normal priors on w and q]
    p = ( 1 - ( 1 / ( 1 + exp( w + q * log(S / A) ))))
  
    # Movement slowdown
    M = 1 / ( 1 + a * b * ( p * S + (1-p) * A )^2 )
    
    # Consumption rates
    f_S = S * H * M * a * p
    f_A = A * H * M * a * (1-p)
    
    # Drift, Kelp, Stomach
  	dS_dt = - U * f_S
  	dA_dt = - U * f_A
  	dF_dt =   f_S + f_A

  
    return(list(c(dS_dt, dA_dt, dF_dt)))
  })
}
