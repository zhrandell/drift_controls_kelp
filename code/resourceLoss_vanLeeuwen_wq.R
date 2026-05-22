resourceLoss_vanLeeuwen_wq <- function(S0, A0, F0, params) {
  with(as.list(c(S0, A0, F0)), {
      
      # Hunger level
      H = exp( - s * F )
      # H = 2 / (1 + exp( ( F / z )^s ))
      
      # vanLeeuwen et al. reformulated
      # We use parameters 'w = q-4' and 'q' for convenience though in the notes we use \nu for w and \psi for q
      p = ( 1 -  ( 1 + exp( (q-4) + log(S / A) )) / ( 1 + exp( log(2) + (q-4) + log(S / A) ) + exp( q + 2 * log(S / A) )) )
      
      # Movement slowdown 
      M = 1 / ( 1 + a * b * ( p * S + (1-p) * A )^2 )
      
      # Consumption rates
      f_S = S * H * M * a * p
      f_A = A * H * M * a * (1 - p)

      # Drift, Kelp, Stomach
      dS_dt = - U * f_S
      dA_dt = - U * f_A
      dF_dt =   f_S + f_A
      
    return(list(c(dS_dt, dA_dt, dF_dt)))
  })
}
