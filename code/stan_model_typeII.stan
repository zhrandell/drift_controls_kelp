// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Multispecies type II functional response with hunger regulation.
// Equivalent to stan_model_Logistic.stan but with the squaring removed from the
// denominator of M, per Reviewer #1's "Equation Rev1":
//   f_S = h * (p * a * S) / (1 + a * b * (p*S + (1-p)*A))
// Priors / bounds are inherited from the Logistic scaffold; review them once
// the form of M has settled.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

functions {
  vector resourceLoss(real t,              		// time
                      vector Y,            		// state variables          
                      array[] real theta, 		// params
                      vector x_r) {			      // urchins

    real dS_dt;				// Drift rate of change
    real dA_dt;				// Kelp rate of change
    real dF_dt;				// Stomach fullness rate of change
    real f_S;         // Feeding rate on drift
    real f_A;         // Feeding rate on kelp
    real H;           // Hunger level
    real M;           // Movement rate
    real p;           // Preference for drift
    real S = fmax(Y[1], 1e-10); 	// Drift (clamped to prevent log(0) or log(negative))
    real A = fmax(Y[2], 1e-10); 	// Kelp  (clamped to prevent log(0) or log(negative))
    real F = Y[3];			// Stomach fullness 
    real a = theta[1]; 			// baseline attack rate
    real b = theta[2];      // velocity reduction
    real w = theta[3];			// preference par 1
    real q = theta[4]; 		  // preference par 2
    real s = theta[5];			// stomach satiation sensitivity
    real z = theta[6];			// stomach clearance rate 
    real U = x_r[1]; 			  // Urchins
  
  // ## Do Not Remove the ODE_BODY_START and ODE_BODY_END flags ###

  // ## ODE_BODY_START ##

    // Hunger level
    H = exp(- s * F);
    
    // Logistic preference - additive log formulation [permitting Normal priors on w and q]
    p = ( 1 - ( 1 / ( 1 + exp( w + q * log(S / A) ))));
  
    // Movement slowdown (type II: no squaring in denominator)
    M = 1 / ( 1 + a * b * ( p * S + (1-p) * A ) );
    
    // Consumption rates
    f_S = S * H * M * a * p;
    f_A = A * H * M * a * (1-p);
    
    // Drift, Kelp, Stomach
  	dS_dt = - U * f_S;
  	dA_dt = - U * f_A;
  	dF_dt =   f_S + f_A;

  
  // ## ODE_BODY_END ##
  
    return [dS_dt, dA_dt, dF_dt]';
  }
}

data {
  int n_subject_1;					// number of urchin cohorts temp. seq. 1 = 62
  int n_subject_2;					// number of urchin cohorts temp. seq. 2 = 59  
  array[n_subject_1, 3] real y1_init_s_a; 			// Initial condition for S & A for period 1; 		
  array[n_subject_1, 3] real y2_init_s_a; 			// the 1st restocked value for S, A for period 2; 	
  array[n_subject_2, 3] real y3_init_s_a; 			// Initial condition for S & A for period 3; 		
  array[n_subject_2, 3] real y4_init_s_a; 			// the 1st restocked value for S, A for period 4; 	
  array[n_subject_2, 3] real y5_init_s_a; 			// the 2nd restocked value for S, A for period 5;
  int <lower = 1> nts1;               			// # of data collection times in period 1; nts1  =  1		
  int <lower = 1> nts2;               			// # of data collection times in period 2; nts2  =  1		
  int <lower = 1> nts3;               			// # of data collection times in period 3; nts3  =  1 	
  int <lower = 1> nts4;               			// # of data collection times in period 4; nts4  =  1		
  int <lower = 1> nts5;               			// # of data collection times in period 5; nts5  =  1  
  real <lower = 1> t0_1;               			// starting time; t0_1  =  1
  real <lower = 1> t0_2;               			// starting time; t0_2  =  1
  array[nts1] real <lower = t0_1> ts1;          		// data collection time in period 1; ts1  =  24[1]			 
  array[nts2] real <lower = t0_1> ts2;          		// data collection time in period 2; ts2  =  48[1]			
  array[nts3] real <lower = t0_2> ts3;          		// data collection time in period 3; ts3  =  44[1]
  array[nts4] real <lower = t0_2> ts4;          		// data collection time in period 4; ts4  =  89[1]			
  array[nts5] real <lower = t0_2> ts5;          		// data collection time in period 5; ts5  =  134[1]			
  array[n_subject_1, nts1 + nts2] real <lower = 0> S_obs_1; 		     // drift consumed observations period 1
  array[n_subject_2, nts3 + nts4 + nts5] real <lower = 0> S_obs_2;   // drift consumed observations period 1
  array[n_subject_1, nts1 + nts2] real <lower = 0> A_obs_1;          // kelp consumed observations period 2
  array[n_subject_2, nts3 + nts4 + nts5] real <lower = 0> A_obs_2;   // kelp consumed observations period 2
}

transformed data {
  int x_i;
  int x_r;
  int n_total_1;			// n_total_1 = nts1 + nts2 = 2
  int n_total_2;			// n_total_2 = nts3 + nts4 + nts5 = 3
  int N_obs;			// total observations for log_lik / y_rep
  n_total_1 = nts1 + nts2;
  n_total_2 = nts3 + nts4 + nts5;
  N_obs     = 2 * n_subject_1 * n_total_1 + 2 * n_subject_2 * n_total_2;
}

// Narrow down limits to increase sampling efficiency, 
// but keep wide enough to not affect accepted priors
// parameters {
//   real <lower =  0, upper = 0.05> a;
//   real <lower =  0, upper = 0.1>  b;
//   real <lower = -6, upper = 6>    w;
//   real <lower = -5, upper = 5>    q;
//   real <lower =  0, upper = 3>    s;
//   real <lower =  0, upper = 0.1>  z;
//   real <lower = 10, upper = 20>   sigma;
// }

parameters {
  real <lower =  0, upper = 0.05> a;
  real <lower =  0, upper = 0.1>  b;
  real <lower = -6, upper = 6>    w;
  real <lower = -5, upper = 5>    q;
  real <lower =  0, upper = 0.5>  s;
  real <lower =  0, upper = 0.1>  z;
  real <lower = 10, upper = 20>   sigma;
}



transformed parameters {
  array[6] real theta;
  array[nts1] vector[3] y1;					// two-dimensional container of size (nts1, 3) i.e. y1[1, 3] 
  array[nts2] vector[3] y2;					// two-dimensional container of size (nts2, 3)   
  array[nts3] vector[3] y3;					// two-dimensional container of size (nts3, 3)   
  array[nts4] vector[3] y4;					// two-dimensional container of size (nts4, 3)     
  array[nts5] vector[3] y5;					// two-dimensional container of size (nts5, 3)  
  vector[3] init_1;
  vector[3] init_2;
  vector[1] U;			
  
  // obs 
  array[n_subject_1, n_total_1] real drift_1; 	// Drift remaining
  array[n_subject_1, n_total_1] real kelp_1; 		// Kelp remaining
  array[n_subject_2, n_total_2] real drift_2; 	// Drift remaining
  array[n_subject_2, n_total_2] real kelp_2; 		// Kelp remaining
  
  theta[1] = a; 
  theta[2] = b;
  theta[3] = w;
  theta[4] = q;
  theta[5] = s;
  theta[6] = z;

  // Temporal sequence 2 -----------------------------------
  for (i in 1:n_subject_1) {

    // period 1
    init_1[1] = y1_init_s_a[i, 1]; // initial drift
    init_1[2] = y1_init_s_a[i, 2]; // initial kelp
    init_1[3] = 0;                 // initial urchin gut fullness
    U[1] = y1_init_s_a[i, 3];      // initial urchin
    y1 = ode_rk45(resourceLoss, init_1, t0_1, ts1, theta, U);
    drift_1[i, 1:nts1] = y1[, 1];
    kelp_1[i, 1:nts1] = y1[, 2];


    // period 2
    init_1[1] = y2_init_s_a[i, 1];
    init_1[2] = y2_init_s_a[i, 2];
    init_1[3] = y1[nts1, 3];
    U[1] = y2_init_s_a[i, 3];
    y2 = ode_rk45(resourceLoss, init_1, ts1[nts1], ts2, theta, U);
    drift_1[i, (nts1 + 1):(nts1 + nts2)] = y2[, 1];
    kelp_1[i, (nts1 + 1):(nts1 + nts2)] = y2[, 2];
    
    // print(i, " prm[a,b,w,q,s]: ", theta);
    // print(i, " y1: ", y1);
    // print(i, " y2: ", y2);
    // print(i, " U: ", U);
    // print(i, " dXdt: ", resourceLoss(t0_1, init_1, theta, U));
    
    // print(i, " drift_1: ", drift_1[i, ]);
    // print(i, " kelp_1: ", kelp_1[i, ]);
    // print(i, " sigma: ", sigma);
  }


  // Temporal sequence 1 -----------------------------------
  for (i in 1:n_subject_2) {      
   
    // period 3
    init_2[1] = y3_init_s_a[i, 1];
    init_2[2] = y3_init_s_a[i, 2];
    init_2[3] = 0;
    U[1] = y3_init_s_a[i, 3];
    y3 = ode_rk45(resourceLoss, init_2, t0_2, ts3, theta, U);
    drift_2[i, 1:nts3] = y3[, 1];
    kelp_2[i, 1:nts3] = y3[, 2];

    // period 4
    init_2[1] = y4_init_s_a[i, 1];
    init_2[2] = y4_init_s_a[i, 2];
    init_2[3] = y3[nts3, 3];
    U[1] = y4_init_s_a[i, 3];
    y4 = ode_rk45(resourceLoss, init_2, ts3[nts3], ts4, theta, U);
    drift_2[i, (nts3 + 1):(nts3 + nts4)] = y4[, 1];
    kelp_2[i, (nts3 + 1):(nts3 + nts4)] = y4[, 2];

    // period 5
    init_2[1] = y5_init_s_a[i, 1];
    init_2[2] = y5_init_s_a[i, 2];
    init_2[3] = y4[nts4, 3];
    U[1] = y5_init_s_a[i, 3];
    y5 = ode_rk45(resourceLoss, init_2, ts4[nts4], ts5, theta, U); 	
    drift_2[i, (nts3 + nts4 + 1):(nts3 + nts4 + nts5)] = y5[, 1];
    kelp_2[i, (nts3 + nts4 + 1):(nts3 + nts4 + nts5)] = y5[, 2];

    // print(i, " init_2: ", init_2);
    // print(i, " prm[a,v,p,w,q]: ", theta);
    // print(i, " y3: ", y3);
    // print(i, " y4: ", y4);
    // print(i, " y5: ", y5);
    // print(i, " U: ", U);
    // print(i, " dXdt: ", resourceLoss(t0_1, init_2, theta, U));
  
    // print(i, " drift_2: ", drift_2[i, ]);
    // print(i, " kelp_2: ", kelp_2[i, ]);
    // print(i, " sigma: ", sigma);
  }
}

model {
  // a ~ exponential(10);
  // b ~ exponential(10);
  // // b ~ normal(0, 1);
  // w ~ normal(0, 1.8); // normal(0, 1.8) is ~uniform on logistic scale
  // q ~ normal(0, 10);
  // s ~ exponential(1);
  // // s ~ normal(0, 1);
  // z ~ exponential(10);
  // sigma ~ exponential(0.1);
  
  a ~ exponential(10);
  b ~ exponential(10);
  w ~ normal(0, 1.8); // normal(0, 1.8) is ~uniform on logistic scale
  q ~ normal(0, 10);
  s ~ exponential(1);
  z ~ exponential(10);
  sigma ~ exponential(0.1);


  for (i in 1:n_subject_1) {
    if(y1_init_s_a[i, 1] > 0){
      target += normal_lpdf(S_obs_1[i, ] | drift_1[i, ], sigma);
    }
    if(y1_init_s_a[i, 2] > 0){
      target += normal_lpdf(A_obs_1[i, ] | kelp_1[i, ], sigma);
    }
  }

  for (i in 1:n_subject_2) {
    if(y3_init_s_a[i, 1] > 0){
      target += normal_lpdf(S_obs_2[i, ] | drift_2[i, ], sigma);
    }
    if(y3_init_s_a[i, 2] > 0){
      target += normal_lpdf(A_obs_2[i, ] | kelp_2[i, ], sigma);
    }
  }
}

// Per-observation log_lik for PSIS-LOO (loo::loo) and y_rep for posterior predictive checks (bayesplot::ppc_*).
// log_lik / y_rep ordering: drift_1 (subjects x periods row-major), kelp_1, drift_2, kelp_2.
generated quantities {
  vector[N_obs] log_lik;
  array[N_obs] real y_rep;
  {
    int idx = 1;
    for (i in 1:n_subject_1) {
      for (t in 1:n_total_1) {
        log_lik[idx] = normal_lpdf(S_obs_1[i, t] | drift_1[i, t], sigma);
        y_rep[idx]   = normal_rng(drift_1[i, t], sigma);
        idx += 1;
      }
    }
    for (i in 1:n_subject_1) {
      for (t in 1:n_total_1) {
        log_lik[idx] = normal_lpdf(A_obs_1[i, t] | kelp_1[i, t], sigma);
        y_rep[idx]   = normal_rng(kelp_1[i, t], sigma);
        idx += 1;
      }
    }
    for (i in 1:n_subject_2) {
      for (t in 1:n_total_2) {
        log_lik[idx] = normal_lpdf(S_obs_2[i, t] | drift_2[i, t], sigma);
        y_rep[idx]   = normal_rng(drift_2[i, t], sigma);
        idx += 1;
      }
    }
    for (i in 1:n_subject_2) {
      for (t in 1:n_total_2) {
        log_lik[idx] = normal_lpdf(A_obs_2[i, t] | kelp_2[i, t], sigma);
        y_rep[idx]   = normal_rng(kelp_2[i, t], sigma);
        idx += 1;
      }
    }
  }
}

// End of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

