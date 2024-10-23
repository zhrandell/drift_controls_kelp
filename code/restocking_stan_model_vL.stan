// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ODE system for Drift loss, Kelp loss, and Urchin stomach fullness
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
      real S = Y[1]; 			// Drift
      real A = Y[2]; 			// Kelp
      real F = Y[3];			// Stomach fullness 
      real a = theta[1]; 			// baseline attack rate
      real v = theta[2];			// stomach satiation sensitivity
      real w = theta[3];			// relative preference
      real q = theta[4]; 		  // switching rate
      real U = x_r[1]; 			  // Urchins
      
      // For control treatments
      if(S == 0){
        f_S = 0;
      }
      else{
        f_S = S * a;
      }
      if(A == 0){
        f_A = 0;
      }
      else{
        f_A = A * a;
      }
      if(S > 0 && A > 0){
        // Yodzis preference formulation [requiring constrained 0-1 prior on w]
        //   f_S = S * a * (( w  *   pow(S, q)) / ( (w * pow(S, q)) + ((1-w) * pow(A, q)) ));
        //   f_A = A * a * (((1-w) * pow(A, q)) / ( (w * pow(S, q)) + ((1-w) * pow(A, q)) ));
        
        // Logistic preference formulation [requiring constrained 0-1 prior on w]
        // f_S = S * a * ( 1 - ( 1 / ( 1 + (w / (1 - w)) * pow((S / A), q)) ));
        // f_A = A * a * (     ( 1 / ( 1 + (w / (1 - w)) * pow((S / A), q)) ));
        
        // Logistic preference - multiplicative log formulation [permitting Normal prior on w]
        // f_S = S * a * ( 1 - ( 1 / ( 1 + exp(w) * pow((S / A), q)) ));
        // f_A = A * a * (     ( 1 / ( 1 + exp(w) * pow((S / A), q)) ));
        
        // Logistic preference - additive log formulation [permitting Normal priors on w and q]
        // f_S = S * a * ( 1 - ( 1 / ( 1 + exp( w + q * log(S / A) ))));
        // f_A = A * a * (     ( 1 / ( 1 + exp( w + q * log(S / A) ))));
        
        // vanLeeuwen et al.
        // f_S = S * a * ( 1 - ( 1 + w * (1 / q) * (S / A) ) / ( 1 + w * (1 / q) * (S / A) * 2 + ( w * (S / A) )^2 ) );
        // f_A = A * a * (     ( 1 + w * (1 / q) * (S / A) ) / ( 1 + w * (1 / q) * (S / A) * 2 + ( w * (S / A) )^2 ) )

        // vanLeeuwen et al. reformulated
        f_S = S * a * ( 1 -  ( 1 + exp( w + log(S / A) )) / ( 1 + exp( log(2) + w + log(S / A) ) + exp( q + 2 * log(S / A) )) );
        f_A = A * a * (      ( 1 + exp( w + log(S / A) )) / ( 1 + exp( log(2) + w + log(S / A) ) + exp( q + 2 * log(S / A) )) );      
      }
      
      // Hunger level
      H = exp(- v * F);
      // H = 1 - v * F;
      
      // Drift, Kelp, Stomach
      dS_dt = - U * f_S * H;
      dA_dt = - U * f_A * H;
      dF_dt =   f_S + f_A;
      
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
  n_total_1 = nts1 + nts2;
  n_total_2 = nts3 + nts4 + nts5;
}

// Narrow down limits to increase sampling efficiency, 
// but keep wide enough to not affect accepted priors
parameters {
  real <lower = 0, upper = 0.01> a;
  real <lower = 0, upper = 0.3> v;
  real <lower = -20, upper = 20> w;
  real <lower = 0, upper = 10> q;
  real <lower = 10, upper = 30> sigma;
}

transformed parameters {
  array[4] real theta;
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
  theta[2] = v;
  theta[3] = w;
  theta[4] = q;


  // Temporal sequence 1 -----------------------------------
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
    
    // print(i, " init_1: ", init_1);
    // print(i, " prm[a,v,w,q]: ", theta);
    // print(i, " y1: ", y1);
    // print(i, " y2: ", y2);
    // print(i, " U: ", U);
    // print(i, " dXdt: ", resourceLoss(t0_1, init_1, theta, U));
    
    // print(i, " drift_1: ", drift_1[i, ]);
    // print(i, " kelp_1: ", kelp_1[i, ]);
    // print(i, " sigma: ", sigma);
  }


  // Temporal sequence 2 -----------------------------------
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
  a ~ exponential(1);
  v ~ exponential(1);
  w ~ normal(0, 10);
  q ~ normal(0, 10);
  sigma ~ exponential(0.1);

  for (i in 1:n_subject_1) {
    if(y1_init_s_a[i, 1] > 0){
      target += normal_lpdf(S_obs_1[i, ] | drift_1[i, ], sigma);
    }
    if(y1_init_s_a[i, 2] > 0){
      target += normal_lpdf(S_obs_1[i, ] | drift_1[i, ], sigma);
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

// End of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

