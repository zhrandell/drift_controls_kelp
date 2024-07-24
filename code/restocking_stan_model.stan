// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ODE system for Drift loss, Kelp loss, and Urchin stomach fullness

// created April 24th, 2021; updated July 17th, 2024
// zhr
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Originally created April 24th, 2021; updated October 7th, 2021
// zhr 
// Updated using StanGTP (via ChatGTP) on July 20th, 2024 by MN
// StanGTP-stated key changes:
// 1. Changed `real[] theta` to `array[] real theta`.
// 2. Updated array declarations in the `data`, `transformed data`, and `transformed parameters` blocks to the new array syntax `array[dimensions] real`.
// 3. Ensured all vectors and arrays are declared with the correct new syntax.
// Additional changes and Improvements:
// 1. **Consistent Naming**: Ensured consistent naming conventions and alignment.
// 2. **Correct `array` Syntax**: Used `array` syntax appropriately for multi-dimensional arrays.
// 3. **General Code Cleanliness**: Improved readability and organization of the code.

functions {
  vector resourceLoss(real t,              		// time
                      vector Y,            		// state variables          
                      array[] real theta, 		// params
                      vector x_r) {			      // urchins

    real dS_dt;				// Drift remaining through time
    real dA_dt;				// Kelp remaining through time
    real dF_dt;				// Stomach fullness through time
    real S = Y[1]; 			// Drift
    real A = Y[2]; 			// Kelp
    real F = Y[3];			// Stomach fullness 
    real a = theta[1]; 			// encounter rate
    real v = theta[2];			// max Stomach volume
    real q = theta[3];			// switching param for Low data
    real p = theta[4];			// stomach clearance
    real w = theta[5]; 		  // baseline preference for drift over kelp
    real U = x_r[1]; 			  // Urchins

// Yodzis formulation
// 	dS_dt = - U * a * S * (1 - F / v) * (( w  *   pow(S, q)) / ( (w * pow(S, q)) + ((1-w) * pow(A, q)) ));
// 	dA_dt = - U * a * A * (1 - F / v) * (((1-w) * pow(A, q)) / ( (w * pow(S, q)) + ((1-w) * pow(A, q)) ));
// 	dF_dt =   a * S * (1 - F / v) * (( w  *   pow(S, q)) / ( (w * pow(S, q)) + ((1-w) * pow(A, q)) ))
// 	        + a * A * (1 - F / v) * (((1-w) * pow(A, q)) / ( (w * pow(S, q)) + ((1-w) * pow(A, q)) ))
//           - p * F;
          
// Logistic formulation
	dS_dt = - U * a * S * (1 - F / v) * ( 1 - ( 1 / ( 1 + (w / (1 - w)) * pow((S / A), q)) ));
	dA_dt = - U * a * A * (1 - F / v) * (     ( 1 / ( 1 + (w / (1 - w)) * pow((S / A), q)) ));
	dF_dt =   a * S * (1 - F / v) * ( 1 - ( 1 / ( 1 + (w / (1 - w)) * pow((S / A), q)) ))
	        + a * A * (1 - F / v) * (     ( 1 / ( 1 + (w / (1 - w)) * pow((S / A), q)) ))
          - p * F;

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
  int <lower=1> nts1;               			// # of data collection times in period 1; nts1 = 1		
  int <lower=1> nts2;               			// # of data collection times in period 2; nts2 = 1		
  int <lower=1> nts3;               			// # of data collection times in period 3; nts3 = 1 	
  int <lower=1> nts4;               			// # of data collection times in period 4; nts4 = 1		
  int <lower=1> nts5;               			// # of data collection times in period 5; nts5 = 1  
  real <lower=1> t0_1;               			// starting time; t0_1 = 1
  real <lower=1> t0_2;               			// starting time; t0_2 = 1
  array[nts1] real <lower=t0_1> ts1;          		// data collection time in period 1; ts1 = 24[1]			 
  array[nts2] real <lower=t0_1> ts2;          		// data collection time in period 2; ts2 = 48[1]			
  array[nts3] real <lower=t0_2> ts3;          		// data collection time in period 3; ts3 = 44[1]
  array[nts4] real <lower=t0_2> ts4;          		// data collection time in period 4; ts4 = 89[1]			
  array[nts5] real <lower=t0_2> ts5;          		// data collection time in period 5; ts5 = 134[1]			
  array[n_subject_1, nts1 + nts2] real S_obs_1;            	// experimental observations; [59, nts1 + nts2]  	       					
  array[n_subject_2, nts3 + nts4 + nts5] real S_obs_2;      // experimental drift consumed observations; [59, nts3 + nts4 + nts5]  
  array[n_subject_1, nts1 + nts2] real A_obs_1;            	// experimental kelp consumed observations; [59, nts1 + nts2]
  array[n_subject_2, nts3 + nts4 + nts5] real A_obs_2;      // experimental kelp consumed observations; [59, nts3 + nts4 + nts5]
}

transformed data {
  int x_i; 				 
  int x_r;              
  int n_total_1;			// n_total_1 = nts1 + nts2 = 2
  int n_total_2;			// n_total_2 = nts3 + nts4 + nts5 = 3	
  n_total_1 = nts1 + nts2;
  n_total_2 = nts3 + nts4 + nts5;
}

parameters {
   real <lower=0, upper=.1> a; 
   real <lower=0, upper=10> v;
   real <lower=0, upper=10> q;
   real <lower=0, upper=1> p;
   real <lower=0, upper=1> w;
   real <lower=0, upper=40> sigma;      
}

transformed parameters {
  array[5] real theta;
  array[nts1] vector[3] y1;					// two-dimensional container of size (nts1, 3) i.e. y1[1, 3] 
  array[nts2] vector[3] y2;					// two-dimensional container of size (nts2, 3)   
  array[nts3] vector[3] y3;					// two-dimensional container of size (nts3, 3)   
  array[nts4] vector[3] y4;					// two-dimensional container of size (nts4, 3)     
  array[nts5] vector[3] y5;					// two-dimensional container of size (nts5, 3)  
  vector[3] init_1;
  vector[3] init_2;
  vector[1] U;			
  
  // obs 
  array[n_subject_1, n_total_1] real drift_loss_1; 		// Drift consumed; [59, n_total_1]
  array[n_subject_1, n_total_1] real kelp_loss_1; 		// Kelp consumed; [59, n_total_1]
  array[n_subject_2, n_total_2] real drift_loss_2; 		// Drift consumed; [59, n_total_2]
  array[n_subject_2, n_total_2] real kelp_loss_2; 		// Kelp consumed; [59, n_total_2]
    
  // drift alpha, beta 
  array[n_subject_1, n_total_1] real alphaS_1;  // reparameterized shape for gamma distribution; drift likelihood; [59, n_total_1] 
  array[n_subject_1, n_total_1] real betaS_1;		// reparameterized scale for gamma distribution; drift likelihood; [59, n_total_1]
  array[n_subject_2, n_total_2] real alphaS_2;  // reparameterized shape for gamma distribution; drift likelihood; [59, n_total_2] 
  array[n_subject_2, n_total_2] real betaS_2;  

  // kelp alpha, beta
  array[n_subject_1, n_total_1] real alphaA_1;  // reparameterized shape for gamma distribution; kelp likelihood; [59, n_total_1] 
  array[n_subject_1, n_total_1] real betaA_1;		// reparameterized scale for gamma distribution; kelp likelihood; [59, n_total_1]
  array[n_subject_2, n_total_2] real alphaA_2;  // reparameterized shape for gamma distribution; kelp likelihood; [59, n_total_2] 
  array[n_subject_2, n_total_2] real betaA_2;
  
  theta[1] = a; 
  theta[2] = v;
  theta[3] = q;
  theta[4] = p;
  theta[5] = w;


  // Temporal sequence 1 -----------------------------------
  for (i in 1:n_subject_1) {

    // period 1
    init_1[1] = y1_init_s_a[i, 1]; // initial drift
    init_1[2] = y1_init_s_a[i, 2]; // initial kelp
    init_1[3] = 0;                 // initial urchin gut fullness
    U[1] = y1_init_s_a[i, 3];      // initial urchin
    y1 = ode_rk45(resourceLoss, init_1, t0_1, ts1, theta, U);
    drift_loss_1[i, 1:nts1] = y1[, 1];
    kelp_loss_1[i, 1:nts1] = y1[, 2];
	
    // period 2
    init_1[1] = y2_init_s_a[i, 1]; 		
    init_1[2] = y2_init_s_a[i, 2];
    init_1[3] = y1[nts1, 3];
    U[1] = y2_init_s_a[i, 3];
    y2 = ode_rk45(resourceLoss, init_1, ts1[nts1], ts2, theta, U);
    drift_loss_1[i, (nts1+1):(nts1+nts2)] = y2[, 1];
    kelp_loss_1[i, (nts1+1):(nts1+nts2)] = y2[, 2];
    
    // print(i, " init_1: ", init_1);
    // print(i, " prm[a,v,q,p]: ", theta);
    // print(i, " y1: ", y1);
    // print(i, " y2: ", y2);
    // print(i, " U: ", U);
    // print(i, " dXdt: ", resourceLoss(t0_1, init_1, theta, U));
    
    // print(i, " drift_loss_1: ", drift_loss_1[i, ]);
    // print(i, " kelp_loss_1: ", kelp_loss_1[i, ]);
    // print(i, " sigma: ", sigma);

    for (j in 1:n_total_1) {
      alphaS_1[i, j] = pow(drift_loss_1[i, j], 2) / pow(sigma, 2);
      betaS_1[i, j] = (1 / (pow(sigma, 2) / drift_loss_1[i, j]));
      alphaA_1[i, j] = pow(kelp_loss_1[i, j], 2) / pow(sigma, 2);
      betaA_1[i, j] = (1 / (pow(sigma, 2) / kelp_loss_1[i, j]));            
    }
    
    // print(i, " alphaS_1: ", alphaS_1[i, ]);
    // print(i, " betaS_1: ", betaS_1[i, ]);
    // print(i, " alphaA_1: ", alphaA_1[i, ]);
    // print(i, " betaA_1: ", betaA_1[i, ]);

  }


  // Temporal sequence 2 -----------------------------------
  for (i in 1:59) {      
   
    // period 3
    init_2[1] = y3_init_s_a[i, 1];
    init_2[2] = y3_init_s_a[i, 2];
    init_2[3] = 0;
    U[1] = y3_init_s_a[i, 3];
    y3 = ode_rk45(resourceLoss, init_2, t0_2, ts3, theta, U);
    drift_loss_2[i, 1:nts3] = y3[, 1];
    kelp_loss_2[i, 1:nts3] = y3[, 2];
	
    // period 4
    init_2[1] = y4_init_s_a[i, 1]; 		
    init_2[2] = y4_init_s_a[i, 2];
    init_2[3] = y3[nts3, 3];
    U[1] = y4_init_s_a[i, 3];
    y4 = ode_rk45(resourceLoss, init_2, ts3[nts3], ts4, theta, U);
    drift_loss_2[i, (nts3+1):(nts3+nts4)] = y4[, 1];
    kelp_loss_2[i, (nts3+1):(nts3+nts4)] = y4[, 2];
    
    // period 5
    init_2[1] = y5_init_s_a[i, 1];
    init_2[2] = y5_init_s_a[i, 2];
    init_2[3] = y4[nts4, 3];
    U[1] = y5_init_s_a[i, 3];
    y5 = ode_rk45(resourceLoss, init_2, ts4[nts4], ts5, theta, U); 	
    drift_loss_2[i, (nts3+nts4+1):(nts3+nts4+nts5)] = y5[, 1];
    kelp_loss_2[i, (nts3+nts4+1):(nts3+nts4+nts5)] = y5[, 2];

    // print(i, " init_2: ", init_2);
    // print(i, " prm[a,v,q,p]: ", theta);
    // print(i, " y3: ", y3);
    // print(i, " y4: ", y4);
    // print(i, " y5: ", y5);
    // print(i, " U: ", U);
    // print(i, " dXdt: ", resourceLoss(t0_1, init_2, theta, U));
  
    // print(i, " drift_loss_2: ", drift_loss_2[i, ]);
    // print(i, " kelp_loss_2: ", kelp_loss_2[i, ]);
    // print(i, " sigma: ", sigma);
    
    for (j in 1:n_total_2) {
      alphaS_2[i, j] = pow(drift_loss_2[i, j], 2) / pow(sigma, 2);
      betaS_2[i, j] = (1 / (pow(sigma, 2) / drift_loss_2[i, j]));
      alphaA_2[i, j] = pow(kelp_loss_2[i, j], 2) / pow(sigma, 2);
      betaA_2[i, j] = (1 / (pow(sigma, 2) / kelp_loss_2[i, j]));            
    }
    
    // print(i, " alphaS_2: ", alphaS_2[i, ]);
    // print(i, " betaS_2: ", betaS_2[i, ]);
    // print(i, " alphaA_2: ", alphaA_2[i, ]);
    // print(i, " betaA_2: ", betaA_2[i, ]);
    
  }
  
    // print(" alphaS_1: ", alphaS_1);
    // print(" betaS_1: ", betaS_1);
    // print(" alphaA_1: ", alphaA_1);
    // print(" betaA_1: ", betaA_1);
    
    // print(" alphaS_2: ", alphaS_2);
    // print(" betaS_2: ", betaS_2);
    // print(" alphaA_2: ", alphaA_2);
    // print(" betaA_2: ", betaA_2);
    
}

model { 
  a ~ exponential(0.1);                
  v ~ exponential(0.1); 
  q ~ lognormal(1, 1);
  p ~ beta(1,1);
  w ~ uniform(0,1);
  sigma ~ exponential(0.1);

  for (i in 1:n_subject_1) { 
    S_obs_1[i, ] ~ gamma(alphaS_1[i, ], betaS_1[i, ]);
    A_obs_1[i, ] ~ gamma(alphaA_1[i, ], betaA_1[i, ]);
  }

  for (i in 1:n_subject_2) {
    S_obs_2[i, ] ~ gamma(alphaS_2[i, ], betaS_2[i, ]);
    A_obs_2[i, ] ~ gamma(alphaA_2[i, ], betaA_2[i, ]);
  }
}

// End of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

