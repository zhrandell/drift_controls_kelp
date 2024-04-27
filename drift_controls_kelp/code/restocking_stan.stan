// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ODE system for Drift loss, Kelp loss, and Urchin stomach fullness
// created April 24th, 2021; updated October 7th, 2021
// zhr 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
functions {
  vector resourceLoss(real t,              		// time
                      vector Y,            		// state variables          
                      real[] theta, 			// params
		      vector x_r) {			// urchin #

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
    //real e = theta[5];
    real U = x_r[1]; 			// urchin #

	dS_dt = (-(a * U * S) * (1 - (1 / pow (1 + S/A, q))) * ((v - F) / v));
	dA_dt = (-(a * U * A) * (1 / pow (1 + S/A, q)) * ((v - F) / v));
	dF_dt = ((((a * S) * (1 - (1 / pow (1 + S/A, q))) * ((v - F) / v)) + 
        	((a * A) * (1 / pow (1 + S/A, q)) * ((v - F) / v))) - (p * F));

    return [dS_dt, dA_dt, dF_dt]';
  }
}


data {
  int n_subject_1;					// number of urchin cohorts temp. seq. 1 = 62
  int n_subject_2;					// number of urchin cohorts temp. seq. 2 = 59  
  real y1_init_s_a[n_subject_1, 3]; 			// Initial condiiton for S & A for period 1; 		
  real y2_init_s_a[n_subject_1, 3]; 			// the 1st restocked value for S, A for period 2; 	
  real y3_init_s_a[n_subject_2, 3]; 			// Initial condiiton for S & A for period 3; 		
  real y4_init_s_a[n_subject_2, 3]; 			// the 1st restocked value for S, A for period 4; 	
  real y5_init_s_a[n_subject_2, 3]; 			// the 2nd restocked value for S, A for period 5;
  int <lower=1> nts1;               			// # of data collection times in period 1; nts1 = 1		
  int <lower=1> nts2;               			// # of data collection times in period 2; nts2 = 1		
  int <lower=1> nts3;               			// # of data collection times in period 3; nts3 = 1 	
  int <lower=1> nts4;               			// # of data collection times in period 2; nts2 = 1		
  int <lower=1> nts5;               			// # of data collection times in period 3; nts3 = 1  
  real <lower=1> t0_1;               			// starting time; t0 = 1
  real <lower=1> t0_2;               			// starting time; t0 = 1
  real <lower=t0_1> ts1[nts1];          		// data collection time in period 1; ts1 = 24[1]			 
  real <lower=t0_1> ts2[nts2];          		// data collection time in period 2; ts2 = 48[1]			
  real <lower=t0_2> ts3[nts3];          		// data collection time in period 3; ts3 = 44[1]
  real <lower=t0_2> ts4[nts4];          		// data collection time in period 4; ts5 = 89[1]			
  real <lower=t0_2> ts5[nts5];          		// data collection time in period 5; ts6 = 134[1]			
  real S_obs_1[n_subject_1, nts1 + nts2];            	// experimental observations; [59, 3]  	       					
  real S_obs_2[n_subject_2, nts3 + nts4 + nts5];        // experimental drift consumed observations; [59, 3]  
  real A_obs_1[n_subject_1, nts1 + nts2];            	// experimental kelp consumed observations; [59, 3]
  real A_obs_2[n_subject_2, nts3 + nts4 + nts5];        // experimental kelp consumed observations; [59, 3]
}

transformed data {
  int x_i[0]; 				 
  real x_r[3];              
  int n_total_1 = nts1 + nts2;				// n_total = 1 + 1 = 2
  int n_total_2 = nts3 + nts4 + nts5;			// n_total = 1 + 1 + 1 = 3	
}

parameters {
   real <lower=0, upper=.1> a; 
   real <lower=0, upper=10> v;
   real <lower=0, upper=10> q;
   real <lower=0, upper=1> p;
   //real <lower=0, upper=1> e;
   real <lower=0, upper=40> sigma;      
}

transformed parameters {
  real<lower=0> theta[4];
  vector<lower=0>[3] y1[nts1];					// two-dimensional container of size (nts1, 3) i.e. y1[1, 3] 
  vector<lower=0>[3] y2[nts2];					// two-dimensional container of size (nts2, 3)   
  vector<lower=0>[3] y3[nts3];					// two-dimensional container of size (nts3, 3)   
  vector<lower=0>[3] y4[nts4];					// two-dimensional container of size (nts4, 3)     
  vector<lower=0>[3] y5[nts5];					// two-dimensional container of size (nts5, 3)  
  vector<lower=0>[3] init_1;
  vector<lower=0>[3] init_2;
  vector<lower=0>[1] U;			
  
  // obs 
  real <lower=0> drift_loss_1[n_subject_1, n_total_1]; 		// Drift consumed; [59, 3]
  real <lower=0> kelp_loss_1[n_subject_1, n_total_1]; 		// Kelp consumed; [59, 3]
  real <lower=0> drift_loss_2[n_subject_2, n_total_2]; 		// Drift consumed; [59, 3]
  real <lower=0> kelp_loss_2[n_subject_2, n_total_2]; 		// Kelp consumed; [59, 3]
    
  // drift alpha, beta 
  real <lower=0> alphaS_1[n_subject_1, n_total_1];   		// reparameterized shape for gamma distribution; drift likelihood; [59, 3] 
  real <lower=0> betaS_1[n_subject_1, n_total_1];		// reparameterized scale for gamma distribution; drift likelihood; [59, 3]
  real <lower=0> alphaS_2[n_subject_2, n_total_2];   		// reparameterized shape for gamma distribution; drift likelihood; [59, 3] 
  real <lower=0> betaS_2[n_subject_2, n_total_2];  

  // kelp alpha, beta
  real <lower=0> alphaA_1[n_subject_1, n_total_1];   		// reparameterized shape for gamma distribution; kelp likelihood; [59, 3] 
  real <lower=0> betaA_1[n_subject_1, n_total_1];		// reparameterized scale for gamma distribution; kelp likelihood;[59, 3]
  real <lower=0> alphaA_2[n_subject_2, n_total_2];   		// reparameterized shape for gamma distribution; kelp likelihood; [59, 3] 
  real <lower=0> betaA_2[n_subject_2, n_total_2];
  
  theta[1] = a; 
  theta[2] = v;
  theta[3] = q;
  theta[4] = p;
  //theta[5] = e;


  // Temporal sequence 1 -----------------------------------
  for (i in 1:n_subject_1) {

    // period 1
    init_1[1] = y1_init_s_a[i, 1];
    init_1[2] = y1_init_s_a[i, 2];
    init_1[3] = 0;
    U[1] = y1_init_s_a[i, 3];
    y1 = ode_rk45(resourceLoss, init_1, t0_1, ts1, theta, U);
    drift_loss_1[i, 1:nts1] = y1[, 1];
    kelp_loss_1[i, 1:nts1] = y1[, 2];
	
    // period 2
    init_1[1] = y2_init_s_a[i, 1]; 		
    init_1[2] = y2_init_s_a[i, 2];
    init_1[3] = y1[nts1, 3];
    U[1] = y2_init_s_a[i, 3];
    y2 = ode_rk45(resourceLoss, init_1, ts1[nts1], ts2, theta, U);
    //print(i, " ", init_1, " prm[a,v,q,p]: ", theta, " y1: ", y1, " y2: ", y2, " dXdt: ", resourceLoss(t0_1, init_1, theta, U), " U: ", U);
    drift_loss_1[i, (nts1+1):(nts1+nts2)] = y2[, 1];
    kelp_loss_1[i, (nts1+1):(nts1+nts2)] = y2[, 2];
   
    for (j in 1:n_total_1) {
      alphaS_1[i, j] = pow (drift_loss_1[i, j], 2) / pow (sigma, 2);
      betaS_1[i, j] = (1 / (pow (sigma, 2) / drift_loss_1[i, j]));
      alphaA_1[i, j] = pow (kelp_loss_1[i, j], 2) / pow (sigma, 2);
      betaA_1[i, j] = (1 / (pow (sigma, 2) / kelp_loss_1[i, j]));            
    }
  }

  // Temporal sequence 2 -----------------------------------
  for (i in 1:n_subject_2) {      
   
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

    for (j in 1:n_total_2) {
      alphaS_2[i, j] = pow (drift_loss_2[i, j], 2) / pow (sigma, 2);
      betaS_2[i, j] = (1 / (pow (sigma, 2) / drift_loss_2[i, j]));
      alphaA_2[i, j] = pow (kelp_loss_2[i, j], 2) / pow (sigma, 2);
      betaA_2[i, j] = (1 / (pow (sigma, 2) / kelp_loss_2[i, j]));            
    }
  }
}

model { 
  a ~ exponential(0.1);                
  v ~ exponential(0.1); 
  q ~ lognormal(1, 1);
  p ~ beta(1,1);
  //e ~ uniform(0,1) 
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