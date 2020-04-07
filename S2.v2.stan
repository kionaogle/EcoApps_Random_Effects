data {
  // Dimensions
  int<lower=1> m;     // # of units (i.e., # groups = 10)
  int<lower=1> ni[m]; // # obs in each group (ni=10, all i)
  int<lower=1> N;     // total # of obs (sum(ni) = 100)
  int<lower=1> p;     // # of fixed effects: intercept and slope (p=2)
  int<lower=1> q;     // # of random effects per unit: intercept (q=1)

  // Data
  vector[N] y;   // observed response values
  matrix[N,p] X; // fixed effects regression matrix (1's augmented
		 // with covariate (p=2))
  matrix[N,q] Z; // random effects regression matrix (q=1 column of
		 // 1's)

  // Fixed effects prior hyper-parameters
  vector[p] mu_b;
  cov_matrix[p] Sig_b;

  // Random effects prior hyper-parameters
  vector<lower=0>[q] tau_eps_a;
  vector<lower=0>[q] tau_eps_b;
  vector[q] mu_eps;

  // Error precision prior hyper-parameters
  real<lower=0> tau_a;
  real<lower=0> tau_b;
}

parameters {
  vector[p] b;                // fixed effects b0 and b1 (p=2)
  matrix[q,m] eps;            // random effects epsj j=1 to m=10 (q=1)
  vector<lower=0>[q] tau_eps; // random effects precision
  real<lower=0> tau;          // obs error precision

}

transformed parameters {
  real<lower=0> sig = pow(sqrt(tau),-1);             // obs error sd
  vector<lower=0>[q] sig_eps = 1.0 ./ sqrt(tau_eps); // ran. eff. sd
}

model{
  int pos = 1;
  matrix[q,q] Sig_eps = diag_matrix(square(sig_eps)); // ran. eff. variance
  vector[N] mu_y = X*b; // fixed effects contribution to mean
  
  b ~ multi_normal(mu_b, Sig_b);         // fixed effect prior
  tau_eps ~ gamma(tau_eps_a, tau_eps_b); // random effects precision
					 // prior
  tau ~ gamma(tau_a, tau_b);             // obs error precision prior
  for(i in 1:m) {
    eps[,i] ~ multi_normal(mu_eps, Sig_eps); // random effects
    segment(y,pos,ni[i]) ~ normal(segment(mu_y,pos,ni[i]) +
				  Z[pos:(pos+ni[i]-1),]*eps[,i],
				  sig);
    pos += ni[i];
  }
}
