// ZIP-Hidden Markov Model 
// Asumiendo todos los estados con cero inflacion

data {
  int<lower=0> N;    // length of chain
  int<lower=0> y[N]; // emissions
  int<lower=1> m;    // num states
}

parameters {
  simplex[m] start_pos;         // initial pos probs
  real<lower=0, upper=1> theta; // zero-inflation parameter
  positive_ordered[m] lambda;   // emission poisson params
  simplex[m] Gamma[m];          // transition prob matrix
}

model {
  vector[m] log_Gamma_tr[m];
  vector[m] lp;
  vector[m] lp_p1;

  // transposing tpm and taking the log of each entry
  for (i in 1:m) {
    for (j in 1:m) { 
      log_Gamma_tr[j, i] = log(Gamma[i, j]);
    }
  }

  // initial position log-lik
  lp = log(start_pos);

  for (n in 1:N) {
    for (j in 1:m) {
      // log-lik for state
      lp_p1[j] = log_sum_exp(log_Gamma_tr[j] + lp);

      // log-lik for emission
      // assuming all states has zero-inflation
        if (y[n] == 0) {
          lp_p1[j] += log_mix(theta, 0, poisson_lpmf(0 | lambda[j]));
        } else {
          lp_p1[j] += log1m(theta) + poisson_lpmf(y[n] | lambda[j]);
        }
    }
    lp = lp_p1; // log-lik for next position
  }
  target += log_sum_exp(lp);
}
