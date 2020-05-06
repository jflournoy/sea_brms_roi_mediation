// generated with brms 2.12.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  int<lower=1> N_Y;  // number of observations
  vector[N_Y] Y_Y;  // response variable
  int<lower=1> K_Y;  // number of population-level effects
  matrix[N_Y, K_Y] X_Y;  // population-level design matrix
  int<lower=1> N_brain;  // number of observations
  vector[N_brain] Y_brain;  // response variable
  int<lower=1> K_brain;  // number of population-level effects
  matrix[N_brain, K_brain] X_brain;  // population-level design matrix
  int<lower=1> nresp;  // number of responses
  int nrescor;  // number of residual correlations
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1_Y[N_Y];  // grouping indicator per observation
  int<lower=1> J_1_brain[N_brain];  // grouping indicator per observation
  // group-level predictor values
  vector[N_Y] Z_1_Y_1;
  vector[N_Y] Z_1_Y_2;
  vector[N_Y] Z_1_Y_3;
  vector[N_brain] Z_1_brain_4;
  vector[N_brain] Z_1_brain_5;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_Y = K_Y - 1;
  matrix[N_Y, Kc_Y] Xc_Y;  // centered version of X_Y without an intercept
  vector[Kc_Y] means_X_Y;  // column means of X_Y before centering
  int Kc_brain = K_brain - 1;
  matrix[N_brain, Kc_brain] Xc_brain;  // centered version of X_brain without an intercept
  vector[Kc_brain] means_X_brain;  // column means of X_brain before centering
  vector[nresp] Y[N];  // response array
  for (i in 2:K_Y) {
    means_X_Y[i - 1] = mean(X_Y[, i]);
    Xc_Y[, i - 1] = X_Y[, i] - means_X_Y[i - 1];
  }
  for (i in 2:K_brain) {
    means_X_brain[i - 1] = mean(X_brain[, i]);
    Xc_brain[, i - 1] = X_brain[, i] - means_X_brain[i - 1];
  }
  for (n in 1:N) {
    Y[n] = [Y_Y[n], Y_brain[n]]';
  }
}
parameters {
  vector[Kc_Y] b_Y;  // population-level effects
  real Intercept_Y;  // temporary intercept for centered predictors
  real<lower=0> sigma_Y;  // residual SD
  vector[Kc_brain] b_brain;  // population-level effects
  real Intercept_brain;  // temporary intercept for centered predictors
  real<lower=0> sigma_brain;  // residual SD
  cholesky_factor_corr[nresp] Lrescor;  // parameters for multivariate linear models
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_Y_1;
  vector[N_1] r_1_Y_2;
  vector[N_1] r_1_Y_3;
  vector[N_1] r_1_brain_4;
  vector[N_1] r_1_brain_5;
  // compute actual group-level effects
  r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)';
  r_1_Y_1 = r_1[, 1];
  r_1_Y_2 = r_1[, 2];
  r_1_Y_3 = r_1[, 3];
  r_1_brain_4 = r_1[, 4];
  r_1_brain_5 = r_1[, 5];
}
model {
  // initialize linear predictor term
  vector[N_Y] mu_Y = Intercept_Y + Xc_Y * b_Y;
  // initialize linear predictor term
  vector[N_brain] mu_brain = Intercept_brain + Xc_brain * b_brain;
  // multivariate predictor array
  vector[nresp] Mu[N];
  vector[nresp] sigma = [sigma_Y, sigma_brain]';
  // cholesky factor of residual covariance matrix
  matrix[nresp, nresp] LSigma = diag_pre_multiply(sigma, Lrescor);
  for (n in 1:N_Y) {
    // add more terms to the linear predictor
    mu_Y[n] += r_1_Y_1[J_1_Y[n]] * Z_1_Y_1[n] + r_1_Y_2[J_1_Y[n]] * Z_1_Y_2[n] + r_1_Y_3[J_1_Y[n]] * Z_1_Y_3[n];
  }
  for (n in 1:N_brain) {
    // add more terms to the linear predictor
    mu_brain[n] += r_1_brain_4[J_1_brain[n]] * Z_1_brain_4[n] + r_1_brain_5[J_1_brain[n]] * Z_1_brain_5[n];
  }
  // combine univariate parameters
  for (n in 1:N) {
    Mu[n] = [mu_Y[n], mu_brain[n]]';
  }
  // priors including all constants
  target += normal_lpdf(b_Y | 0, 20);
  target += normal_lpdf(Intercept_Y | 0, 20);
  target += weibull_lpdf(sigma_Y | 2, 1);
  target += normal_lpdf(b_brain | 0, 20);
  target += normal_lpdf(Intercept_brain | 0, 20);
  target += weibull_lpdf(sigma_brain | 2, 1);
  target += lkj_corr_cholesky_lpdf(Lrescor | 1);
  target += cauchy_lpdf(sd_1 | 0,2.5)
    - 5 * cauchy_lccdf(0 | 0,2.5);
  target += normal_lpdf(to_vector(z_1) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
  // likelihood including all constants
  if (!prior_only) {
    target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Y_Intercept = Intercept_Y - dot_product(means_X_Y, b_Y);
  // actual population-level intercept
  real b_brain_Intercept = Intercept_brain - dot_product(means_X_brain, b_brain);
  // residual correlations
  corr_matrix[nresp] Rescor = multiply_lower_tri_self_transpose(Lrescor);
  vector<lower=-1,upper=1>[nrescor] rescor;
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:nresp) {
    for (j in 1:(k - 1)) {
      rescor[choose(k - 1, 2) + j] = Rescor[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
