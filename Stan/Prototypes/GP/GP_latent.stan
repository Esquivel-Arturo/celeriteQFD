functions {
  real generalized_inverse_gaussian_lpdf(real x, int p,
                                        real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
 }
}

data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
  //real p;// for gig prior
  //real a;
  //real b;
}
transformed data {
  real delta = 1e-9;
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  //vector[N] eta;
  vector[N] f;
}

model {
  //vector[N] f;

    //matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(x, alpha, rho);

    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;

    //L_K = cholesky_decompose(K);
    //f = L_K * eta;
  f ~ multi_normal(0 * y, K);
  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  // alternative
  // target += generalized_inverse_gaussian_lpdf(rho, p, a, b)
  sigma ~ std_normal();
  //eta ~ std_normal();

  y ~ normal(f, sigma);
}

generated quantities{
  vector[N] y_detrend;
  y_detrend = y - f;
}
