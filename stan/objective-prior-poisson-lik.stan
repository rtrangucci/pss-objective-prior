functions {
  real[] dz_dt(real t, real[] z, real[] theta,
               real[] x_r, int[] x_i) {
    real exp_m_u_theta = z[1];
  
    return { -sqrt(x_r[1] / exp_m_u_theta - 2 * (1 - log(exp_m_u_theta)))*exp_m_u_theta };
  }
}
data {
  real<lower=0> c;
  real u_0;
  int N;
  int y[N];
  real<lower=0> acc;
}
transformed data{
  real x_r[1] = {c};
  real z_init[1] = {exp(-u_0)};
  real dummy_theta[0];
}
parameters {
  real<lower=0> theta;
}
model {
  real ts[1] = {theta};
   real z[1, 1]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, rep_array(0.0,0),
                         x_r, rep_array(0, 0),
                         acc, acc, 1e5);
 target += log(z[1,1]);
 if (N > 0)
    y ~ poisson(theta);
}

