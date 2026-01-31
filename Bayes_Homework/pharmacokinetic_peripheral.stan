functions {
  vector pck(real t, vector u, real cl, real vc, real ka, real sigma, real q, real vp){
    vector[3] du_dt;
    
    du_dt[1] = -ka * u[1];
    du_dt[2] = ka * u[1] - (cl / vc + q / vc) * u[2] + q / vp * u[3];
    du_dt[3] = q / vc * u[2] - q / vp * u[3];
    
    return du_dt;
  }
}

data {
  real t0;
  vector[3] u0;
  int<lower = 1> n;
  array[n] real y;
  array[n] real t;
  
}

parameters {
  real<lower = 0> cl;
  real<lower = 0> vc;
  real<lower = 0> ka;
  real<lower = 0> sigma;
  real<lower = 0> q;
  real<lower = 0> vp;
  
}

transformed parameters {
  array[n] vector[3] u = ode_rk45(pck, u0, t0, t, cl, vc, ka, sigma, q, vp);
  
}

model {
  cl ~ lognormal(log(10), 0.25);
  vc ~ lognormal(log(35), 0.25);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ normal(0, 1);
  q ~ lognormal(log(15), sqrt(0.5));
  vp ~ lognormal(log(105), sqrt(0.5));
  
  for (i in 1:n)
    y[i] ~ lognormal(log(u[i, 2] / vc), sigma);
  
}

generated quantities {
  array[n] vector[20] y_pred_arr;
  
  for (i in 1:n){
    y_pred_arr[i] = rep_vector(0, 20);
    for(j in 1:20)
      y_pred_arr[i, j] = lognormal_rng(log(u[i, 2] / vc), sigma);
  }
  
  array[n] real y_pred;
  for (i in 1:n)
    y_pred[i] = mean(y_pred_arr[i]);
    
  vector[n] loglik;
  for(i in 1:n)
    loglik[i] = lognormal_lpdf(y[i] | log(u[i, 2] / vc), sigma);
  
}
