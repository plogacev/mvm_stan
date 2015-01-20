data{
  int<lower=1> n_obs;
  int<lower=1> n_subj;
  int<lower=1> n_item;
  int<lower=0,upper=1> data_response[n_obs];
  real<lower=0> data_reading_time[n_obs];
  real<lower=0> data_response_time[n_obs];
  real<lower=0,upper=1> ind_np1[n_obs];
  real<lower=0,upper=1> ind_np2[n_obs];
}

parameters{
  real<lower=0,upper=1> par_p_err;
  real<lower=0,upper=1> par_p_u;
  real<lower=0,upper=1> par_p_np2;
}

model {
  matrix[7,n_obs] p_trial_type;
  
  // type 1A: attach N1 + correct response
  p_trial_type[1] <- (1-par_p_u)*(1-p_np2)*(1-par_p_err);
  
  for(i in 1:n_obs)
    data_response[i] ~ bernoulli(par_p_err);
  par_p_err ~ uniform(0,1);
  par_p_u ~ uniform(0,1);
  par_p_np2 ~ uniform(0,1);
}