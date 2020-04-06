data {
  int nAges;
}
parameters {
  real<lower=0> alpha; //rate of disease progression from exposed to infected
  real<lower=0> tau; //relative transmission from asymptomatics vs symptomatics
  real<lower=0> sigma;
  real<lower=0> gamma; //rate of removal/recovery/hospitalisation
  real<lower=0,upper=1> H; //fraction of first detected case in household quarantined
  real<lower=0,upper=1> phi; //compliance level
  real<lower=0,upper=1> det_rate[nAges]; //rate of detection across the population; should constrain such that this is less than 80%
  real<lower=0> delta; //rate of progression from hospitalisation to death
  real<lower=0> obs_noise; //observation noise
}
model {
  //priors
  alpha ~ normal(1,1);
  tau ~ normal(0,1);
  sigma ~ normal(0,1);
  gamma ~ normal(0.33,0.1);
  H ~ beta(3,1); // beta expectation: a/a+b
  phi ~ beta(2.5,1);
  delta ~ normal(1,1);
  obs_noise ~ normal(1,1);
  to_vector(det_rate) ~ normal(0.2,0.2);
}
