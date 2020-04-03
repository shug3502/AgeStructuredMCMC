functions {

real[] age_structed_seir(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    int nAges=x_i[1];
    int nCompartments=x_i[2];
    real N[nAges];
    real sigma;
    real prop_hospitalised[nAges];
    real prop_die[nAges];
    matrix[4,nAges] lambda;
    matrix[nAges,nAges] beta[4];
    vector[nCompartments*nAges] dydt;
//pass input from other environment via x_r and x_i
for (i in 1:4){
  for (b in 1:nAges){
    for (a in 1:nAges){
      beta[i][a,b] = x_r[(i-1)*nAges*nAges + (b-1)*nAges + a];
    }
  }
}
for (a in 1:nAges){
 N[a] = x_r[4*nAges*nAges+a];
 prop_hospitalised[a] = x_r[4*nAges*nAges+nAges+a];
 prop_die[a] = x_r[4*nAges*nAges+2*nAges+a];
}
sigma = x_r[4*nAges*nAges + 3*nAges + 1];

//compute lambda which depends on the compartments
//TODO:check if beta scaling needs to go in here
    lambda = rep_matrix(0,4,nAges);
    for (a in 1:nAges){
      for (b in 1:nAges){
        lambda[1,a] += sigma*(y[4*nAges+b] + y[5*nAges+b] + y[6*nAges+b] + theta[2]*y[7*nAges+b] + theta[2]*y[8*nAges+b])*(beta[2][b,a]+beta[3][b,a]+beta[4][b,a])*theta[3];
        lambda[2,a] += sigma*(y[4*nAges+b]*beta[1][b,a])*theta[3];
        lambda[3,a] += sigma*(theta[2]*y[7*nAges+b]*beta[1][b,a])*theta[3];
        lambda[4,a] += sigma*(y[10*nAges+b]*beta[1][b,a])*theta[3];
      }
    }

    //age structured model: differential equations
    for (a in 1:nAges){
    //susceptibles
    dydt[a] = -(lambda[1][a] + lambda[2][a] + lambda[3][a] + lambda[4][a])*y[a]/N[a];
    //exposed
    dydt[nAges+a] = lambda[1][a]*y[a]/N[a] - theta[1]*y[nAges+a];
    dydt[2*nAges+a] = lambda[2][a]*y[a]/N[a] - theta[1]*y[2*nAges+a];
    dydt[3*nAges+a] = lambda[3][a]*y[a]/N[a] - theta[1]*y[3*nAges+a]; //TODO: check that the typo is in the equations written in supplementary
    //detected
    dydt[4*nAges+a] = theta[6+a]*(1-theta[5])*theta[1]*y[nAges+a] - theta[4]*y[4*nAges+a];
    dydt[5*nAges+a] = theta[6+a]*theta[1]*y[2*nAges+a] - theta[4]*y[5*nAges+a];
    dydt[6*nAges+a] = theta[6+a]*(1-theta[5])*theta[1]*y[3*nAges+a] - theta[4]*y[6*nAges+a];
    //undetected
    dydt[7*nAges+a] = (1-theta[6+a])*theta[1]*y[nAges+a] - theta[4]*y[7*nAges+a];
    dydt[8*nAges+a] = (1-theta[6+a])*theta[1]*(y[2*nAges+a]+y[3*nAges+a]) - theta[4]*y[8*nAges+a];
    //quarantined
    dydt[9*nAges+a] = lambda[4][a]*y[a]/N[a] - theta[1]*y[9*nAges+a];
    dydt[10*nAges+a] = theta[6+a]*theta[5]*theta[1]*y[nAges+a] - theta[4]*y[10*nAges+a];
    dydt[11*nAges+a] = theta[6+a]*theta[5]*theta[1]*y[3*nAges+a] + theta[6+a]*theta[1]*y[9*nAges+a] - theta[4]*y[11*nAges+a];
    dydt[12*nAges+a] = (1-theta[6+a])*theta[1]*y[9*nAges+a] - theta[4]*y[12*nAges+a];
    //hospitalised
    dydt[13*nAges+a] = theta[4]*prop_hospitalised[a]*(y[4*nAges+a] + y[5*nAges+a] + y[6*nAges+a] + y[7*nAges+a] + y[8*nAges+a] + y[10*nAges+a] + y[11*nAges+a] + y[12*nAges+a]) - theta[7+nAges]*y[13*nAges+a];
    //recovered
    dydt[14*nAges+a] = theta[4]*(1-prop_hospitalised[a])*(y[4*nAges+a] + y[5*nAges+a] + y[6*nAges+a] + y[7*nAges+a] + y[8*nAges+a] + y[10*nAges+a] + y[11*nAges+a] + y[12*nAges+a]) + (1-prop_die[a])*theta[7+nAges]*y[13*nAges+a];
    //deaths
    dydt[15*nAges+a] = theta[7+nAges]*prop_die[a]*y[13*nAges+a];
    }

//notes: use symptoms to hospital to inform the rate gamma, assuming that progression is either to recovery or hospitalisation

    return to_array_1d(dydt);
  }
}

data {
  int<lower=1> T; //number of time points (days) data available for
  real ts[T];
  int<lower=1> nAges; //21
  int nCompartments; //16; //1 Susceptibles, 4 exposed, 5 detected, 3 undetected; plus hospitalised and removed and recovered
  matrix[nAges,nAges] beta[4]; //transmission rates: these depend on contacts between ages, and vary between household, school, work and other. Get these from PREM et al 2017 Plos Comp Biol
  real N[nAges]; //total population size of each age compartment
  real prop_hospitalised[nAges];
  real prop_die[nAges];
  real theta[7+nAges];
  real obs_noise;
  real y0[nCompartments*nAges]; //initial condition age stratified
}

transformed data {
  real prop_service=0.3; //proportion of service industries, assumed as 0.3
  real sigma=1; //assume age specific differences in transmission
  real x_r[4*nAges*nAges + 3*nAges + 1]; //can use these to input into DE solver
  int x_i[2];
  x_i[1] = nAges;
  x_i[2] = nCompartments;

for (i in 1:4){
  for (b in 1:nAges){
    for (a in 1:nAges){
      x_r[(i-1)*nAges*nAges + (b-1)*nAges + a] = beta[i][a,b];
    }
  }
}
for (a in 1:nAges){
  x_r[4*nAges*nAges+a] = N[a];
  x_r[4*nAges*nAges+nAges+a] = prop_hospitalised[a];
  x_r[4*nAges*nAges+2*nAges+a] = prop_die[a];
}
x_r[4*nAges*nAges+3*nAges+1] = sigma;

/*
  //initial condition: assume 1 case in each age category
  //may want to input this as data later
  for (a in 1:nAges){
    y0[a] = 1 - 1/N[a];
    y0[nAges+a] = 1/N[a];
    for (i in 3:nCompartments){
      y0[(i-1)*nAges+a] = 0;
    }
  }
*/

}

parameters {
  /*
  real<lower=0> alpha; //rate of disease progression from exposed to infected
  real<lower=0> tau; //relative transmission from asymptomatics vs symptomatics
  real<lower=0> beta_scaling; //check if you need this; may involve reading Prem et al 2017 and or asking louise
  real<lower=0> gamma; //rate of removal/recovery/hospitalisation
  real<lower=0,upper=1> H; //fraction of first detected case in household quarantined
  real<lower=0,upper=1> det_rate[nAges]; //rate of detection across the population; should constrain such that this is less than 80%
  real<lower=0> delta; //rate of progression from hospitalisation to death
  real<lower=0> obs_noise; //observation noise
  */
}

transformed parameters {
  /*
  real theta[7+nAges];
  theta[1] = alpha;
  theta[2] = tau;
  theta[3] = beta_scaling;
  theta[4] = gamma;
  theta[5] = H;
  theta[6] = phi;
  for (a in 1:nAges){
    theta[6+a] = det_rate[a];
  }
  theta[7+nAges] = delta;
*/
}
model {
  //need to use lag tables for symptoms to hospital, and hospital to death to inform gamma and delta
}

generated quantities {
  real z[T,nCompartments*nAges];
  int hospitalised[T,nAges];
  int deaths[T,nAges];
  /*
  //priors
  alpha ~ normal(0,1);
  tau ~ normal(0,1);
  beta_scaling ~ normal(1,0.01);
  gamma ~ normal(0,1);
  H ~ beta(1,1);
  phi ~ beta(1,1);
  delta ~ normal(0,1);
  to_vector(det_rate) ~ normal(0.5,0.1);
  obs_noise ~ normal(1,1);
  */
/////////////////////////////
//to incorporate lockdown and social distancing, split up the time series and integrate the ode in two portions
  z = integrate_ode_rk45(age_structed_seir, y0, 0, ts, theta, x_r, x_i);

//observations only available for hospitalisations and deaths
  for (t in 1:T){
    for (a in 1:nAges){
      hospitalised[t,a] = neg_binomial_2_rng(z[t,13*nAges+a], obs_noise);
      deaths[t,a] = neg_binomial_2_rng(z[t,15*nAges+a], obs_noise); //may want to convert to daily observed deaths in sim and inference
    }
  }
}
