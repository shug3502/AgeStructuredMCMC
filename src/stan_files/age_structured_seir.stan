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
  int<lower=1> Tlockdown; //time point (day) at which social distancing measures introduced
  int<lower=1> nAges; //21
  int nCompartments; //16; //1 Susceptibles, 4 exposed, 5 detected, 3 undetected; plus hospitalised and removed and recovered
  matrix[nAges,nAges] beta[4]; //transmission rates: these depend on contacts between ages, and vary between household, school, work and other. Get these from PREM et al 2017 Plos Comp Biol
  real N[nAges]; //total population size of each age compartment
  real prop_hospitalised[nAges];
  real prop_die[nAges];
  int hospitalised[T,nAges];
  int deaths[T,nAges];
  real y0[nCompartments*nAges]; //initial condition age stratified
}

transformed data {
  real prop_service=0.3; //proportion of service industries, assumed as 0.3
  matrix[4,nAges] Q; //max compliance benefits from social distancing
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

  Q[1] = [1.25, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25];
  Q[2] = [0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
  Q[3] = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
  Q[4] = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
}

parameters {
  real<lower=0> alpha; //rate of disease progression from exposed to infected
  real<lower=0> tau; //relative transmission from asymptomatics vs symptomatics
  real<lower=0> beta_scaling; //check if you need this; may involve reading Prem et al 2017 and or asking louise
  real<lower=0> gamma; //rate of removal/recovery/hospitalisation
  real<lower=0,upper=1> H; //fraction of first detected case in household quarantined
  real<lower=0,upper=1> phi; //compliance level
  real<lower=0,upper=1> det_rate[nAges]; //rate of detection across the population; should constrain such that this is less than 80%
  real<lower=0> delta; //rate of progression from hospitalisation to death
  real<lower=0> obs_noise; //observation noise
}

transformed parameters {
  matrix[nAges,nAges] beta_with_social_distancing[4];
  matrix[4,nAges] distancing_compliance;
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


  for (j in 1:nAges){
    for (i in 1:4){
      //linearly interpolate between 1 (no change) and the maximum effect of social distancing according to the compliance level, phi
      distancing_compliance[i,j] = 1 + phi*(Q[i,j] - 1);
    }
  }
  for (b in 1:nAges){
    for (a in 1:nAges){
      beta_with_social_distancing[1][a,b] = beta[1][a,b]*distancing_compliance[1,a];
      beta_with_social_distancing[2][a,b] = (1-prop_service)*beta[2][a,b]*distancing_compliance[2,a] + prop_service*beta[2][a,b].*distancing_compliance[2,a]*distancing_compliance[4,b];
      beta_with_social_distancing[3][a,b] = beta[3][a,b]*distancing_compliance[3,a];
      beta_with_social_distancing[4][a,b] = beta[4][a,b]*distancing_compliance[4,a]*distancing_compliance[4,b];
      //TODO: check what the bracket-phi notation means to make sure not multiplying by phi twice or scaling for compliance wrongly
    }
  }

}
model {
  real z[T,nCompartments*nAges];
  //priors
  alpha ~ normal(1,0.01);
  tau ~ normal(0.5,0.01);
  beta_scaling ~ normal(0.1,0.01);
  gamma ~ normal(0.33,0.01);
  H ~ beta(3,1); // beta expectation: a/a+b
  phi ~ beta(2.5,1);
  delta ~ normal(1,0.01);
  obs_noise ~ normal(1,0.01);
  for (a in 1:14){
    det_rate[a] ~ normal(0.2,0.01);
  }
  for (a in 15:nAges){
    det_rate ~ normal(0.5,0.01);
  }
/////////////////////////////
//to incorporate lockdown and social distancing, split up the time series and integrate the ode in two portions
  z = integrate_ode_rk45(age_structed_seir, y0, 0, ts, theta, x_r, x_i, 10^-4, 1, 10^3);

//observations only available for hospitalisations and deaths
  for (t in 1:T){
    for (a in 1:nAges){
      hospitalised[t,a] ~ neg_binomial_2(z[t,13*nAges+a], obs_noise);
      deaths[t,a] ~ neg_binomial_2(z[t,15*nAges+a], obs_noise);
    }
  }
  //need to use lag tables for symptoms to hospital, and hospital to death to inform gamma and delta
}

