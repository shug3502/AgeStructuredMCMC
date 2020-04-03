functions {

real[] age_structed_seir(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    int nAges=x_i[1];
    int nCompartments=x_i[2];
    matrix[4,nAges] lambda;
    real N[nAges];
    vector[nCompartments*nAges] dydt;
    for (a in 1:nAges){
      N[a] = x_r[a];
    }
    lambda = rep_matrix(1,4,nAges);
    //age structured model: differential equations
    for (a in 1:nAges){
    //susceptibles
    dydt[a] = -(lambda[1][a] + lambda[2][a] + lambda[3][a] + lambda[4][a])*y[a]/N[a];
    //exposed
    dydt[nAges+a] = 0;
    dydt[2*nAges+a] = 0;
    dydt[3*nAges+a] = 0; //TODO: check that the typo is in the equations written in supplementary
    //detected
    dydt[4*nAges+a] = 0;
    dydt[5*nAges+a] = 0;
    dydt[6*nAges+a] = 0;
    //undetected
    dydt[7*nAges+a] = 0;
    dydt[8*nAges+a] = 0;
    //quarantined
    dydt[9*nAges+a] = 0;
    dydt[10*nAges+a] = 0;
    dydt[11*nAges+a] = 0;
    dydt[12*nAges+a] = 0;
    //hospitalised
    dydt[13*nAges+a] = 0;
    //recovered
    dydt[14*nAges+a] = 0;
    //deaths
    dydt[15*nAges+a] = 0;
    }
    return to_array_1d(dydt);
  }
}

data {
  int<lower=1> T; //number of time points (days) data available for
  real ts[T];
  int<lower=1> nAges; //21
  matrix[nAges,nAges] beta[4]; //transmission rates: these depend on contacts between ages, and vary between household, school, work and other. Get these from PREM et al 2017 Plos Comp Biol
  real N[nAges]; //total population size of each age compartment
  real prop_hospitalised[nAges];
  real prop_die[nAges];
  int hospitalised[T,nAges];
  int deaths[T,nAges];
  real theta[7+nAges];
  real obs_noise;
}

transformed data {
  int nCompartments=16; //1 Susceptibles, 4 exposed, 5 detected, 3 undetected; plus hospitalised and removed and recovered
  real prop_service=0.3; //proportion of service industries, assumed as 0.3
  real sigma=1; //assume age specific differences in transmission
  real y0[nCompartments*nAges]; //initial condition age stratified
  real x_r[nAges]; //can use these to input into DE solver
  int x_i[2];
  x_i[1] = nAges;
  x_i[2] = nCompartments;

  for (a in 1:nAges){
    x_r[a] = N[a];
    y0[a] = 1 - 1/N[a];
    y0[nAges+a] = 1/N[a];
    for (i in 3:nCompartments){
      y0[(i-1)*nAges+a] = 0;
    }
  }
}

parameters {
}

transformed parameters {
}

model {
}

generated quantities {
  real z[T,nCompartments*nAges];

  z = integrate_ode_rk45(age_structed_seir, y0, 0, ts, theta, x_r, x_i);
}
