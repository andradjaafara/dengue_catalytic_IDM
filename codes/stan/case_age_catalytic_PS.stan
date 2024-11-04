//--- Time-varying dengue catalytic model ---//
// assumes constant endemic FOI prior to data
// assumes complete immunity after 2nd infection
// assumes equal transmissability of 4 serotypes

data {
  
  int nA; // N age groups
  int nT; // N time points
  int cases[nT,nA]; // reported case data
  matrix[nT,nA] pop; // population data
  int ageLims[2,nA]; // lower & upper bounds of age groups
  int hist_length; // length of historical lambda
  // row_vector[100] age;
  
}


parameters {
  
  real<lower=0,upper=0.15> lam_H; // historic average FOI
  real<lower=0,upper=0.15> lam_t[nT]; // time varying FOI
  real<lower=0,upper=1> rho; // reporting rate of 2nd infections
  real<lower=0,upper=1> gamma; // relative reporting rate of 1st infections

}

transformed parameters {
  
  matrix<lower=0,upper=1>[nT+1,100] susc; // proportion susceptible
  matrix<lower=0,upper=1>[nT+1,100] mono; // proportion monotypic
  // matrix<lower=0,upper=1>[nT,100] multi; // proportion multitypic
  matrix<lower=0,upper=1>[nT,100] inc1; // incidence of primary infections
  matrix<lower=0,upper=1>[nT,100] inc2; // incidence of secondary infections
  vector<lower=0>[nA] Ecases[nT]; // expected reported cases
  vector<lower=0>[nA] Epropcases[nT]; // expected proportion of reported cases by age group
  real<lower=0,upper=0.15> lam_H_vec[hist_length];
  real sum_lamH;
  
  
  //--- immune profiles at beginning of time series
  for (i in 1:100){
    
    if (i > hist_length) {
      sum_lamH = sum(lam_H_vec[1:hist_length]);
    } else {
      lam_H_vec[i] = lam_H;
      sum_lamH = sum(lam_H_vec[1:i]);
    }
    
    susc[1,i] = exp(-4*sum_lamH);
    mono[1,i] = 4*exp(-3*sum_lamH)*(1-exp(-sum_lamH));
      
  }
  
  //--- subsequent time steps
  for (t in 2:(nT+1)){
    
    susc[t,1] = exp(-4*lam_t[t-1]);
    mono[t,1] = 4*exp(-3*lam_t[t-1])*(1-exp(-lam_t[t-1]));
    susc[t,2:100] = susc[t-1,1:99] - 4*lam_t[t-1]*susc[t-1,1:99];
    mono[t,2:100] = mono[t-1,1:99] + 4*lam_t[t-1]*susc[t-1,1:99] - 3*lam_t[t-1]*mono[t-1,1:99];
    inc1[t-1,1] = 4*lam_t[t-1]*1;
    inc2[t-1,1] = 3*lam_t[t-1]*0;
    inc1[t-1,2:100] = 4*lam_t[t-1]*susc[t-1,1:99];
    inc2[t-1,2:100] = 3*lam_t[t-1]*mono[t-1,1:99];
    
  }
  
  //--- reported cases
  // index of age + 1 as the age matrix starts from 0 and ends at 99
  for(t in 1:nT) for(a in 1:nA){
    Ecases[t,a] = rho*(mean(inc2[t,(ageLims[1,a]+1):(ageLims[2,a]+1)]) + gamma*mean(inc1[t,(ageLims[1,a]+1):(ageLims[2,a]+1)]))*pop[t,a];
  }
  
  for(t in 1:nT){
    Epropcases[t,] = Ecases[t,] ./ sum(Ecases[t,]);
  }

}

model {
  
  //--- priors
  lam_H ~ normal(0,0.05);
  lam_t ~ normal(0,0.05);
  rho ~ normal(0,0.25);
  gamma ~ normal(0,0.25);
  
  //--- likelihood 
  for(t in 1:nT) cases[t,] ~ poisson(Ecases[t,]); // poisson likelihood
  // for(t in 1:nT) cases[t,] ~ multinomial(Epropcases[t,]); // multinomial likelihood
  // for(t in 1:nT) cases[t,] ~ neg_binomial_2(Ecases[t,], phi); // negative binomial likelihood with estimated phi
  
}



