//--- Time-varying dengue catalytic model ---//
// assumes constant endemic FOI prior to data
// assumes complete immunity after 2nd infection
// assumes equal transmissability of 4 serotypes
// partial pooling at province-district level

data {
  
  int nA; // N age groups
  int nT; // N time points
  int nD; // N admin2 within a province
  int cases[nD,nT,nA]; // reported case data
  int pop[nD,nT,nA]; // population data
  int ageLims[2,nA]; // lower & upper bounds of age groups
  int hist_length; // length of historical lambda
  
}

parameters {
  
  // real<lower=0,upper=0.15> lam_H_prov[hist_length]; // prov level historic FOI
  real<lower=0,upper=0.15> lam_H_prov; // prov level historic FOI
  vector<lower=0,upper=0.15>[hist_length] lam_H[nD]; // historic FOI
  real<lower=0,upper=0.15> lam_t_prov[nT]; // prov level time varying FOI
  // real<lower=0,upper=0.15> lam_t_prov; // prov level time varying FOI
  vector<lower=0,upper=0.15>[nT] lam_t[nD]; // time varying FOI
  real<lower=0,upper=1> rho_prov; // prov level reporting rate of 2nd infections
  real<lower=0,upper=1> rho[nD]; // reporting rate of 2nd infections
  real<lower=0,upper=1> gamma_prov; // prov level relative reporting rate of 1st infections
  real<lower=0,upper=1> gamma[nD]; // relative reporting rate of 1st infections
  real<lower=0> kappa_lam_H; // standard deviation for district level sampling of lam_H
  real<lower=0> kappa_lam_t; // standard deviation for district level sampling of lam_t
  real<lower=0> kappa_rho; // standard deviation for district level sampling of rho
  real<lower=0> kappa_gamma; // standard deviation for district level sampling of gamma

}

transformed parameters {
  
  matrix<lower=0,upper=1>[nT+1,100] susc[nD]; // proportion susceptible
  matrix<lower=0,upper=1>[nT+1,100] mono[nD]; // proportion monotypic
  matrix<lower=0,upper=1>[nT,100] inc1[nD]; // incidence of primary infections
  matrix<lower=0,upper=1>[nT,100] inc2[nD]; // incidence of secondary infections
  matrix<lower=0>[nT,nA] Ecases[nD]; // expected reported cases
  real sum_lamH[nD];
  
  //--- immune profiles at beginning of time series
  for(d in 1:nD) for (i in 1:100){
    
    if (i > hist_length) {
      sum_lamH[d] = sum(lam_H[d][1:hist_length]);
    } else {
      // lam_H_vec[i] = lam_H;
      sum_lamH[d] = sum(lam_H[d][1:i]);
    }
    
    susc[d][1,i] = exp(-4*sum_lamH[d]);
    mono[d][1,i] = 4*exp(-3*sum_lamH[d])*(1-exp(-sum_lamH[d]));
      
  }
  
  //--- subsequent time steps
  for(d in 1:nD) for (t in 2:(nT+1)){
    
    susc[d][t,1] = exp(-4*lam_t[d][t-1]);
    mono[d][t,1] = 4*exp(-3*lam_t[d][t-1])*(1-exp(-lam_t[d][t-1]));
    susc[d][t,2:100] = susc[d][t-1,1:99] - 4*lam_t[d][t-1]*susc[d][t-1,1:99];
    mono[d][t,2:100] = mono[d][t-1,1:99] + 4*lam_t[d][t-1]*susc[d][t-1,1:99] - 3*lam_t[d][t-1]*mono[d][t-1,1:99];
    inc1[d][t-1,1] = 4*lam_t[d][t-1]*1;
    inc2[d][t-1,1] = 3*lam_t[d][t-1]*0;
    inc1[d][t-1,2:100] = 4*lam_t[d][t-1]*susc[d][t-1,1:99];
    inc2[d][t-1,2:100] = 3*lam_t[d][t-1]*mono[d][t-1,1:99];
    
  }
  
  //--- reported cases
  // index of age + 1 as the age matrix starts from 0 and ends at 99
  for(d in 1:nD) for(t in 1:nT) for(a in 1:nA){
    Ecases[d][t,a] = rho[d]*(mean(inc2[d][t,(ageLims[1,a]+1):(ageLims[2,a]+1)]) + gamma[d]*mean(inc1[d][t,(ageLims[1,a]+1):(ageLims[2,a]+1)]))*pop[d][t,a];
  }

}

model {
  
  // think about variations in reporting rate for adults, using 0-4 years old reporting as a reference
  
  //--- priors
  // kappa for standard deviation
  // lam_H as long as historical length
  lam_H_prov ~ uniform(0,1); // normal(0,0.05);
  kappa_lam_H ~ normal(0,1);
  // for (d in 1:nD) for (t in 1:hist_length) lam_H[d][t] ~ normal(lam_H_prov[t],kappa_lam_H);
  for (d in 1:nD) lam_H[d] ~ normal(lam_H_prov,kappa_lam_H);
  lam_t_prov ~ uniform(0,1); // normal(0,0.05);
  kappa_lam_t ~ normal(0,1);
  for (d in 1:nD) for (t in 1:nT) lam_t[d][t] ~ normal(lam_t_prov[t],kappa_lam_t);
  // for (d in 1:nD) lam_t[d] ~ normal(lam_t_prov,0.01);
  rho_prov ~ beta(10,100); // normal(0,0.25)
  kappa_rho ~ normal(0,1);
  rho ~ normal(rho_prov,kappa_rho);
  gamma_prov ~ beta(30,40); // normal(0,0.25);
  kappa_gamma ~ normal(0,1);
  gamma ~ normal(gamma_prov,kappa_gamma);
  
  //--- likelihood 
  for(d in 1:nD) for(t in 1:nT) {
    cases[d][t,] ~ poisson(Ecases[d][t,]); // poisson likelihood
  }
  
}

// generated quantities{
//   matrix<lower=0>[nT,nA] Estim_cases[nD];
//   
//   for(d in 1:nD) for(t in 1:nT){
//     Estim_cases[d][t,] = to_row_vector(poisson_rng(Ecases[d][t,]));
//   }
//   
// }


