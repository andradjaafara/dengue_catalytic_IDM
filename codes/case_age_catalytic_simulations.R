library(tidyverse)
library(janitor)
library(rstan)
library(cowplot)
library(bayesplot)
source("codes/functions.R")
source("codes/pop_structure_prep.R")

options(mc.cores = parallel::detectCores())

#### simulate cases based on historical lambda and current lambda scenarios

#### current lambda scenarios
#### 1 total FoI < 0.1 with small variations for 10 years
#### 2 total FoI > 0.1 with small variations for 10 years
#### 3 total FoI can be anywhere between 0.04 to 0.24

#### historical lambda scenarios: 80 years
#### 1 constant low with total FoI < 0.1
#### 2 constant high with total FoI > 0.1
#### 3 low to high changing every 40 years
#### 4 low to high changing every 20 years
#### 5 low to high changing every 10 years

#### combinations of scenarios: 11 simulation scenarios
#### 1 - current with 1 - historical - scenario 1
#### 2 - current with 1,2,3,4,5 - historical - scenario 2-6
#### 3 - current with 1,2,3,4,5 - historical - scenario 7-11

#### set scenario parameters
#### current lambda
set.seed(12345)
lam_curr1 <- runif(10,0.005,0.025)
set.seed(29819)
lam_curr2 <- runif(10,0.025,0.06)
set.seed(98089)
lam_curr3 <- runif(10,0.01,0.06)

#### historical lambda
lam_H1 <- rep(0.02,80)
lam_H2 <- rep(0.03,80)
lam_H3 <- c(rep(0.02,40),rep(0.03,40))
lam_H4 <- c(rep(0.015,20),rep(0.02,20),rep(0.025,20),rep(0.03,20))
lam_H5 <- c(rep(0.010,10),rep(0.015,10),rep(0.020,10),rep(0.025,10),
            rep(0.0275,10),rep(0.03,10),rep(0.035,10),rep(0.04,10))

lam_curr_scenario_list <- list(lam_curr1,lam_curr2,lam_curr2,lam_curr2,lam_curr2,lam_curr2,
                               lam_curr3,lam_curr3,lam_curr3,lam_curr3,lam_curr3)
lam_H_scenario_list <- list(lam_H1,lam_H1,lam_H2,lam_H3,lam_H4,lam_H5,
                            lam_H1,lam_H2,lam_H3,lam_H4,lam_H5)

#### other parameters
#### jakarta pop structure
rho <- 0.15
gamma <- 0.35
amin <- c(0,5,10,15,20,45,55,65,75) # lower bound of age groups
amax <- c(4,9,14,19,44,54,64,74,99) # upper bound of age groups
pop <- matrix(round(pop_jakarta_sim_vec * 3000000),nrow=10,ncol=9,byrow=TRUE) # population
hist_length <- 80

#### simulating using jakarta pop structures
jkt_cases_sim_list <- list()
for (i in 1:11){
  
  jkt_cases_sim_list[[i]] <- simcases(nT=length(lam_curr_scenario_list[[i]]), 
                                      nA=length(amin), 
                                      lamH=lam_H_scenario_list[[i]], 
                                      lam=lam_curr_scenario_list[[i]], 
                                      rho=rho, gamma=gamma, 
                                      pop=pop, 
                                      amin=amin, 
                                      amax=amax)
  
}

#### load stan model - poisson likelihood
catalytic_mod1 <- stan_model(file="codes/stan/case_age_catalytic_PS.stan")
catalytic_mod2 <- stan_model(file="codes/stan/case_age_catalytic_multinomial_PS.stan")
catalytic_mod3 <- stan_model(file="codes/stan/case_age_catalytic_multi_lamH_PS.stan")
catalytic_mod4 <- stan_model(file="codes/stan/case_age_catalytic_multinomial_multi_lamH_PS.stan")

#### prep data for stan fitting for each scenario
data_pois_scenario <- list()
data_multi_scenario <- list()
for (i in 1:11){
  
  data_pois_scenario[[i]] <- list(nA=length(amin), nT=length(lam_curr_scenario_list[[i]]), 
                                  cases=jkt_cases_sim_list[[i]]$cases_PS, 
                                  pop=pop, 
                                  ageLims=rbind(amin,amax),
                                  hist_length=length(lam_H_scenario_list[[i]]))
  
  data_multi_scenario[[i]] <- list(nA=length(amin), nT=length(lam_curr_scenario_list[[i]]), 
                                   cases=jkt_cases_sim_list[[i]]$cases_PS, 
                                   totalcases=jkt_cases_sim_list[[i]]$totalcases_PS, 
                                   pop=pop, 
                                   ageLims=rbind(amin,amax),
                                   hist_length=length(lam_H_scenario_list[[i]]))
  
}

#### fitting all scenarios based on different models
params1 <- c("lam_H","lam_t","rho","gamma")
params2 <- c("lam_H","lam_t","gamma")

stanfit_mod1_scenario_list <- list()
stanfit_mod2_scenario_list <- list()
stanfit_mod3_scenario_list <- list()
stanfit_mod4_scenario_list <- list()
for (i in 1:11){
  
  # model 1: poisson
  stanfit_mod1_scenario_list[[i]] <- sampling(catalytic_mod1, 
                                              data=data_pois_scenario[[i]], 
                                              warmup = 500, 
                                              iter = 2500, 
                                              chains = 4,
                                              pars = params1,
                                              seed = 12345)
  
}
stanfit_mod1_scenario_list[[1]] <- sampling(catalytic_mod1, 
                                            data=data_pois_scenario[[1]], 
                                            warmup = 1000, 
                                            iter = 5000, 
                                            chains = 4,
                                            pars = params1,
                                            seed = 123)
saveRDS(stanfit_mod1_scenario_list,"output/stanfit_mod1_scenario_list.rds")

for (i in 1:11){
  
  # model 2: multinomial
  stanfit_mod2_scenario_list[[i]] <- sampling(catalytic_mod2, 
                                              data=data_multi_scenario[[i]], 
                                              warmup = 500, 
                                              iter = 2500, 
                                              chains = 4,
                                              pars = params2,
                                              seed = 12345)
  
}
stanfit_mod2_scenario_list[[1]] <- sampling(catalytic_mod2, 
                                            data=data_multi_scenario[[1]], 
                                            warmup = 1000, 
                                            iter = 5000, 
                                            chains = 4,
                                            pars = params2,
                                            seed = 8799)
saveRDS(stanfit_mod2_scenario_list,"output/stanfit_mod2_scenario_list.rds")

for (i in 1:11){
  
  # model 3: poisson
  stanfit_mod3_scenario_list[[i]] <- sampling(catalytic_mod3, 
                                              data=data_pois_scenario[[i]], 
                                              warmup = 500, 
                                              iter = 2500, 
                                              chains = 4,
                                              pars = params1,
                                              seed = 12345)
  
}
saveRDS(stanfit_mod3_scenario_list,"output/stanfit_mod3_scenario_list.rds")

for (i in 1:11){
  
  # model 4: multinomial
  stanfit_mod4_scenario_list[[i]] <- sampling(catalytic_mod4, 
                                              data=data_multi_scenario[[i]], 
                                              warmup = 500, 
                                              iter = 2500, 
                                              chains = 4,
                                              pars = params2,
                                              seed = 12345)
  
}
saveRDS(stanfit_mod4_scenario_list,"output/stanfit_mod4_scenario_list.rds")

#### post-processing and fit check
posterior_samples_mod1_scenario_list <- list()
summary_mod1_scenario_list <- list()
simulated_cases_mod1_scenario_list <- list()
plot_summary_mod1_scenario_list <- list()

posterior_samples_mod2_scenario_list <- list()
summary_mod2_scenario_list <- list()
simulated_cases_mod2_scenario_list <- list()
plot_summary_mod2_scenario_list <- list()

posterior_samples_mod3_scenario_list <- list()
summary_mod3_scenario_list <- list()
simulated_cases_mod3_scenario_list <- list()
plot_summary_mod3_scenario_list <- list()

posterior_samples_mod4_scenario_list <- list()
summary_mod4_scenario_list <- list()
simulated_cases_mod4_scenario_list <- list()
plot_summary_mod4_scenario_list <- list()

plot_compare_scenario_list <- list()

sample_sizes <- 1000

for (i in 1:11){
  
  # model 1: poisson
  posterior_samples_mod1_scenario_list[[i]] <- 
    posterior_samples_output(stanfit_mod1_scenario_list[[i]],params1,sample_sizes)
  summary_mod1_scenario_list[[i]] <- 
    summarise_posterior(stanfit_mod1_scenario_list[[i]],
                        params=params1,
                        original_values = list(lam_H=lam_H_scenario_list[[i]],
                                               lam_t=lam_curr_scenario_list[[i]],
                                               rho=rho,
                                               gamma=gamma))
  simulated_cases_mod1_scenario_list[[i]] <- 
    simulate_cases_from_posterior(posterior_samples_mod1_scenario_list[[i]],
                                  jkt_cases_sim_list[[i]]$cases_PS,
                                  pop=pop,
                                  amin=amin,
                                  amax=amax,
                                  hist_length=hist_length,
                                  model_type="cases_PS",
                                  dist="Poisson")
  plot_summary_mod1_scenario_list[[i]] <- 
    plot_model_summary(summary_mod1_scenario_list[[i]],simulated_cases_mod1_scenario_list[[i]],
                       dist="Poisson")
  
  # model 2: multinomial
  posterior_samples_mod2_scenario_list[[i]] <- 
    posterior_samples_output(stanfit_mod2_scenario_list[[i]],params2,sample_sizes)
  summary_mod2_scenario_list[[i]] <- 
    summarise_posterior(stanfit_mod2_scenario_list[[i]],
                        params=params2,
                        original_values = list(lam_H=lam_H_scenario_list[[i]],
                                               lam_t=lam_curr_scenario_list[[i]],
                                               gamma=gamma))
  simulated_cases_mod2_scenario_list[[i]] <- 
    simulate_cases_from_posterior(posterior_samples_mod2_scenario_list[[i]],
                                  jkt_cases_sim_list[[i]]$cases_PS,
                                  pop=pop,
                                  amin=amin,
                                  amax=amax,
                                  hist_length=hist_length,
                                  model_type="cases_PS",
                                  dist="multinomial")
  plot_summary_mod2_scenario_list[[i]] <- 
    plot_model_summary(summary_mod2_scenario_list[[i]],simulated_cases_mod2_scenario_list[[i]],
                       dist="multinomial")
  
  # model 3: Poisson with multiple lam_H
  posterior_samples_mod3_scenario_list[[i]] <- 
    posterior_samples_output_new(stanfit_mod3_scenario_list[[i]],params1,sample_sizes)
  summary_mod3_scenario_list[[i]] <- 
    summarise_posterior_new(stanfit_mod3_scenario_list[[i]],
                            params=params1,
                            original_values = list(lam_H=lam_H_scenario_list[[i]],
                                                   lam_t=lam_curr_scenario_list[[i]],
                                                   rho=rho,
                                                   gamma=gamma))
  simulated_cases_mod3_scenario_list[[i]] <- 
    simulate_cases_from_posterior_new(posterior_samples_mod3_scenario_list[[i]],
                                      jkt_cases_sim_list[[i]]$cases_PS,
                                      pop=pop,
                                      amin=amin,
                                      amax=amax,
                                      hist_length=hist_length,
                                      model_type="cases_PS",
                                      dist="Poisson")
  
  plot_summary_mod3_scenario_list[[i]] <- 
    plot_model_summary(summary_mod3_scenario_list[[i]],simulated_cases_mod3_scenario_list[[i]],
                       dist="Poisson")
  
  # model 4: multinomial with multiple lam_H
  posterior_samples_mod4_scenario_list[[i]] <- 
    posterior_samples_output_new(stanfit_mod4_scenario_list[[i]],params2,sample_sizes)
  summary_mod4_scenario_list[[i]] <- 
    summarise_posterior_new(stanfit_mod4_scenario_list[[i]],
                            params=params2,
                            original_values = list(lam_H=lam_H_scenario_list[[i]],
                                                   lam_t=lam_curr_scenario_list[[i]],
                                                   gamma=gamma))
  simulated_cases_mod4_scenario_list[[i]] <- 
    simulate_cases_from_posterior_new(posterior_samples_mod4_scenario_list[[i]],
                                      jkt_cases_sim_list[[i]]$cases_PS,
                                      pop=pop,
                                      amin=amin,
                                      amax=amax,
                                      hist_length=hist_length,
                                      model_type="cases_PS",
                                      dist="multinomial")
  
  plot_summary_mod4_scenario_list[[i]] <- 
    plot_model_summary(summary_mod4_scenario_list[[i]],simulated_cases_mod4_scenario_list[[i]],
                       dist="multinomial")
  
  # compare models
  # plot_compare_scenario_list[[i]] <- 
  #   compare_estimates_new(list(summary_mod1_scenario_list[[i]]$summary,
  #                              summary_mod2_scenario_list[[i]]$summary,
  #                              summary_mod3_scenario_list[[i]]$summary,
  #                              summary_mod4_scenario_list[[i]]$summary),
  #                         c("Poisson 1","multinomial 1","Poisson 2","multinomial 2"),
  #                         params1,80)
  
}

for (i in 1:11){
  # compare models
  plot_compare_scenario_list[[i]] <- 
    compare_estimates_new(list(summary_mod1_scenario_list[[i]]$summary,
                               summary_mod2_scenario_list[[i]]$summary,
                               summary_mod3_scenario_list[[i]]$summary,
                               summary_mod4_scenario_list[[i]]$summary),
                          c("Poisson 1","multinomial 1","Poisson 2","multinomial 2"),
                          params1,
                          list(lam_H=lam_H_scenario_list[[i]],
                               lam_t=lam_curr_scenario_list[[i]],
                               rho=rho,
                               gamma=gamma),
                          80)
  
  ggsave(paste0("output/plot_combined_scenario_",i,".jpg"),
         plot_compare_scenario_list[[i]]$plot_all,height=14,width=25,unit="cm")
}

# for (i in 1:11){
#   
#   plot_combined <- plot_grid(plot_summary_mod1_scenario_list[[i]],
#                              plot_summary_mod2_scenario_list[[i]],
#                              labels=c("a","b"))
#   ggsave(paste0("output/plot_combined_mod1_mod2_scenario_",i,".jpg"),
#          plot_combined,height=14,width=25,unit="cm")
#   
# }
