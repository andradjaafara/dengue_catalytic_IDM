library(tidyverse)
library(janitor)
library(rstan)
library(cowplot)
library(bayesplot)
source("codes/functions.R")

options(mc.cores = parallel::detectCores())

#### load stan model
catalytic_mod6 <- stan_model(file="codes/stan/case_age_catalytic_multi_lamH_PS_pooling.stan")

#### load data
dengue_jkt_list <- readRDS("data/dengue_jkt_list.rds")
pop_jkt_list <- readRDS("data/pop_jkt_list.rds")

jkt_admin2_v2 <- c("KOTA JAKARTA PUSAT","KOTA JAKARTA UTARA","KOTA JAKARTA BARAT",
                   "KOTA JAKARTA SELATAN","KOTA JAKARTA TIMUR","KEPULAUAN SERIBU")

#### prep stan data
cases_jkt <- aperm(simplify2array(dengue_jkt_list),c(3,2,1))
pop_jkt <- aperm(simplify2array(pop_jkt_list),c(3,2,1))
data_jkt_prov <- list(nA=dim(cases_jkt)[3], nT=dim(cases_jkt)[2], nD=dim(cases_jkt)[1],
                      cases=cases_jkt, 
                      pop=pop_jkt, 
                      ageLims=rbind(amin,amax),
                      hist_length=60)

params_prov <- c("lam_H_prov","kappa_lam_H","lam_H",
                 "lam_t_prov","kappa_lam_t","lam_t",
                 "rho_prov","kappa_rho","rho",
                 "gamma_prov","kappa_gamma","gamma")

#### stan fitting
set.seed(87987987)
seed_rand <- round(runif(1,1,1000000),0)
stanfit_jkt_mod6_pooling <- sampling(catalytic_mod6, 
                                     data=data_jkt_prov, 
                                     warmup = 800, 
                                     iter = 8000, 
                                     chains = 4,
                                     pars = params_prov,
                                     seed = seed_rand)

View(round(summary(stanfit_jkt_mod6_pooling)$summary,3))

traceplot(stanfit_jkt_mod6_pooling,pars="lam_t")
traceplot(stanfit_jkt_mod6_pooling,pars="lam_t_prov")
traceplot(stanfit_jkt_mod6_pooling,pars="rho")
traceplot(stanfit_jkt_mod6_pooling,pars="kappa_rho")
traceplot(stanfit_jkt_mod6_pooling,pars="rho_prov")
traceplot(stanfit_jkt_mod6_pooling,pars="gamma")
traceplot(stanfit_jkt_mod6_pooling,pars="kappa_gamma")
traceplot(stanfit_jkt_mod6_pooling,pars="gamma_prov")
traceplot(stanfit_jkt_mod6_pooling,pars="phi")

#### summarising posterior
stanfit_jkt_mod6_pooling_samples <- extract(stanfit_jkt_mod6_pooling,pars=params_prov)
n_samples <- dim(stanfit_jkt_mod6_pooling_samples$lam_H_prov)[1]



colMeans(stanfit_jkt_mod6_pooling_samples$gamma)
colMeans(stanfit_jkt_mod6_pooling_samples$rho)

no_of_admin2 <- length(jkt_admin2_v2)

#### lam_H
lam_H_prov <- stanfit_jkt_mod6_pooling_samples$lam_H_prov
lam_H_jkt_list <- list()
for (i in 1:no_of_admin2){
  lam_H_jkt_list[[i]] <- as_tibble(stanfit_jkt_mod6_pooling_samples$lam_H[,i,]) 
  colnames(lam_H_jkt_list[[i]]) <- 2017 + -1:-60
  lam_H_jkt_list[[i]] <- lam_H_jkt_list[[i]] %>% 
    mutate(admin2=jkt_admin2_v2[i],
           rep=1:n_samples) %>% 
    select(admin2,rep,everything()) %>% 
    pivot_longer(-c(admin2,rep),names_to="t",values_to="lam_H")
}

#### lam_t
lam_t_prov <- stanfit_jkt_mod6_pooling_samples$lam_t_prov
lam_t_jkt_list <- list()
for (i in 1:no_of_admin2){
  lam_t_jkt_list[[i]] <- as_tibble(stanfit_jkt_mod6_pooling_samples$lam_t[,i,]) 
  colnames(lam_t_jkt_list[[i]]) <- 2017:2023
  lam_t_jkt_list[[i]] <- lam_t_jkt_list[[i]] %>% 
    mutate(admin2=jkt_admin2_v2[i],
           rep=1:n_samples) %>% 
    select(admin2,rep,everything()) %>% 
    pivot_longer(-c(admin2,rep),names_to="t",values_to="lam_t")
}

#### rho
rho_prov <- stanfit_jkt_mod6_pooling_samples$rho_prov
rho_jkt_list <- list()
for (i in 1:no_of_admin2){
  rho_jkt_list[[i]] <- tibble(admin2=jkt_admin2_v2[i],t=1:n_samples,
                              rho=stanfit_jkt_mod6_pooling_samples$rho[,i])
}

#### gamma
gamma_prov <- stanfit_jkt_mod6_pooling_samples$gamma_prov
gamma_jkt_list <- list()
for (i in 1:no_of_admin2){
  gamma_jkt_list[[i]] <- tibble(admin2=jkt_admin2_v2[i],t=1:n_samples,
                                gamma=stanfit_jkt_mod6_pooling_samples$gamma[,i])
}

lam_H_jkt <- bind_rows(lam_H_jkt_list) %>% 
  group_by(admin2,t) %>% 
  summarise(median=median(lam_H,na.rm=TRUE),
            lo1=quantile(lam_H,probs=0.025,na.rm=TRUE),
            up1=quantile(lam_H,probs=0.975,na.rm=TRUE),
            lo2=quantile(lam_H,probs=0.25,na.rm=TRUE),
            up2=quantile(lam_H,probs=0.75,na.rm=TRUE)) %>% 
  mutate(t=as.numeric(t)) %>% 
  arrange(admin2,t)
lam_t_jkt <- bind_rows(lam_t_jkt_list) %>% 
  group_by(admin2,t) %>% 
  summarise(median=median(lam_t,na.rm=TRUE),
            lo1=quantile(lam_t,probs=0.025,na.rm=TRUE),
            up1=quantile(lam_t,probs=0.975,na.rm=TRUE),
            lo2=quantile(lam_t,probs=0.25,na.rm=TRUE),
            up2=quantile(lam_t,probs=0.75,na.rm=TRUE)) %>% 
  mutate(t=as.numeric(t)) %>% 
  arrange(admin2,t)
rho_jkt <- bind_rows(rho_jkt_list) %>% 
  group_by(admin2) %>% 
  summarise(median=median(rho,na.rm=TRUE),
            lo1=quantile(rho,probs=0.025,na.rm=TRUE),
            up1=quantile(rho,probs=0.975,na.rm=TRUE),
            lo2=quantile(rho,probs=0.25,na.rm=TRUE),
            up2=quantile(rho,probs=0.75,na.rm=TRUE))
gamma_jkt <- bind_rows(gamma_jkt_list) %>% 
  group_by(admin2) %>% 
  summarise(median=median(gamma,na.rm=TRUE),
            lo1=quantile(gamma,probs=0.025,na.rm=TRUE),
            up1=quantile(gamma,probs=0.975,na.rm=TRUE),
            lo2=quantile(gamma,probs=0.25,na.rm=TRUE),
            up2=quantile(gamma,probs=0.75,na.rm=TRUE))

lam_all_jkt <- bind_rows(lam_H_jkt,lam_t_jkt) %>% 
  arrange(admin2,t) %>% 
  ungroup()

lambda_estim <- readRDS("lambda_estim.rds") %>% 
  filter(study %in% c("2_Kali Deres","2_Pesanggrahan","2_Pulo Gadung")) %>% 
  mutate(admin2=c("KOTA JAKARTA BARAT","KOTA JAKARTA SELATAN","KOTA JAKARTA TIMUR"),
         t=2014)

plot_lam_all_jkt <- lam_all_jkt %>% 
  ggplot(aes(x=t,group=admin2,col=admin2,fill=admin2)) +
  theme_classic(base_size=24) +
  theme(panel.grid.major.y = element_line(linewidth = 0.25,colour="black")) +
  geom_line(aes(y=4*median)) +
  geom_ribbon(aes(ymin=4*lo1,ymax=4*up1),alpha=0.25) +
  geom_ribbon(aes(ymin=4*lo2,ymax=4*up2),alpha=0.50) +
  geom_point(data=lambda_estim,aes(y=lambda_median),size=5) +
  geom_linerange(data=lambda_estim,aes(ymin=lambda_lo,ymax=lambda_up),linewidth=1) +
  labs(x="Year",y=expression(4~lambda)) +
  expand_limits(y=0) +
  facet_grid(admin2~.) +
  geom_vline(aes(xintercept=2016.5),linetype=2) +
  theme(strip.background = element_blank(),strip.text = element_blank(),
        legend.position = "none") +
  scale_colour_grafify(palette="vibrant") +
  scale_fill_grafify(palette="vibrant")

plot_rho_jkt <- rho_jkt %>% 
  ggplot(aes(x=admin2,group=admin2,col=admin2)) +
  geom_point(aes(y=median),size=5,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lo1,ymax=up1),linewidth=1,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lo2,ymax=up2),linewidth=2,position=position_dodge(width=0.5)) +
  theme_classic(base_size=24) +
  labs(x=NULL,y=expression(rho)) +
  expand_limits(y=0) +
  scale_colour_grafify(palette="vibrant") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

plot_gamma_jkt <- gamma_jkt %>% 
  ggplot(aes(x=admin2,group=admin2,col=admin2)) +
  geom_point(aes(y=median),size=5,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lo1,ymax=up1),linewidth=1,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lo2,ymax=up2),linewidth=2,position=position_dodge(width=0.5)) +
  theme_classic(base_size=24) +
  labs(x=NULL,y=expression(gamma),col=NULL) +
  expand_limits(y=0) +
  scale_colour_grafify(palette="vibrant") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position=c(0.7,0.8))

plot_1_jkt <- plot_grid(plot_rho_jkt,plot_gamma_jkt)
plot_2_jkt <- plot_grid(plot_lam_all_jkt,plot_1_jkt,ncol=1,rel_heights = c(1,0.4))

ggsave("plot_fitting.png",plot_2_jkt,width=37.8,height=25,units="cm")

#### visualise
plot_lam_H_jkt <- lam_H_jkt %>% 
  ggplot(aes(x=t,group=admin2,col=admin2,fill=admin2)) +
  geom_line(aes(y=4*median)) +
  geom_ribbon(aes(ymin=4*lo1,ymax=4*up1),alpha=0.25) +
  geom_ribbon(aes(ymin=4*lo2,ymax=4*up2),alpha=0.50) +
  theme_classic(base_size=12) +
  labs(x="Year",y="4*lam_H") +
  expand_limits(y=0)

plot_lam_t_jkt <- lam_t_jkt %>% 
  ggplot(aes(x=t,group=admin2,col=admin2)) +
  geom_point(aes(y=4*median),size=2.5,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=4*lo1,ymax=4*up1),linewidth=1,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=4*lo2,ymax=4*up2),linewidth=2,position=position_dodge(width=0.5)) +
  theme_classic(base_size=12) +
  labs(x="Year",y="4*lam_t") +
  expand_limits(y=0)

plot_rho_jkt <- rho_jkt %>% 
  ggplot(aes(x=admin2,group=admin2,col=admin2)) +
  geom_point(aes(y=median),size=2.5,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lo1,ymax=up1),linewidth=1,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lo2,ymax=up2),linewidth=2,position=position_dodge(width=0.5)) +
  theme_classic(base_size=12) +
  labs(x="Year",y="rho") +
  expand_limits(y=0)

plot_gamma_jkt <- gamma_jkt %>% 
  ggplot(aes(x=admin2,group=admin2,col=admin2)) +
  geom_point(aes(y=median),size=2.5,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lo1,ymax=up1),linewidth=1,position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=lo2,ymax=up2),linewidth=2,position=position_dodge(width=0.5)) +
  theme_classic(base_size=12) +
  labs(x="Year",y="gamma") +
  expand_limits(y=0)

plot_jkt_1 <- plot_grid(plot_lam_H_jkt,plot_lam_t_jkt,nrow=2)
plot_jkt_2 <- plot_grid(plot_rho_jkt,plot_gamma_jkt,nrow=1)
plot_jkt_3 <- plot_grid(plot_jkt_1,plot_jkt_2,nrow=2)

#### data and model comparisons
sample_idx <- sample(1:n_samples,1000)
nT_jkt_mod3 <- ncol(stanfit_jkt_mod6_pooling_samples[["lam_t"]])
nA_jkt_mod3 <- length(amin)
age_groups <- paste(amin,amax,sep='-')

simulated_cases_jkt_mod3_list <- list()
for (i in seq_len(ncol(stanfit_jkt_mod6_pooling_samples[["gamma"]]))){
  simulated_cases_jkt_mod3_iter_list <- list()
  for (j in seq_len(length(sample_idx))){
    
    lamH <- stanfit_jkt_mod6_pooling_samples[["lam_H"]][j,i,]
    lam <- stanfit_jkt_mod6_pooling_samples[["lam_t"]][j,i,]
    rho <- stanfit_jkt_mod6_pooling_samples[["rho"]][j,i]
    gamma <- stanfit_jkt_mod6_pooling_samples[["gamma"]][j,i]
    phi <- stanfit_jkt_mod6_pooling_samples[["phi"]][j]
    pop_jkt_mod3 <- t(pop_jkt_list[[i]])
    simulated_output <- simcases_sto(nT=nT_jkt_mod3, nA=nA_jkt_mod3, lamH=lamH, lam=lam, 
                                     rho=rho, gamma=gamma, phi=phi, dist="Poisson",
                                     pop=pop_jkt_mod3, amin=amin, amax=amax)
    simulated_cases_df <- simulated_output[["cases_PS"]]
    colnames(simulated_cases_df) <- age_groups
    simulated_cases_jkt_mod3_iter_list[[j]] <- simulated_cases_df %>% 
      as_tibble() %>% 
      mutate(admin2=jkt_admin2_v2[i],t=1:nT_jkt_mod3,iter=j) %>% 
      pivot_longer(-c(t,iter,admin2),names_to="age_group",values_to="simulated_cases") %>% 
      mutate(age_group=factor(age_group,levels=age_groups))
    
  }
  simulated_cases_jkt_mod3_list[[i]] <- bind_rows(simulated_cases_jkt_mod3_iter_list)
}

simulated_cases_jkt_mod3 <- bind_rows(simulated_cases_jkt_mod3_list)
simulated_cases_jkt_mod3_summary <- simulated_cases_jkt_mod3 %>% 
  group_by(admin2,t,age_group) %>% 
  summarise(median=median(simulated_cases,na.rm=TRUE),
            lo1=quantile(simulated_cases,probs=0.025,na.rm=TRUE),
            lo2=quantile(simulated_cases,probs=0.25,na.rm=TRUE),
            up1=quantile(simulated_cases,probs=0.975,na.rm=TRUE),
            up2=quantile(simulated_cases,probs=0.75,na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(year=t+2016) %>% 
  select(admin2,t,year,age_group,everything())

reported_cases_jkt_list <- list()
for (i in seq_len(length(jkt_admin2_v2))){
  
  reported_cases_jkt_list[[i]] <- as_tibble(dengue_jkt_list[[i]]) %>% 
    mutate(age_group=factor(age_groups,levels=age_groups),admin2=jkt_admin2_v2[i]) %>% 
    pivot_longer(-c(age_group,admin2),names_to="year",values_to="reported_cases") %>% 
    mutate(year=as.numeric(year)) %>% 
    select(admin2,year,age_group,reported_cases) %>% 
    arrange(admin2,year,age_group)
  
}
reported_cases_jkt <- bind_rows(reported_cases_jkt_list)

simulated_cases_jkt_mod3_summary <- simulated_cases_jkt_mod3_summary %>% 
  left_join(reported_cases_jkt)

plot_compare_jkt_mod3_list <- list()
for (i in seq_len(length(jkt_admin2_v2))){
  
  plot_compare_jkt_mod3_list[[i]] <- simulated_cases_jkt_mod3_summary %>% 
    filter(admin2==jkt_admin2_v2[i]) %>% 
    ggplot(aes(x=t+2016)) +
    geom_linerange(aes(ymin=lo1,ymax=up1),col="#be0032",linewidth=1) +
    geom_linerange(aes(ymin=lo2,ymax=up2),col="#be0032",linewidth=2) +
    geom_point(aes(y=median),col="#be0032",size=1.5) +
    geom_point(aes(y=reported_cases),size=1.5) +
    theme_classic(base_size=12) +
    facet_grid(age_group~.) +
    scale_x_continuous(breaks=1:nT_jkt_mod3) +
    labs(x="Year",y="Cases")
  
}

plot_fitting_data_jkt <- simulated_cases_jkt_mod3_summary %>% 
  filter(admin2 %in% c("KOTA JAKARTA TIMUR","KEPULAUAN SERIBU")) %>%
  # pivot_longer(-c(admin2,t,year,age_group,lo1,lo2,up1,up2),names_to="var",values="val") +
  ggplot(aes(x=t+2016)) +
  geom_linerange(aes(ymin=lo1,ymax=up1,col="Modelled"),linewidth=1) +
  geom_linerange(aes(ymin=lo2,ymax=up2,col="Modelled"),linewidth=2) +
  geom_point(aes(y=median,col="Modelled"),size=5) +
  geom_point(aes(y=reported_cases,col="Reported"),size=5) +
  scale_colour_manual(breaks=c("Modelled","Reported"),values=c("#be0032","black")) +
  theme_classic(base_size=24) +
  facet_grid(admin2~age_group,scales = "free_y") +
  scale_x_continuous(breaks=(1:nT_jkt_mod3)+2016) +
  labs(x=NULL,y="Dengue cases",colour=NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.925,0.875),
        strip.background.y = element_blank(),
        strip.text.y = element_blank())

ggsave("output/plot_fitting_data.png",plot_fitting_data_jkt,width=37.8,height=13.46,units="cm")


