#### read admin data
admin2_list <- read_csv("data/admin2_list.csv") %>% 
  clean_names() %>%
  mutate(region=replace(region, region=="SUMATERA", "SUMATRA"))

#### prepare data for simulation study
#### read pop data
#### 2010 & 2020 population data based on censuses
#### other years data from linear interpolation and extrapoliation
pop_census_2020 <- readRDS("data/pop_age_simplegrowth_2009_2023.rds") %>% 
  bind_rows() %>% 
  clean_names() %>% 
  left_join(admin2_list) %>% 
  dplyr::select(region,idadmin1,idadmin2,admin1,admin2,everything()) %>% 
  filter(year==2020)

#### using jakarta pop structure
pop_jakarta <- pop_census_2020 %>% 
  filter(admin1=="DKI JAKARTA") %>% 
  group_by(admin1,age_group) %>% 
  summarise(pop=sum(pop)) %>% 
  ungroup()

pop_age_census <- as.character(pop_jakarta$age_group)
pop_age_jkt <- c("0-4","5-9","10-14","15-19","20-44","20-44","20-44","20-44","20-44",
                 "45-54","45-54","55-64","55-64","65-74","65-74","75-99")
pop_age_jkt_group <- unique(pop_age_jkt)

pop_age_jkt_df <- tibble(age_group=pop_age_census,age_group2=pop_age_jkt)

pop_jakarta_sim <- pop_jakarta %>% 
  left_join(pop_age_jkt_df) %>% 
  group_by(age_group2) %>% 
  summarise(pop=sum(pop)) %>% 
  ungroup() %>% 
  mutate(prop=pop/sum(pop)) %>% 
  select(age_group=age_group2,pop,prop) %>% 
  mutate(age_group=factor(age_group,levels=unique(pop_age_jkt))) %>% 
  arrange(age_group)

pop_jakarta_sim_vec <- pop_jakarta_sim$prop

#### prepare data for real data fitting
pop_jkt_2017_2023 <- readRDS("data/pop_age_simplegrowth_2009_2023.rds") %>% 
  bind_rows() %>% 
  clean_names() %>% 
  left_join(admin2_list) %>% 
  dplyr::select(region,idadmin1,idadmin2,admin1,admin2,everything()) %>% 
  filter(year %in% 2017:2023) %>% 
  filter(admin1 == "DKI JAKARTA") %>% 
  left_join(pop_age_jkt_df) %>% 
  group_by(idadmin2,admin2,year,age_group2) %>% 
  summarise(pop=sum(pop)) %>% 
  ungroup() %>% 
  mutate(age_group=factor(age_group2,levels=pop_age_jkt_group)) %>% 
  select(idadmin2,admin2,year,age_group,pop) %>% 
  arrange(idadmin2,year,age_group) %>% 
  pivot_wider(names_from=age_group,values_from=pop)

jkt_admin2 <- unique(pop_jkt_2017_2023$admin2)

pop_jkt_admin2_2017_2023_list <- list()
for(i in seq_len(length(jkt_admin2))){
  
  pop_data_jkt <- pop_jkt_2017_2023 %>% filter(admin2==jkt_admin2[i])
  pop_jkt_admin2_2017_2023_list[[i]] <- as.matrix(pop_data_jkt[,-(1:3)])
  
}

#### pop total jkt
pop_total_jkt_2017_2023 <- pop_jkt_2017_2023 %>% 
  mutate(pop=`0-4`+`5-9`+`10-14`+`15-19`+`20-44`+`45-54`+`55-64`+`65-74`+`75-99`) %>% 
  select(admin2,year,pop)
