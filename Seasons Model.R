### Just keep building up the complexity.
# Now I have a temperature*depth interaction that differs for each of the seasons.
# Next add main effect of distance to shore that differs across the seasons.
# Is the ultimate goal to build up to a three way interaction that differs over seasons? 

library(rstan)
library(tidyverse)
library(lubridate)
library(matrixStats)
library(patchwork)
#library(cmdstanr)



########### Trialling model ############

## Combine deltatemp_data with ssf data.
# head(ssf_data$t2_h)
# 
# ssf_data<-left_join(ssf_data,deltatemp_data,by=c("t2_h"="DateTimeUTC"))

temp<-readRDS("available_temp_2020-2021_clean.RDS")

temp_range<-temp%>%filter(depth_m ==0 | depth_m == 10)
temp_range_daily<-temp_range%>%group_by(DateTimeUTC,depth_m)%>%summarise(temp_c=mean(temp_c))

temp_range_daily_w<-temp_range_daily%>%pivot_wider(names_from = depth_m, values_from = temp_c,id_cols=c(DateTimeUTC))

deltatemp_data<-temp_range_daily_w%>%mutate(delta_temp=`0`-`10`)

deltatemp_data$month_num<-month(deltatemp_data$DateTimeUTC)
deltatemp_data$year<-year(deltatemp_data$DateTimeUTC)

## Creating the homogenous and heterogenous layers
deltatemp_data_sum<-deltatemp_data%>%group_by(month_num,year)%>%summarise(mean(`0`),mean(`10`,na.rm=T),meanDT=mean(delta_temp,na.rm=T))
#deltatemp_data_sum$season<-ifelse(deltatemp_data_sum$meanDT<1.5, "cold", 
                                  #ifelse(deltatemp_data_sum$meanDT>2 ,"hetero", "hot"))

# if greater than 3.5 then hetero. If less then 2.5 then cold homo

deltatemp_data_sum$season<-NA
deltatemp_data_sum$season<-ifelse(deltatemp_data_sum$meanDT>2.5, "hetero", "cold")

# Combine delta temp summary data
ssf_data<-readRDS("10randomstps_withtempdist.RDS")
ssf_data$month_num <-month(ssf_data$t2_)
ssf_data$year<-year(ssf_data$t2_)

ssf_data<-left_join(ssf_data,deltatemp_data_sum,by=c("month_num", "year"))

ssf_data<-ssf_data%>%ungroup()%>%mutate(st.WT=base::scale(watertemp_),st.WT2=base::scale(watertemp_^2),st.depth=base::scale(depth_m),st.depth2=base::scale(depth_m^2),st.distshore=base::scale(dist_to_shore_), logsl_=log(sl_),cos_ta_=cos(ta_)) #,year=year(t2_)) #,st.deltaT=scale(delta_temp))

ssf_data<-drop_na(ssf_data,cos_ta_)

ssf_data$season_num<-as.numeric(as.factor(ssf_data$season))

model <- "
data {
  int<lower=1> N; // no data points
  int<lower=1> I; // no steps (over all individuals)
  int<lower=1> J; // no individuals
  int<lower=1> K; // number of predictors

  
  int<lower=0> y[N]; // response
  matrix[N,K]  X;
  int<lower=1, upper=I> stepid[N]; // step id
  int<lower=1, upper=J> indid[N]; // individual id
  int<lower=1,upper=3> season[N];
}

parameters {
  vector[3] beta_temp;
  vector[3] beta_depth;
  vector[3] beta_dist;
  vector[3] beta_td_int;
  vector[3] beta_tdist_int;
  vector[3] beta_ddist_int;
  vector[3] beta_threeway_int;
  vector[3] beta_sl;
  vector[3] beta_ta;
  vector<lower=0>[9] sigma;
  vector[I] a_re; // RE for steps
  vector[J] i_re1; // RE effects for each individual (temp)
  vector[J] i_re2; // RE effects for each individual (deltaT)
  vector[J] i_re3; // RE effects for each individual (interaction)
  vector[J] i_re4; // RE effects for each individual (interaction)
  vector[J] i_re5; // RE effects for each individual (interaction)
  vector[J] i_re6;
  vector[J] i_re7;
  vector[J] i_re8;
  vector[J] i_re9;
}

model {
  vector[N] mu;
  
  // priors
  sigma ~ normal(0, 1);
  a_re ~ normal(0, 1); // The paper has this as a_re ~ normal(0, 1000000)
  i_re1 ~ normal(0, sqrt(sigma[1]));
  i_re2 ~ normal(0, sqrt(sigma[2]));
  i_re3 ~ normal(0, sqrt(sigma[3]));
  i_re4 ~ normal(0, sqrt(sigma[4]));
  i_re5 ~ normal(0, sqrt(sigma[5]));
  i_re6 ~ normal(0, sqrt(sigma[6]));
  i_re7 ~ normal(0, sqrt(sigma[7]));
  i_re8 ~ normal(0, sqrt(sigma[8]));
  i_re9 ~ normal(0, sqrt(sigma[9]));
  beta_temp ~ normal(0,2);
  beta_depth ~ normal(0,2);
  beta_dist ~ normal(0,2);
  beta_td_int ~ normal(0,2);
  beta_tdist_int ~ normal(0,2);
  beta_ddist_int ~ normal(0,2);
  beta_threeway_int ~ normal(0,2);
  beta_sl ~ normal(0,2);
  beta_ta ~ normal(0,2);

  //likelihood
  mu = a_re[stepid]+
    (beta_temp[season] + i_re1[indid]) .* X[,1]+ //temp effect 
    (beta_depth[season] + i_re2[indid]) .* X[,2] + // depth effect
    (beta_dist[season] + i_re3[indid]) .* X[,3] + // distance to shore effect
    (beta_td_int[season] + i_re4[indid]) .* X[,1] .* X[,2] + //temperature*depth
    (beta_tdist_int[season] + i_re5[indid]) .* X[,1] .* X[,3] + //temperature*distance
    (beta_ddist_int[season] + i_re6[indid]) .* X[,2] .* X[,3] + //distance*depth
    (beta_threeway_int[season] + i_re7[indid]) .* X[,1] .* X[,2] .* X[,3]+
    (beta_sl[season] + i_re8[indid]) .* X[,4] +
    (beta_ta[season] + i_re9[indid]) .* X[,5];
  y ~ poisson_log(mu);
}
"

X <- model.matrix(y ~ st.WT +st.depth + st.distshore + logsl_ + cos_ta_, data = ssf_data)
head(X)

stan_dat <- list(N = nrow(ssf_data), I = length(unique(ssf_data$step_id)), 
                 J = length(unique(ssf_data$id)), 
                 K = 5,
                 y = ssf_data$y, 
                 X=X[,2:6],
                 stepid = as.numeric(factor(ssf_data$step_id)), 
                 indid = ssf_data$id,
                 year=ssf_data$year,
                 season=as.numeric(as.factor(ssf_data$season)))

mod1<-stan(model_code = model,data=stan_dat,iter=500,chains=2,cores = 2)

file <- file.path("seasonmod.stan")
mod <- cmdstan_model(file)

fit_mcmc <- mod$sample(
  data = stan_dat,
  seed = 123,
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 40,
  iter_sampling = 20,
  init = 1
)



traceplot(mod1,pars=c("sigma"))
traceplot(mod1,pars=c("beta"))
print(mod1,pars=c("sigma"))
print(mod1,pars=c("beta"))
pairs(mod1,pars=c("beta"))
pairs(mod1,pars=c("sigma"))

print(mod1)


#Relative selection strength for Temperature across seasons, depths and distances to shore
# 1 = cold homo; 2 = hetero, 3 = hot homo
beta_temp<-rstan::extract(mod1,pars="beta_temp")
beta_depth<-rstan::extract(mod1,pars="beta_depth")
beta_dist<-rstan::extract(mod1,pars="beta_dist")
beta_td_int<-rstan::extract(mod1,pars="beta_td_int")
beta_tdist_int<-rstan::extract(mod1,pars="beta_tdist_int")
beta_ddist_int<-rstan::extract(mod1,pars="beta_ddist_int")
beta_3_int<-rstan::extract(mod1,pars="beta_threeway_int")


reftemp<-(17-mean(ssf_data$watertemp_))/sd(ssf_data$watertemp_)
temp<-seq(11.9,23,length=25)
st.temp<-(temp-mean(ssf_data$watertemp_))/sd(ssf_data$watertemp_)
depth<-seq(0,10,by=2)
st.depth<-(depth-mean(ssf_data$depth_m))/sd(ssf_data$depth_m)
dist<-c(10,210,410,610,810,1010)
st.dist<-(dist-mean(ssf_data$dist_to_shore_))/sd(ssf_data$dist_to_shore_)

dim(beta_temp$beta_temp)
dim(beta_temp$beta_temp)

rss.pred<-array(dim=c(8000,25,3,6,6))
for(i in 1:8000){
  for(t in 1:25){
    for(s in 1:3){
      for(d in 1:6){
        for(dist in 1:6){
      rss.pred[i,t,s,d,dist]<-beta_temp$beta_temp[i,s]*st.temp[t]+beta_depth$beta_depth[i,s]*st.depth[d]+beta_dist$beta_dist[i,s]*st.dist[dist]+
        beta_td_int$beta_td_int[i,s]*st.temp[t]*st.depth[d]+beta_ddist_int$beta_ddist_int[i,s]*st.temp[t]*st.dist[dist]+beta_tdist_int$beta_tdist_int[i,s]*st.temp[t]*st.dist[dist]+
        beta_3_int$beta_threeway_int[i,s]*st.temp[t]*st.dist[dist]*st.depth[d]
    }}}}}

rss.ref<-array(dim=c(8000,3,6,6))
for(i in 1:8000){
  for(s in 1:3){
    for(d in 1:6){
      for(dist in 1:6){
    rss.ref[i,s,d,dist]<-beta_temp$beta_temp[i,s]*reftemp+beta_depth$beta_depth[i,s]*st.depth[d]+
      beta_dist$beta_dist[i,s]*st.dist[dist]+
      beta_td_int$beta_td_int[i,s]*reftemp*st.depth[d]+
      beta_ddist_int$beta_ddist_int[i,s]*st.depth[d]*st.dist[dist]+
      beta_tdist_int$beta_tdist_int[i,s]*reftemp*st.dist[dist]+
      beta_3_int$beta_threeway_int[i,s]*reftemp*st.depth[d]*st.dist[dist]
  }}}}

rss<-array(dim=c(8000,25,3,6,6))
for(i in 1:8000){
  for(t in 1:25){
    for(s in 1:3){
      for(d in 1:6){
        for(dist in 1:6){
      rss[i,t,s,d,dist]<-exp(rss.pred[i,t,s,d,dist]-rss.ref[i,s,d,dist])  
    }}
  }}}

## Plots fill horizontally A-G, then H-M. its increase in depth on the x-axis then distance to shore on y

plot_list<-list()


depthindex<-rep(1:6,6)
distindex<-rep(1:6,each=6)

for(i in 1:36){
  var_list[[i]]<-c(depthindex[i],distindex[i])
}


for(i in 1:36){
 p = ggplot()+
      geom_line(aes_string(y=colMeans(rss[,,1,var_list[[i]][1],var_list[[i]][2]]),x=st.temp),col="blue")+
     # geom_line(aes_string(y=colMeans(rss[,,2,var_list[[i]][1],var_list[[i]][2]]),x=st.temp),col="purple")+
  #    geom_line(aes_string(y=colMeans(rss[,,3,var_list[[i]][1],var_list[[i]][2]]),x=st.temp),col="red")+
      geom_ribbon(aes_string(x=st.temp,ymin=colQuantiles(rss[,,1,var_list[[i]][1],var_list[[i]][2]],prob=0.025),ymax=colQuantiles(rss[,,1,var_list[[i]][1],var_list[[i]][2]],prob=0.975)),alpha=0.25,fill="blue")+
      #geom_ribbon(aes_string(x=st.temp,ymin=colQuantiles(rss[,,2,var_list[[i]][1],var_list[[i]][2]],prob=0.025),ymax=colQuantiles(rss[,,2,var_list[[i]][1],var_list[[i]][2]],prob=0.975)),alpha=0.25,fill="purple")+
    #    geom_ribbon(aes_string(x=st.temp,ymin=colQuantiles(rss[,,3,var_list[[i]][1],var_list[[i]][2]],prob=0.025),ymax=colQuantiles(rss[,,3,var_list[[i]][1],var_list[[i]][2]],prob=0.975)),alpha=0.25,fill="red")+
      theme_classic()+
   geom_hline(aes(yintercept=1),lty=2)

 
 plot_list[[i]] = p
}



(plot_list[[1]]+ggtitle("0m Depth")+ylab("10m Offshore

RSS")+
    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())
  |plot_list[[2]]+ggtitle("2m Depth")+
    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)
  |plot_list[[3]]+ggtitle("4m Depth")+
    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)
  |plot_list[[4]]+ggtitle("6m Depth")+
    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)
  |plot_list[[5]]+ggtitle("8m Depth")+
    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)
  |plot_list[[6]]+ggtitle("10m Depth")+
    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL))/
  (plot_list[[7]]+ylab("210m Offshore

RSS")+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())|plot_list[[8]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[9]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[10]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[11]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[12]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL))/
  (plot_list[[13]]+ylab("410m Offshore

RSS")+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())|plot_list[[14]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[15]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[16]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[17]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[18]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL))/
  (plot_list[[19]]+ylab("610m Offshore

RSS")+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())|plot_list[[20]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[21]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[22]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[23]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[24]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL))/
  (plot_list[[25]]+ylab("810m Offshore

RSS")+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())|plot_list[[26]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[27]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[28]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[29]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[30]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL))/
  (plot_list[[31]]+ylab("1010m Offshore

RSS")+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))|plot_list[[32]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL)
   |plot_list[[33]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL)|plot_list[[34]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL)|plot_list[[35]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL)|plot_list[[36]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL))










A<-ggplot()+
  geom_line(aes(y=colMeans(rss[,,1,1,1]),x=st.temp),col="blue")+
  geom_line(aes(y=colMeans(rss[,,2,1,1]),x=st.temp),col="purple")+
  geom_line(aes(y=colMeans(rss[,,3,1,1]),x=st.temp),col="red")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,1,1,1],prob=0.025),ymax=colQuantiles(rss[,,1,1,1],prob=0.975)),alpha=0.25,fill="blue")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,2,1,1],prob=0.025),ymax=colQuantiles(rss[,,2,1,1],prob=0.975)),alpha=0.25,fill="purple")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,3,1,1],prob=0.025),ymax=colQuantiles(rss[,,3,1,1],prob=0.975)),alpha=0.25,fill="red")+
  theme_classic()+
  ggtitle("0m Depth, 10m Offshore")
A

B<-ggplot()+
  geom_line(aes(y=colMeans(rss[,,1,2,1]),x=st.temp),col="blue")+
  geom_line(aes(y=colMeans(rss[,,2,2,1]),x=st.temp),col="purple")+
  geom_line(aes(y=colMeans(rss[,,3,2,1]),x=st.temp),col="red")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,1,2,1],prob=0.025),ymax=colQuantiles(rss[,,1,2,1],prob=0.975)),alpha=0.25,fill="blue")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,2,2,1],prob=0.025),ymax=colQuantiles(rss[,,2,2,1],prob=0.975)),alpha=0.25,fill="purple")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,3,2,1],prob=0.025),ymax=colQuantiles(rss[,,3,2,1],prob=0.975)),alpha=0.25,fill="red")+
  theme_classic()+
  ggtitle("2m Depth, 10m Offshore")
B

#4m  
C<-ggplot()+
  geom_line(aes(y=colMeans(rss[,,1,3]),x=st.temp),col="blue")+
  geom_line(aes(y=colMeans(rss[,,2,3]),x=st.temp),col="purple")+
  geom_line(aes(y=colMeans(rss[,,3,3]),x=st.temp),col="red")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,1,3],prob=0.025),ymax=colQuantiles(rss[,,1,3],prob=0.975)),alpha=0.25,fill="blue")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,2,3],prob=0.025),ymax=colQuantiles(rss[,,2,3],prob=0.975)),alpha=0.25,fill="purple")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,3,3],prob=0.025),ymax=colQuantiles(rss[,,3,3],prob=0.975)),alpha=0.25,fill="red")+
  theme_classic()+
  ggtitle("4m Depth")
C

#6m  
D<-ggplot()+
  geom_line(aes(y=colMeans(rss[,,1,4]),x=st.temp),col="blue")+
  geom_line(aes(y=colMeans(rss[,,2,4]),x=st.temp),col="purple")+
  geom_line(aes(y=colMeans(rss[,,3,4]),x=st.temp),col="red")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,1,4],prob=0.025),ymax=colQuantiles(rss[,,1,4],prob=0.975)),alpha=0.25,fill="blue")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,2,4],prob=0.025),ymax=colQuantiles(rss[,,2,4],prob=0.975)),alpha=0.25,fill="purple")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,3,4],prob=0.025),ymax=colQuantiles(rss[,,3,4],prob=0.975)),alpha=0.25,fill="red")+
  theme_classic()+
  ggtitle("6m Depth")
D

#8m  
E<-ggplot()+
  geom_line(aes(y=colMeans(rss[,,1,5]),x=st.temp),col="blue")+
  geom_line(aes(y=colMeans(rss[,,2,5]),x=st.temp),col="purple")+
  geom_line(aes(y=colMeans(rss[,,3,5]),x=st.temp),col="red")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,1,5],prob=0.025),ymax=colQuantiles(rss[,,1,5],prob=0.975)),alpha=0.25,fill="blue")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,2,5],prob=0.025),ymax=colQuantiles(rss[,,2,5],prob=0.975)),alpha=0.25,fill="purple")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,3,5],prob=0.025),ymax=colQuantiles(rss[,,3,5],prob=0.975)),alpha=0.25,fill="red")+
  theme_classic()+
  ggtitle("8m Depth")
E

#10m
G<-ggplot()+
  geom_line(aes(y=colMeans(rss[,,1,6]),x=st.temp),col="blue")+
  geom_line(aes(y=colMeans(rss[,,2,6]),x=st.temp),col="purple")+
  geom_line(aes(y=colMeans(rss[,,3,6]),x=st.temp),col="red")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,1,6],prob=0.025),ymax=colQuantiles(rss[,,1,6],prob=0.975)),alpha=0.25,fill="blue")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,2,6],prob=0.025),ymax=colQuantiles(rss[,,2,6],prob=0.975)),alpha=0.25,fill="purple")+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss[,,3,6],prob=0.025),ymax=colQuantiles(rss[,,3,6],prob=0.975)),alpha=0.25,fill="red")+
  theme_classic()+
  ggtitle("10m Depth")
G
(A|B|C)/(D|E|G)



#Relative selection strength for Depth across seasons
# 1 = cold homo; 2 = hetero, 3 = hot homo
beta<-rstan::extract(mod1,pars="beta_depth")
refdepth<-(4-mean(ssf_data$depth_m))/sd(ssf_data$depth_m)
depth<-seq(0,10,length=25)
st.depth<-(depth-mean(ssf_data$depth_m))/sd(ssf_data$depth_m)

dim(beta$beta_depth)

rss.pred<-array(dim=c(100,25,3))
for(i in 1:100){
  for(t in 1:25){
    for(s in 1:3){
      rss.pred[i,t,s]<-beta$beta_depth[i,s]*st.depth[t]
    }}}

rss.ref<-matrix(ncol=3,nrow=100)
for(i in 1:100){
  for(s in 1:3){
    rss.ref[i,s]<-beta$beta_depth[i,s]*refdepth
  }}

rss<-array(dim=c(100,25,3))
for(i in 1:100){
  for(t in 1:25){
    for(s in 1:3){
      rss[i,t,s]<-exp(rss.pred[i,t,s]-rss.ref[i,s])  
    }
  }}

ggplot()+
  geom_line(aes(y=colMeans(rss[,,1]),x=st.depth),col="blue")+
  geom_line(aes(y=colMeans(rss[,,2]),x=st.depth),col="purple")+
  geom_line(aes(y=colMeans(rss[,,3]),x=st.depth),col="red")+
  geom_ribbon(aes(x=st.depth,ymin=colQuantiles(rss[,,1],prob=0.025),ymax=colQuantiles(rss[,,1],prob=0.975)),alpha=0.25,fill="blue")+
  geom_ribbon(aes(x=st.depth,ymin=colQuantiles(rss[,,2],prob=0.025),ymax=colQuantiles(rss[,,2],prob=0.975)),alpha=0.25,fill="purple")+
  geom_ribbon(aes(x=st.depth,ymin=colQuantiles(rss[,,3],prob=0.025),ymax=colQuantiles(rss[,,3],prob=0.975)),alpha=0.25,fill="red")
