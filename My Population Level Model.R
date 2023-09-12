## Population level model for Juvenile White Sharks

library(cmdstanr)
library(tidyverse)
library(amt)
library(here)
library(tictoc)
library(bayesplot)
library(rgdal)

# Read in the dataset
# create dataset of datetime reference levels
dateindex<-temp%>%select(DateTimeUTC,DateTime_Level)%>%distinct(DateTime_Level,.keep_all = T)
dateindex$DateTimeUTC<-as_datetime(dateindex$DateTimeUTC)
# join this to the shark data
vps<-readRDS("vps_2020-2021_full.RDS")
vps$DateTimeUTC<-as_datetime(vps$DateTimeUTC)
vps_datein<-left_join(vps,dateindex,by="DateTimeUTC" )


# put shark dataset into format we need - ID, long, lat timestamp 
alldata<-vps%>%filter(is.na(temp_c)==F)%>%select(x = "lon",y="lat", t = "DateTimeUTC", id="shark",depth=depth_m)
alldata$t<-as_datetime(alldata$t)

# there are 10 duplicate moves that need to be removed
alldata%>%group_by(id,t)%>%summarise(n=n())%>%filter(n > 1)
# remove the second of each
alldata<-alldata[-which(alldata$id=="2020-20" & alldata$t==as_datetime("2020-07-21 09:33:35"))[2],]
alldata<-alldata[-which(alldata$id=="2020-20" & alldata$t==as_datetime("2020-11-10 10:34:57"))[2],]
alldata<-alldata[-which(alldata$id=="2020-21" & alldata$t==as_datetime("2020-06-19 04:34:36"))[2],]
alldata<-alldata[-which(alldata$id=="2020-21" & alldata$t==as_datetime("2020-06-19 11:08:09"))[2],]
alldata<-alldata[-which(alldata$id=="2020-21" & alldata$t==as_datetime("2020-08-18 08:15:29"))[2],]
alldata<-alldata[-which(alldata$id=="2020-21" & alldata$t==as_datetime("2020-10-13 03:34:36"))[2],]
alldata<-alldata[-which(alldata$id=="2020-33" & alldata$t==as_datetime("2020-09-11 04:40:54"))[2],]
alldata<-alldata[-which(alldata$id=="2020-32" & alldata$t==as_datetime("2020-07-17 05:48:09"))[2],]
alldata<-alldata[-which(alldata$id=="2020-32" & alldata$t==as_datetime("2020-09-21 16:47:53"))[2],]
alldata<-alldata[-which(alldata$id=="2020-32" & alldata$t==as_datetime("2020-10-03 08:29:50"))[2],]





## Create a dataframe of nested dataframes
alldata_list<-alldata%>%group_by(id)%>%nest()

# add sex column
alldata_list$sex<-vps%>%filter(is.na(temp_c)==F)%>%group_by(shark)%>%filter(row_number()==1)%>%select(sex)

#make tracks and transform coordinates
alldata_list<-alldata_list%>%mutate(trk=lapply(data, function(d){
  amt::make_track(d,x,y,t, crs=sp::CRS("+proj=longlat +datum=WGS84"))%>%
    amt::transform_coords(sp::CRS("+init=epsg:5070"))
}))


print(alldata_list%>%mutate(sr=lapply(trk,summarize_sampling_rate))%>% select(id,sr)%>%unnest(cols=c(sr)),n=22)

# resample so that moves are evenly spaced. 
dat1<-alldata_list%>%mutate(dat_clean=map(trk, ~ {
  .x %>% track_resample(rate = minutes(10), tolerance = seconds(120))
}))

# create 10 random steps per individual (this will need to be increased, but for now 10 keeps time to run down)
dat_ssf <- dat1 %>% 
  mutate(stps = map(dat_clean, ~ .x %>% filter_min_n_burst(min_n=3) %>%steps_by_burst() %>% 
                      random_steps() ))
 

# create depth dataframe
depth.d<-vps%>%filter(is.na(temp_c)==F)%>%select(DateTimeUTC,depth_m,shark)

# remove the second of each
depth.d<-depth.d[-which(depth.d$shark=="2020-20" & depth.d$DateTimeUTC==as_datetime("2020-07-21 09:33:35"))[2],]
depth.d<-depth.d[-which(depth.d$shark=="2020-20" & depth.d$DateTimeUTC==as_datetime("2020-11-10 10:34:57"))[2],]
depth.d<-depth.d[-which(depth.d$shark=="2020-21" & depth.d$DateTimeUTC==as_datetime("2020-06-19 04:34:36"))[2],]
depth.d<-depth.d[-which(depth.d$shark=="2020-21" & depth.d$DateTimeUTC==as_datetime("2020-06-19 11:08:09"))[2],]
depth.d<-depth.d[-which(depth.d$shark=="2020-21" & depth.d$DateTimeUTC==as_datetime("2020-08-18 08:15:29"))[2],]
depth.d<-depth.d[-which(depth.d$shark=="2020-21" & depth.d$DateTimeUTC==as_datetime("2020-10-13 03:34:36"))[2],]
depth.d<-depth.d[-which(depth.d$shark=="2020-33" & depth.d$DateTimeUTC==as_datetime("2020-09-11 04:40:54"))[2],]
depth.d<-depth.d[-which(depth.d$shark=="2020-32" & depth.d$DateTimeUTC==as_datetime("2020-07-17 05:48:09"))[2],]
depth.d<-depth.d[-which(depth.d$shark=="2020-32" & depth.d$DateTimeUTC==as_datetime("2020-09-21 16:47:53"))[2],]
depth.d<-depth.d[-which(depth.d$shark=="2020-32" & depth.d$DateTimeUTC==as_datetime("2020-10-03 08:29:50"))[2],]


# For each shark get the depth it was recorded at
for(l in 1:22){
  dat_ssf$stps[[l]]$t2_h<-round(dat_ssf$stps[[l]]$t2_,units="hours")
  dat_ssf$stps[[l]]<-left_join(dat_ssf$stps[[l]],dateindex,by=c("t2_h"="DateTimeUTC"))
  shark.ob<- dat_ssf$id[l]
  dep.fil<-depth.d%>%filter(shark == shark.ob)
  dat_ssf$stps[[l]]<-left_join(dat_ssf$stps[[l]],dep.fil, by=c("t2_"="DateTimeUTC"))
  for(i in 1:length(dat_ssf$stps[[l]]$x1_)){
    dat_ssf$stps[[l]]$depth_m[i]<-ifelse(dat_ssf$stps[[l]]$case_[i]=="FALSE",runif(1,0,10),dat_ssf$stps[[l]]$depth_m[i])
  }
}

# Now extract temperatures at each of the six depth levels for all of the movements
cov.full <- vector(mode='list', length=22)
for(l in 1:22){
  for(DI in 1:10095){
    a<-dat_ssf$stps[[l]]%>%filter(DateTime_Level==DI)
    covs <- extract_covariates(a, get(paste(DI)))
    covs$id<-rep(dat_ssf$id[l],nrow(covs))
    cov.full[[l]]<-rbind(cov.full[[l]],covs)
  }}


# Only keep the watertemperature in the depth range that the shark is in / simulated to be in
for(l in 1:22){
  cov.full[[l]]$watertemp_<-ifelse(cov.full[[l]]$depth_m<1,cov.full[[l]]$X0,
                                   ifelse(cov.full[[l]]$depth_m>=1 & cov.full[[l]]$depth_m <3, cov.full[[l]]$X2, 
                                          ifelse(cov.full[[l]]$depth_m >= 3 & cov.full[[l]]$depth_m <5, cov.full[[l]]$X4,
                                                 ifelse(cov.full[[l]]$depth_m>=5 &cov.full[[l]]$depth_m<7, cov.full[[l]]$X6,
                                                        ifelse(cov.full[[l]]$depth_m>=7 &cov.full[[l]]$depth_m<9, cov.full[[l]]$X8,cov.full[[l]]$X10)))))
}

#unlist the object
ssf_data<-bind_rows(cov.full)

# add id column
ssf_data<-ssf_data%>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    step_id = paste0(id, step_id_, sep = "-"))

## Remove temp na's
ssf_data<-ssf_data%>%drop_na(watertemp_)

# Distance to shore calculations from Kaylen's code - but distances arent on the correct scale. 
# ## Find distance to shore for each move and proposed move
# california <- readOGR("ca-state-boundary/CA_State_TIGER2016.shp")
# # change the CRS for california data to have the same CRS with shark data
# california2 <- spTransform(california, CRS("+proj=longlat +datum=WGS84"))
# # find the shortest distance between the location of the shark and the coastline of California
# 
# ssf_data$xtran<-terra::project(as.matrix(ssf_data[,c(3,5)]), from="+init=epsg:5070", to="+proj=longlat +datum=WGS84")[,1]
# ssf_data$ytran<-terra::project(as.matrix(ssf_data[,c(3,5)]), from="+init=epsg:5070", to="+proj=longlat +datum=WGS84")[,2]
# 
# distfromshore<-dist2Line(ssf_data[,c(28,29)],california2,distfun=distGeo)
# hist(distfromshore)
# ssf_data<-cbind(ssf_data,distfromshore=distfromshore[,1])
# 
# # drop distance NAs
# ssf_data<-ssf_data%>%drop_na(distfromshore)


## Create water temperature seasonality
ssf_data$Date<-as.Date(ssf_data$t2_h)
temp<-readRDS("available_temp_2020-2021_clean.rds")
names(temp)

sumtemp_surface<-temp%>%filter(depth_m ==0)%>%group_by(depth_m,Date,year)%>%summarise(meantemp_0 = mean(temp_c))

sumtemp_bottom<-temp%>%filter(depth_m ==10)%>%group_by(depth_m,Date,year)%>%summarise(meantemp_8 = mean(temp_c))

sumtemp<-cbind(sumtemp_surface,sumtemp_bottom)

sumtemp$diff<-sumtemp$meantemp_0-sumtemp$meantemp_8

hist(sumtemp$diff)

plot(diff~Date...2,sumtemp)
library(segmented)



sumtemp$Date<-as.Date(sumtemp$Date...2)

low<-filter(sumtemp,diff<3.5)
high<-filter(sumtemp,diff>3.5)

plot(diff~Date,low,ylim=c(0,8))
points(diff~Date,high,col="red")
# identify(high$Date,high$diff,n=4)
# high[136,2]
# which(sumtemp$Date...2=="2021-09-17")

sumtemp$Season<-rep(0,423)
sumtemp$Season[1:26]<-"Homog"
sumtemp$Season[27:98]<-"Heterog"
sumtemp$Season[99:239]<-"Homog"
sumtemp$Season[240:340]<-"Heterog"
sumtemp$Season[341:423]<-"Homog"

ssf_data<-left_join(ssf_data,sumtemp, by="Date")

## Keep only 4 individuals to build and troubleshoot
print(ssf_data%>%group_by(shark)%>%summarise(n=n(),year=mean(year...7)),n=22)

ssf_data<-ssf_data%>%mutate(st.WT=base::scale(watertemp_),st.WT2=base::scale(watertemp_^2),st.depth=base::scale(depth_m),st.depth2=base::scale(depth_m^2))
unique(ssf_data$id)
subset_ssf<-ssf_data%>%filter(id == 1 | id == 2 | id==3 | id == 4) # all 2020 sharks

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
  int<lower=1, upper=8> month[N];
}

parameters {
  vector[8] beta_temp; // fixed effects
  vector[8] beta_depth;
  vector[8] beta_temp2; // fixed effects
  vector[8] beta_depth2;
  vector[7] beta;
  vector[I] a_re; // RE for steps
  vector[J] i_re1; // RE effects for each individual (temp)
  vector[J] i_re2; // RE effects for each individual (temp^2)
  vector[J] i_re3; // RE effects for each individual (depth)
  vector[J] i_re4; // RE effects for each individual (depth^2)
  vector[J] i_re5; // RE effects for each individual (distance to shore)
  vector[J] i_re6; // RE effects for each individual (distance to shore^2)
  vector[J] i_re7;
  vector[J] i_re8;
  vector[J] i_re9;
    vector[J] i_re10;
     vector[J] i_re11;
  real<lower = 0> sigmaind1;
  real<lower = 0> sigmaind2;
  real<lower = 0> sigmaind3;
  real<lower = 0> sigmaind4;
  real<lower = 0> sigmaind5;
  real<lower = 0> sigmaind6;
  real<lower = 0> sigmaind7;
  real<lower = 0> sigmaind8;
  real<lower = 0> sigmaind9;
    real<lower = 0> sigmaind10;
    real<lower = 0> sigmaind11;
  real<lower = 0> sigmastps;
}


model {
  vector[N] mu;
  
  // priors
  sigmaind1 ~ normal(0, 3);
  sigmaind2 ~ normal(0, 3);
  sigmaind3 ~ normal(0, 3);
  sigmaind4 ~ normal(0, 3);
  sigmaind5 ~ normal(0, 3);
  sigmaind6 ~ normal(0, 3);
  sigmaind7 ~ normal(0, 3);
  sigmaind8 ~ normal(0, 3);
  sigmaind9 ~ normal(0, 3);
    sigmaind10 ~ normal(0, 3);
     sigmaind11 ~ normal(0, 3);
  sigmastps ~ normal(0,100);
  a_re ~ normal(0, sigmastps); // The paper has this as a_re ~ normal(0, 1000000)
  i_re1 ~ normal(0, sqrt(sigmaind1));
  i_re2 ~ normal(0, sqrt(sigmaind2));
  i_re3 ~ normal(0, sqrt(sigmaind3));
  i_re4 ~ normal(0, sqrt(sigmaind4));
  i_re5 ~ normal(0, sqrt(sigmaind5));
  i_re6 ~ normal(0, sqrt(sigmaind6));
  i_re7 ~ normal(0, sqrt(sigmaind7));
  i_re8 ~ normal(0, sqrt(sigmaind8));
  i_re9 ~ normal(0, sqrt(sigmaind9));
    i_re10 ~ normal(0, sqrt(sigmaind10));
    i_re11 ~ normal(0, sqrt(sigmaind11));
  beta_temp ~ normal(0,5);
  beta_depth ~ normal(0,5);
  beta ~normal(0,5);

  //likelihood
 for(i in 1:N){
  
  mu[i] = a_re[stepid[i]] + 
    (beta_temp[month[i]] + i_re1[indid[i]]) * X[i,1] + //temp effect with slope for each month
    (beta_depth[month[i]] + i_re2[indid[i]]) * X[i,2]+
    (beta_temp2[month[i]] + i_re3[indid[i]]) * X[i,3]+
    (beta_depth2[month[i]] + i_re11[indid[i]]) * X[i,4]+
   (beta[1] + i_re4[indid[i]]) * X[i,5] + //dummy month June
   (beta[2] + i_re5[indid[i]]) * X[i,6] + //dummy month July
   (beta[3] + i_re6[indid[i]]) * X[i,7] + //dummy month Aug
   (beta[4] + i_re7[indid[i]]) * X[i,8] + //dummy month Sept
   (beta[5] + i_re8[indid[i]]) * X[i,9] + //dummy month Oct
   (beta[6] + i_re9[indid[i]]) * X[i,10] + //dummy month Nov
   (beta[7] + i_re10[indid[i]]) * X[i,11]; //dummy month Dec
    
    
 }
  y ~ poisson_log(mu);
}
"

X <- model.matrix(y ~ st.WT + st.depth + st.WT2 + st.depth2+factor(Month), data = subset_ssf)
X
# compile the model
tic()
f <- write_stan_file(model)
stan_mod <- cmdstan_model(f)
toc()

subset_ssf$Season[subset_ssf$Season=="Homog"]<-0
subset_ssf$Season[subset_ssf$Season=="Heterog"]<-1

subset_ssf$year...7[subset_ssf$year...7==2020]<-0
subset_ssf$year...7[subset_ssf$year...7==2021]<-1

subset_ssf$Month<-(month(subset_ssf$Date))
summary(subset_ssf$Month)

subset_ssf$Month[subset_ssf$Month==5]<-1
subset_ssf$Month[subset_ssf$Month==6]<-2
subset_ssf$Month[subset_ssf$Month==7]<-3
subset_ssf$Month[subset_ssf$Month==8]<-4
subset_ssf$Month[subset_ssf$Month==9]<-5
subset_ssf$Month[subset_ssf$Month==10]<-6
subset_ssf$Month[subset_ssf$Month==11]<-7
subset_ssf$Month[subset_ssf$Month==12]<-8

summary(as.numeric(subset_ssf$Season))
hist(subset_ssf$distfromshore)
hist(ssf_data$distfromshore)
ssf_data%>%group_by(year...3)%>%summarise(meanT=mean(watertemp_),minT=min(watertemp_),maxT=max(watertemp_))

subset_ssf$watertempst<-(subset_ssf$watertemp_-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
subset_ssf%>%group_by(Season)%>%summarise(meanT=mean(watertempst),minT=min(watertempst),maxT=max(watertempst))

# Standardise the covariates

stan_dat <- list(N = nrow(subset_ssf), I = length(unique(subset_ssf$step_id)), 
                 J = length(unique(subset_ssf$id)), 
                 M = length(unique(subset_ssf$Month)),
                 y = subset_ssf$y, 
                 temp = (subset_ssf$watertemp_-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp), 
                 temp2 = ((subset_ssf$watertemp_)^2-mean((subset_ssf$watertemp_)^2))/sd((subset_ssf$watertemp_)^2), 
                 depth = (subset_ssf$depth_m-mean(subset_ssf$depth_m))/sd(subset_ssf$depth_m),
                 depth2 = ((subset_ssf$depth_m^2)-mean((subset_ssf$depth_m^2)))/sd((subset_ssf$depth_m^2)),
                 distance = (subset_ssf$distfromshore-mean(subset_ssf$distfromshore))/sd(subset_ssf$distfromshore),
                 distance2 = ((subset_ssf$distfromshore^2)-mean((subset_ssf$distfromshore^2)))/sd((subset_ssf$distfromshore^2)),
                 month = as.numeric(subset_ssf$Month),
                 stepid = as.numeric(factor(subset_ssf$step_id)), 
                 indid = subset_ssf$id)

stan_dat <- list(N = nrow(subset_ssf), I = length(unique(subset_ssf$step_id)), 
                 J = length(unique(subset_ssf$id)), 
                 K = 11,
                 y = subset_ssf$y, 
                 X=X[,2:12],
                 month=subset_ssf$Month,
                 stepid = as.numeric(factor(subset_ssf$step_id)), 
                 indid = subset_ssf$id)

tic()
res_stan <- stan_mod$sample(data=stan_dat, chains = 2, iter_warmup=500,iter_sampling = 500,parallel_chains=2)
toc()

print(res_stan$summary(c("sigmaind1","sigmaind2","sigmaind3","sigmaind4","sigmaind8","sigmaind9","sigmaind7","beta_temp","beta_depth","beta_temp2","a_re[1]","a_re[2]","a_re[100]","a_re[200]")),n=30)

pars=c("sigmaind1","sigmaind2","beta_temp2[3]","beta_temp2[4]","beta_temp[5]","beta_temp[6]","beta_temp[7]")
bayesplot::mcmc_trace(res_stan$draws(), pars = pars)


## Plotting RSS
beta_temp<-res_stan$draws(variable="beta_temp",format="df")
beta_temp2<-res_stan$draws(variable="beta_temp2",format="df")
beta_depth<-res_stan$draws(variable="beta_depth",format="df")
beta_depth2<-res_stan$draws(variable="beta_depth2",format="df")
beta<-res_stan$draws(variable="beta",format="df")

#Temperature in the first month (May), only including temperature coefficients here as everything else is held at 0
reftemp<-(15.5-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
reftemp2<-(15.5^2-mean(subset_ssf$watertemp^2))/sd(subset_ssf$watertemp^2)
temp_may<-seq(13.9,22,length=15)
temp_mayst<-(temp_may-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
temp2_may<-(temp_may^2-mean(subset_ssf$watertemp^2))/sd(subset_ssf$watertemp^2)

pred_may<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
 pred_may[i,t]  <- beta_temp$`beta_temp[1]`[i]*temp_mayst[t]+beta_temp2$`beta_temp2[1]`[i]*temp2_may[t]
  }
}

pred_may_ref<-matrix(nrow=200,ncol=1)
  for(i in 1:200){
    pred_may_ref[i,]  <- beta_temp$`beta_temp[1]`[i]*reftemp+beta_temp2$`beta_temp2[1]`[i]*reftemp2
  }

rss_may<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
rss_may[i,t]<-exp(pred_may[i,t]-pred_may_ref[i,1])
  }}

temp_june<-seq(15.5,20.3,length=15)
temp_junest<-(temp_june-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
temp2_june<-(temp_june^2-mean(subset_ssf$watertemp^2))/sd(subset_ssf$watertemp^2)

pred_june<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_june[i,t]  <- beta_temp$`beta_temp[2]`[i]*temp_junest[t]+beta_temp2$`beta_temp2[2]`[i]*temp2_june[t]+beta$`beta[1]`[i]
  }
}

pred_june_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_june_ref[i,]  <- beta_temp$`beta_temp[2]`[i]*reftemp+beta_temp2$`beta_temp2[2]`[i]*reftemp2+beta$`beta[1]`[i]
}

rss_june<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_june[i,t]<-exp(pred_june[i,t]-pred_june_ref[i,1])
  }}

temp_july<-seq(15.4,20.6,length=15)
temp_julyst<-(temp_july-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
temp2_july<-(temp_july^2-mean(subset_ssf$watertemp^2))/sd(subset_ssf$watertemp^2)

pred_july<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_july[i,t]  <- beta_temp$`beta_temp[3]`[i]*temp_julyst[t]+beta_temp2$`beta_temp2[3]`[i]*temp2_july[t]+beta$`beta[2]`[i]
  }
}

pred_july_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_july_ref[i,]  <- beta_temp$`beta_temp[3]`[i]*reftemp+beta_temp2$`beta_temp2[3]`[i]*reftemp2+beta$`beta[2]`[i]}

rss_july<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_july[i,t]<-exp(pred_july[i,t]-pred_july_ref[i,1])
  }}

temp_aug<-seq(15.4,21.5,length=15)
temp_augst<-(temp_aug-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
temp2_aug<-(temp_aug^2-mean(subset_ssf$watertemp^2))/sd(subset_ssf$watertemp^2)

pred_aug<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_aug[i,t]  <- beta_temp$`beta_temp[4]`[i]*temp_augst[t]+beta_temp2$`beta_temp2[4]`[i]*temp2_aug[t]+beta$`beta[3]`[i]
  }
}

pred_aug_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_aug_ref[i,]  <- beta_temp$`beta_temp[4]`[i]*reftemp+beta_temp2$`beta_temp2[4]`[i]*reftemp2
  +beta$`beta[3]`[i]}

rss_aug<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_aug[i,t]<-exp(pred_aug[i,t]-pred_aug_ref[i,1])
  }}

temp_sept<-seq(15.5,20.4,length=15)
temp_septst<-(temp_sept-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
temp2_sept<-(temp_sept^2-mean(subset_ssf$watertemp^2))/sd(subset_ssf$watertemp^2)

pred_sept<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_sept[i,t]  <- beta_temp$`beta_temp[5]`[i]*temp_sept[t]+ beta_temp2$`beta_temp2[5]`[i]*temp2_sept[t]
  }
}

pred_sept_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_sept_ref[i,]  <- beta_temp$`beta_temp[5]`[i]*reftemp+ beta_temp2$`beta_temp2[5]`[i]*reftemp2
}

rss_sept<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_sept[i,t]<-exp(pred_sept[i,t]-pred_sept_ref[i,1])
  }}

temp_oct<-seq(15.2,21.7,length=15)
temp_octst<-(temp_oct-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
temp2_oct<-(temp_oct^2-mean(subset_ssf$watertemp^2))/sd(subset_ssf$watertemp^2)

pred_oct<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_oct[i,t]  <- beta_temp$`beta_temp[6]`[i]*temp_octst[t]+beta_temp2$`beta_temp2[6]`[i]*temp2_oct[t]
  }
}

pred_oct_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_oct_ref[i,]  <- beta_temp$`beta_temp[6]`[i]*reftemp+beta_temp2$`beta_temp2[6]`[i]*reftemp2
}

rss_oct<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_oct[i,t]<-exp(pred_oct[i,t]-pred_oct_ref[i,1])
  }}

temp_nov<-seq(12.9,17.7,length=15)
temp_novst<-(temp_nov-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
temp2_nov<-(temp_nov^2-mean(subset_ssf$watertemp^2))/sd(subset_ssf$watertemp^2)


pred_nov<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_nov[i,t]  <- beta_temp$`beta_temp[7]`[i]*temp_novst[t]+beta_temp2$`beta_temp2[7]`[i]*temp2_nov[t]
  }
}

pred_nov_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_nov_ref[i,]  <- beta_temp$`beta_temp[7]`[i]*reftemp+beta_temp2$`beta_temp2[7]`[i]*reftemp2
}

rss_nov<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_nov[i,t]<-exp(pred_nov[i,t]-pred_nov_ref[i,1])
  }}

temp_dec<-seq(14.2,15.9,length=15)
temp_decst<-(temp_dec-mean(subset_ssf$watertemp))/sd(subset_ssf$watertemp)
temp2_dec<-(temp_dec^2-mean(subset_ssf$watertemp^2))/sd(subset_ssf$watertemp^2)

pred_dec<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_dec[i,t]  <- beta_temp$`beta_temp[8]`[i]*temp_decst[t]+beta_temp2$`beta_temp2[8]`[i]*temp2_dec[t]
  }
}

pred_dec_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_dec_ref[i,]  <- beta_temp$`beta_temp[8]`[i]*reftemp+beta_temp2$`beta_temp2[8]`[i]*reftemp2
}

rss_dec<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_dec[i,t]<-exp(pred_dec[i,t]-pred_dec_ref[i,1])
  }}



ggplot()+
  # geom_line(aes(x=temp_may,y=colQuantiles(rss_may,prob=0.5)), linewidth=1,col="red")+
  # geom_ribbon(aes(x=temp_may,ymin=colQuantiles(rss_may,prob=0.025),ymax=colQuantiles(rss_may, prob=0.975)),alpha=0.25,fill="red")+
  scale_y_continuous(limits=c(0,50))+
  geom_line(aes(x=temp_june,y=colQuantiles(rss_june,prob=0.5)), linewidth=1,col="orange")+
  geom_ribbon(aes(x=temp_june,ymin=colQuantiles(rss_june,prob=0.025),ymax=colQuantiles(rss_june, prob=0.975)),alpha=0.25,fill="orange")+
  geom_line(aes(x=temp_july,y=colQuantiles(rss_july,prob=0.5)), linewidth=1,col="yellow")+
  geom_ribbon(aes(x=temp_july,ymin=colQuantiles(rss_july,prob=0.025),ymax=colQuantiles(rss_july, prob=0.975)),alpha=0.25,fill="yellow")+
  geom_line(aes(x=temp_aug,y=colQuantiles(rss_aug,prob=0.5)), linewidth=1,col="green")+
  geom_ribbon(aes(x=temp_aug,ymin=colQuantiles(rss_aug,prob=0.025),ymax=colQuantiles(rss_aug, prob=0.975)),alpha=0.25,fill="green")
  geom_line(aes(x=temp_sept,y=colQuantiles(rss_sept,prob=0.5)), linewidth=1,col="blue")+
  geom_ribbon(aes(x=temp_sept,ymin=colQuantiles(rss_sept,prob=0.025),ymax=colQuantiles(rss_sept, prob=0.975)),alpha=0.25,fill="blue")+
  geom_line(aes(x=temp_oct,y=colQuantiles(rss_oct,prob=0.5)), linewidth=1,col="purple")+
  geom_ribbon(aes(x=temp_oct,ymin=colQuantiles(rss_oct,prob=0.025),ymax=colQuantiles(rss_oct, prob=0.975)),alpha=0.25,fill="purple")+
  geom_line(aes(x=temp_nov,y=colQuantiles(rss_nov,prob=0.5)), linewidth=1,col="lightblue")+
  geom_ribbon(aes(x=temp_nov,ymin=colQuantiles(rss_nov,prob=0.025),ymax=colQuantiles(rss_nov, prob=0.975)),alpha=0.25,fill="lightblue")+
  geom_line(aes(x=temp_dec,y=colQuantiles(rss_dec,prob=0.5)), linewidth=1,col="pink")+
  geom_ribbon(aes(x=temp_dec,ymin=colQuantiles(rss_dec,prob=0.025),ymax=colQuantiles(rss_dec, prob=0.975)),alpha=0.25,fill="pink")+
  theme_classic()+
  xlab("Temperature")+
  ylab("RSS (16.3 degrees)")+
  geom_hline(aes(yintercept=1),col="red", linewidth=1,lty=2)+
  scale_y_continuous(limits=c(0,150))


#deptherature in the first month (May), only including deptherature coefficients here as everything else is held at 0
refdepth<-0
refdepth2<-0
depth_may<-seq(0,8.5,length=15)
depth_mayst<-(depth_may-mean(subset_ssf$depth_m))/sd(subset_ssf$depth_m)
depth2_may<-(depth_may^2-mean(subset_ssf$depth_m^2))/sd(subset_ssf$depth_m^2)

pred_may<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_may[i,t]  <- beta_depth$`beta_depth[1]`[i]*depth_mayst[t]+beta_depth2$`beta_depth2[1]`[i]*depth2_may[t]
  }
}

pred_may_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_may_ref[i,]  <- beta_depth$`beta_depth[1]`[i]*refdepth+beta_depth2$`beta_depth2[1]`[i]*refdepth2
}

rss_may<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_may[i,t]<-exp(pred_may[i,t]-pred_may_ref[i,1])
  }}

depth_june<-seq(0,9.1,length=15)
depth_junest<-(depth_june-mean(subset_ssf$depth_m))/sd(subset_ssf$depth_m)
depth2_june<-(depth_june^2-mean(subset_ssf$depth_m^2))/sd(subset_ssf$depth_m^2)

pred_june<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_june[i,t]  <- beta_depth$`beta_depth[2]`[i]*depth_junest[t]+beta_depth2$`beta_depth2[2]`[i]*depth2_june[t]
  }
}

pred_june_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_june_ref[i,]  <- beta_depth$`beta_depth[2]`[i]*refdepth+beta_depth2$`beta_depth2[2]`[i]*refdepth2
}

rss_june<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_june[i,t]<-exp(pred_june[i,t]-pred_june_ref[i,1])
  }}

depth_july<-seq(0,10,length=15)
depth_julyst<-(depth_july-mean(subset_ssf$depth_m))/sd(subset_ssf$depth_m)
depth2_july<-(depth_july^2-mean(subset_ssf$depth_m^2))/sd(subset_ssf$depth_m^2)

pred_july<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_july[i,t]  <- beta_depth$`beta_depth[3]`[i]*depth_julyst[t]+beta_depth2$`beta_depth2[3]`[i]*depth2_july[t]
  }
}

pred_july_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_july_ref[i,]  <- beta_depth$`beta_depth[3]`[i]*refdepth+beta_depth2$`beta_depth2[3]`[i]*refdepth2
}

rss_july<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_july[i,t]<-exp(pred_july[i,t]-pred_july_ref[i,1])
  }}

depth_aug<-seq(0,10,length=15)
depth_augst<-(depth_aug-mean(subset_ssf$depth_m))/sd(subset_ssf$depth_m)
depth2_aug<-(depth_aug^2-mean(subset_ssf$depth_m^2))/sd(subset_ssf$depth_m^2)

pred_aug<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_aug[i,t]  <- beta_depth$`beta_depth[4]`[i]*depth_augst[t]+beta_depth2$`beta_depth2[4]`[i]*depth2_aug[t]
  }
}

pred_aug_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_aug_ref[i,]  <- beta_depth$`beta_depth[4]`[i]*refdepth+beta_depth2$`beta_depth2[4]`[i]*refdepth2
}

rss_aug<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_aug[i,t]<-exp(pred_aug[i,t]-pred_aug_ref[i,1])
  }}

depth_sept<-seq(0,10,length=15)
depth_septst<-(depth_sept-mean(subset_ssf$depth_m))/sd(subset_ssf$depth_m)
depth2_sept<-(depth_sept^2-mean(subset_ssf$depth_m^2))/sd(subset_ssf$depth_m^2)

pred_sept<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_sept[i,t]  <- beta_depth$`beta_depth[5]`[i]*depth_sept[t]+ beta_depth2$`beta_depth2[5]`[i]*depth2_sept[t]
  }
}

pred_sept_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_sept_ref[i,]  <- beta_depth$`beta_depth[5]`[i]*refdepth+ beta_depth2$`beta_depth2[5]`[i]*refdepth2
}

rss_sept<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_sept[i,t]<-exp(pred_sept[i,t]-pred_sept_ref[i,1])
  }}

depth_oct<-seq(0,10,length=15)
depth_octst<-(depth_oct-mean(subset_ssf$depth_m))/sd(subset_ssf$depth_m)
depth2_oct<-(depth_oct^2-mean(subset_ssf$depth_m^2))/sd(subset_ssf$depth_m^2)

pred_oct<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_oct[i,t]  <- beta_depth$`beta_depth[6]`[i]*depth_octst[t]+beta_depth2$`beta_depth2[6]`[i]*depth2_oct[t]
  }
}

pred_oct_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_oct_ref[i,]  <- beta_depth$`beta_depth[6]`[i]*refdepth+beta_depth2$`beta_depth2[6]`[i]*refdepth2
}

rss_oct<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_oct[i,t]<-exp(pred_oct[i,t]-pred_oct_ref[i,1])
  }}

depth_nov<-seq(0,10,length=15)
depth_novst<-(depth_nov-mean(subset_ssf$depth_m))/sd(subset_ssf$depth_m)
depth2_nov<-(depth_nov^2-mean(subset_ssf$depth_m^2))/sd(subset_ssf$depth_m^2)


pred_nov<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_nov[i,t]  <- beta_depth$`beta_depth[7]`[i]*depth_novst[t]+beta_depth2$`beta_depth2[7]`[i]*depth2_nov[t]
  }
}

pred_nov_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_nov_ref[i,]  <- beta_depth$`beta_depth[7]`[i]*refdepth+beta_depth2$`beta_depth2[7]`[i]*refdepth2
}

rss_nov<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_nov[i,t]<-exp(pred_nov[i,t]-pred_nov_ref[i,1])
  }}

depth_dec<-seq(0,10,length=15)
depth_decst<-(depth_dec-mean(subset_ssf$depth_m))/sd(subset_ssf$depth_m)
depth2_dec<-(depth_dec^2-mean(subset_ssf$depth_m^2))/sd(subset_ssf$depth_m^2)

pred_dec<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    pred_dec[i,t]  <- beta_depth$`beta_depth[8]`[i]*depth_decst[t]+beta_depth2$`beta_depth2[8]`[i]*depth2_dec[t]
  }
}

pred_dec_ref<-matrix(nrow=200,ncol=1)
for(i in 1:200){
  pred_dec_ref[i,]  <- beta_depth$`beta_depth[8]`[i]*refdepth+beta_depth2$`beta_depth2[8]`[i]*refdepth2
}

rss_dec<-matrix(nrow=200,ncol=15)
for(t in 1:15){
  for(i in 1:200){
    rss_dec[i,t]<-exp(pred_dec[i,t]-pred_dec_ref[i,1])
  }}



ggplot()+
  # geom_line(aes(x=depth_may,y=colQuantiles(rss_may,prob=0.5)), linewidth=1,col="red")+
  # geom_ribbon(aes(x=depth_may,ymin=colQuantiles(rss_may,prob=0.025),ymax=colQuantiles(rss_may, prob=0.975)),alpha=0.25,fill="red")+
  # geom_line(aes(x=depth_june,y=colQuantiles(rss_june,prob=0.5)), linewidth=1,col="orange")+
  # geom_ribbon(aes(x=depth_june,ymin=colQuantiles(rss_june,prob=0.025),ymax=colQuantiles(rss_june, prob=0.975)),alpha=0.25,fill="orange")+
  geom_line(aes(x=depth_july,y=colQuantiles(rss_july,prob=0.5)), linewidth=1,col="yellow")+
  geom_ribbon(aes(x=depth_july,ymin=colQuantiles(rss_july,prob=0.025),ymax=colQuantiles(rss_july, prob=0.975)),alpha=0.25,fill="yellow")+
  # geom_line(aes(x=depth_aug,y=colQuantiles(rss_aug,prob=0.5)), linewidth=1,col="green")+
  # geom_ribbon(aes(x=depth_aug,ymin=colQuantiles(rss_aug,prob=0.025),ymax=colQuantiles(rss_aug, prob=0.975)),alpha=0.25,fill="green")+
  # geom_line(aes(x=depth_sept,y=colQuantiles(rss_sept,prob=0.5)), linewidth=1,col="blue")+
  # geom_ribbon(aes(x=depth_sept,ymin=colQuantiles(rss_sept,prob=0.025),ymax=colQuantiles(rss_sept, prob=0.975)),alpha=0.25,fill="blue")+
  # geom_line(aes(x=depth_oct,y=colQuantiles(rss_oct,prob=0.5)), linewidth=1,col="purple")+
  # # geom_ribbon(aes(x=depth_oct,ymin=colQuantiles(rss_oct,prob=0.025),ymax=colQuantiles(rss_oct, prob=0.975)),alpha=0.25,fill="purple")+
  # geom_line(aes(x=depth_nov,y=colQuantiles(rss_nov,prob=0.5)), linewidth=1,col="lightblue")+
  # geom_ribbon(aes(x=depth_nov,ymin=colQuantiles(rss_nov,prob=0.025),ymax=colQuantiles(rss_nov, prob=0.975)),alpha=0.25,fill="lightblue")+
  # geom_line(aes(x=depth_dec,y=colQuantiles(rss_dec,prob=0.5)), linewidth=1,col="pink")+
  # geom_ribbon(aes(x=depth_dec,ymin=colQuantiles(rss_dec,prob=0.025),ymax=colQuantiles(rss_dec, prob=0.975)),alpha=0.25,fill="pink")+
  theme_classic()+
  xlab("Depth (m)")+
  ylab("RSS (0 m)")+
  geom_hline(aes(yintercept=1),col="red", linewidth=1,lty=2)









#Temperature in the Hetero season
temp<-seq(-1,1,length=15)
season<-rep(1,15)

pred_he<-matrix(nrow=10,ncol=15)
for(t in 1:15){
  for(i in 1:10){
    pred_he[i,t]  <- betas$`beta[1]`[i]*temp[t]+betas$`beta[2]`[i]*(temp[t]^2)+betas$`beta[5]`[i]*season[t]+betas$`beta[6]`[i]*temp[t]*season[t]+betas$`beta[7]`[i]*(temp[t]^2)*season[t]
  }
}

pred_he_ref<-matrix(nrow=10,ncol=1)
for(i in 1:10){
  pred_he_ref[i,]  <- betas$`beta[1]`[i]*0+betas$`beta[2]`[i]*0+betas$`beta[5]`[i]*season[t]+betas$`beta[6]`[i]*0*season[t]+betas$`beta[7]`[i]*0*season[t]
}

rss_he<-matrix(nrow=10,ncol=15)
for(t in 1:15){
  for(i in 1:10){
    rss_he[i,t]<-exp(pred_he[i,t]-pred_he_ref[i,1])
  }}


ggplot()+
  geom_line(aes(x=temp,y=colQuantiles(rss_ho,prob=0.5)), linewidth=1,col="purple")+
  geom_ribbon(aes(x=temp,ymin=colQuantiles(rss_ho,prob=0.025),ymax=colQuantiles(rss_ho, prob=0.975)),alpha=0.25,fill="purple")+
  geom_line(aes(x=temp,y=colQuantiles(rss_he,prob=0.5)), linewidth=1,col="green")+
  geom_ribbon(aes(x=temp,ymin=colQuantiles(rss_he,prob=0.025),ymax=colQuantiles(rss_he, prob=0.975)),alpha=0.25,fill="green")+
  theme_classic()+
  xlab("Temperature")+
  ylab("RSS (16.3 degrees)")+
  geom_hline(aes(yintercept=1),col="red", linewidth=1,lty=2)+
  scale_y_continuous(limits = c(0,500))






## GLMMTMB to compare
library(glmmTMB)
subset_ssf$Season<-as.numeric(subset_ssf$Season)
TMBStruc <- glmmTMB(y ~ -1 + st.WT*factor(Month) +I(st.WT^2)*factor(Month)  +(1|step_id) + 
                      (0 + st.WT*factor(Month)  +I(st.WT^2)*factor(Month)  | id),
                    family=poisson, data = subset_ssf, doFit=FALSE,control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))) 

#' Set the value of the standard deviation of the first random effect (here (1|step_id)):
TMBStruc$parameters$theta[1] <- log(1e3) 

#' Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
TMBStruc$mapArg <- list(theta=factor(c(NA,rep(0,300))))

#' Fit the model and look at the summary:
glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
summary(glmm.TMB.random)

#' 95\% CIs for fixed and random effects (standard deviations) are obtained via the confint() function:
confint(glmm.TMB.random)



subset_ssf%>%filter(y==1)%>%group_by(Season)%>%summarise(meantemp=mean(watertemp_))
subset_ssf%>%filter(y==0)%>%group_by(Season)%>%summarise(meantemp=mean(watertemp_))
