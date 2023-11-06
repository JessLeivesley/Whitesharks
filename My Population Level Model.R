## Population level model for Juvenile White Sharks

library(rstan)
library(tidyverse)
library(amt)
library(here)
library(tictoc)
library(bayesplot)
library(rgdal)
library(terra)
library(raster)
library(sf)
library(tmap, quietly = T)
library(geosphere)

# Read in the dataset
temp<-readRDS("available_temp_2020-2021_clean.RDS")
temp$DateTimePT<-as.factor(temp$DateTimePT)
temp$DateTime_Level<-as.numeric(as.factor(temp$DateTimePT))

# create dataset of datetime reference levels
dateindex<-temp%>%dplyr::select(DateTimeUTC,DateTime_Level)%>%distinct(DateTime_Level,.keep_all = T)
dateindex$DateTimeUTC<-as_datetime(dateindex$DateTimeUTC)


## Create date/time index
temp$DateTimePT<-as.factor(temp$DateTimePT)
temp$DateTime_Level<-as.numeric(as.factor(temp$DateTimePT))
## Set depth as a factor
temp$depth_m<-as.factor(temp$depth_m)
#
## how many levels of date/time are there and how many depth levels are there
DateTime<-levels(temp$DateTimePT)
nDT<-nlevels(temp$DateTimePT)
depthlevels<-levels(temp$depth_m)
#
#
### For each date/time (10095) create a spatial object of temperature available to the shark at #each 2m depth from the surface


for(dt in 1:nDT){
  raster.stack<-stack() # create empty raster
  temp1<-temp%>%filter(DateTimePT==DateTime[dt]) # filter date/time
  nd<-nlevels(droplevels(temp1$depth_m)) # number of depth levels taken
  depthlevels<-levels(droplevels(temp1$depth_m)) # Keep only depth levels taken
  for(d in 1:6){# for each depth
    temp2<-temp1%>%filter(depth_m==depthlevels[d])
    x <- raster(xmn=-119.6, xmx=-119.5, ymn=34.39, ymx=34.41, res=0.003) # Create the raster
    raster<- rasterize(temp2[, c('Center_Longitude_250m', 'Center_Latitude_250m')], x, temp2[, 'temp_c'], fun=mean)
    raster.stack<-stack(raster.stack,raster)
  }
  raster.stack<-projectRaster(raster.stack,crs = CRS("+init=epsg:5070")) # set projection
  names(raster.stack) <- c(depthlevels)
  assign(paste(dt),raster.stack) # assign the date/time index name
}


# join this to the shark data
vps<-readRDS("vps_2020-2021_full.RDS")
vps$DateTimeUTC<-as_datetime(vps$DateTimeUTC)
# vps<-left_join(vps,dateindex,by="DateTimeUTC" )
vps$year<-as.numeric(as.character(vps$year))


vps%>%filter(is.na(depth_m)==F)%>%group_by(shark)%>%summarise(minyr=min(year),maxyr=max(year),diff=maxyr-minyr)%>%filter(diff==1)


# put shark dataset into format we need - ID, long, lat timestamp 
alldata<-vps%>%filter(is.na(depth_m)==F)%>%dplyr::select(x = "lon",y="lat", t = "DateTimeUTC", id="shark",depth="depth_m")
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
alldata_list$sex<-vps%>%filter(is.na(temp_c)==F)%>%group_by(shark)%>%filter(row_number()==1)%>%dplyr::select(sex)

#make tracks and transform coordinates
alldata_list<-alldata_list%>%mutate(trk=lapply(data, function(d){
  amt::make_track(d,x,y,t, crs=sp::CRS("+proj=longlat +datum=WGS84"))%>%
    amt::transform_coords(sp::CRS("+init=epsg:5070"))
}))


print(alldata_list%>%mutate(sr=lapply(trk,summarize_sampling_rate))%>% dplyr::select(id,sr)%>%unnest(cols=c(sr)),n=22)
# Shark 2021-20 only has 9 points. 

alldata_list<-alldata_list[1:21,]

# resample so that moves are evenly spaced. 
dat1<-alldata_list%>%mutate(dat_clean=map(trk, ~ {
  .x %>% track_resample(rate = minutes(20), tolerance = seconds(120))
}))


# Explore the number of movements I have per individual
print(dat1%>% dplyr::select(id,data)%>%unnest(cols=c(data))%>%group_by(id)%>%summarise(mindepth=min(depth),maxdepth=max(depth)),n=21)


# create 10 random steps per individual (this will need to be increased, but for now 10 keeps time to run down)
set.seed(123)
dat_ssf <- dat1%>% 
  mutate(stps = map(dat_clean, ~ .x %>% filter_min_n_burst(min_n=3) %>%steps_by_burst() %>% random_steps(n_control = 500) ))
 
## Save this sampled dataset
saveRDS(dat_ssf,"500randomstps.RDS")
dat_ssf<-readRDS("100randomstps.RDS")


# create depth dataframe
depth.d<-vps%>%filter(is.na(temp_c)==F)%>%dplyr::select(DateTimeUTC,depth_m,shark)

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
for(l in 1:21){
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
cov.full <- vector(mode='list', length=21)
for(l in 1:21){
  for(DI in 1:10095){
    a<-dat_ssf$stps[[l]]%>%filter(DateTime_Level==DI)
    covs <- extract_covariates(a, get(paste(DI)))
    covs$id<-rep(dat_ssf$id[l],nrow(covs))
    cov.full[[l]]<-rbind(cov.full[[l]],covs)
  }}


# Only keep the watertemperature in the depth range that the shark is in / simulated to be in
for(l in 1:21){
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
    step_id = paste(id, step_id_, sep = "-"))

## Remove temp na's
ssf_data<-ssf_data%>%drop_na(watertemp_)

# select the first 10 random steps per movement (11 cases per stratum)
ssf_data%>%group_by(step_id)%>%summarise(meandepth=mean(depth_m))
b<-ssf_data%>%group_by(step_id)%>%filter(row_number()<12)
b%>%group_by(step_id)%>%count()%>%filter(n<11)

ssf_data<-ssf_data%>%group_by(step_id)%>%filter(row_number()<12)%>%ungroup()
#remove the two steps that cant be sampled!?
ssf_data<-ssf_data%>%filter(step_id!="2-740")
ssf_data<-ssf_data%>%filter(step_id!="4-23")

## Save this sampled dataset
#saveRDS(ssf_data,"10randomstps_withtemp.RDS")
ssf_data<-readRDS("10randomstps_withtemp.RDS")

## Transform x,y coordinates into lat/long
points<-cbind(ssf_data$x2_,ssf_data$y2_)
v <- vect(points, crs="+init=epsg:5070")
v

y <- project(v, "+proj=longlat +ellps=WGS84 +datum=WGS84")
y

lonlat <- geom(y)[, c("x", "y")]
head(lonlat, 3)
long<-geom(y)[,"x"]
lat<-geom(y)[,"y"]

ssf_data<-cbind(ssf_data,long=long,lat=lat)

## Distance to Shore 
ca = read_sf("Shoreline/Shoreline.shp")
# Cast the spatial lines into a multi-linestring and sample
samp = st_sample(st_cast(ca$geometry[1], "MULTILINESTRING"), 1000)
tm_shape(samp) + 
  tm_dots()
print(samp)

samp <- st_transform(samp, "+proj=longlat +ellps=WGS84 +datum=WGS84")

# Get California border coordinates from multipoint sf
samp_sf = st_as_sf(samp)
plot(samp_sf)
ca_border_coord = as.data.frame(st_coordinates(samp_sf))

# Function for calculating shark location to closest border point
shore_dist = function(s_lat, s_lon) {
  dists = tibble(id = seq(1, 1000, 1), 
                 dist = numeric(1000))
  for (i in 1:1000) {
    lat_dist = ca_border_coord$lat[i] - s_lat
    lon_dist = ca_border_coord$lon[i] - s_lon
    # dists$dist[i] = sqrt(lat_dist^2 + lon_dist^2)
    dists$dist[i] = distm(c(ca_border_coord$X[i], ca_border_coord$Y[i]), c(s_lon, s_lat), fun = distHaversine)
  }
  #sort_dists = dists[order(dists$dist), ]
  #closest_shore_point_i = sort_dists$id[1]
  #closest_dist = sort_dists$dist[1]
  return(min(dists$dist))
}

shark_shore_dist = function(data) {
  # Use the shark dataset we generated from before
  # put shark dataset into format we need - ID, long, lat timestamp and select one individual.
  # Calculate the shore distance for the movements (true and potential) for specific shark
  s_lat_list = data$lat
  s_lon_list = data$long
  shore_distance = numeric(nrow(data))
  for (i in 1:nrow(data)) {
    shore_distance[i] = shore_dist(s_lat_list[i], s_lon_list[i])
    # print(i) # Just to see how much has it run
  }
  
  return(shore_distance)
}

ssf_data$dist_to_shore_<-shark_shore_dist(ssf_data)

## Save this sampled dataset
# saveRDS(ssf_data,"10randomstps_withtempdist.RDS")
ssf_data<-readRDS("10randomstps_withtempdist.RDS")
# 
# ## Create water temperature seasonality
# ssf_data$Date<-as.Date(ssf_data$t2_h)
# temp<-readRDS("available_temp_2020-2021_clean.rds")
# names(temp)
# 
# sumtemp_surface<-temp%>%filter(depth_m ==0)%>%group_by(depth_m,Date,year)%>%summarise(meantemp_0 = mean(temp_c))
# 
# sumtemp_bottom<-temp%>%filter(depth_m ==10)%>%group_by(depth_m,Date,year)%>%summarise(meantemp_8 = mean(temp_c))
# 
# sumtemp<-cbind(sumtemp_surface,sumtemp_bottom)
# 
# sumtemp$diff<-sumtemp$meantemp_0-sumtemp$meantemp_8
# 
# hist(sumtemp$diff)
# 
# plot(diff~Date...2,sumtemp)
# library(segmented)
# 
# 
# 
# sumtemp$Date<-as.Date(sumtemp$Date...2)
# 
# low<-filter(sumtemp,diff<3.5)
# high<-filter(sumtemp,diff>3.5)
# 
# plot(diff~Date,low,ylim=c(0,8))
# points(diff~Date,high,col="red")
# # identify(high$Date,high$diff,n=4)
# # high[136,2]
# # which(sumtemp$Date...2=="2021-09-17")
# 
# sumtemp$Season<-rep(0,423)
# sumtemp$Season[1:26]<-"Homog"
# sumtemp$Season[27:98]<-"Heterog"
# sumtemp$Season[99:239]<-"Homog"
# sumtemp$Season[240:340]<-"Heterog"
# sumtemp$Season[341:423]<-"Homog"
# 
# ssf_data<-left_join(ssf_data,sumtemp, by="Date")

## Keep only 4 individuals to build and troubleshoot
# print(ssf_data%>%group_by(shark)%>%summarise(n=n()))
# 
# ssf_data<-ssf_data%>%mutate(st.WT=base::scale(watertemp_),st.WT2=base::scale(watertemp_^2),st.depth=base::scale(depth_m),st.depth2=base::scale(depth_m^2),st.distshore=base::scale(dist_to_shore_))
# unique(ssf_data$id)
# subset_ssf<-ssf_data%>%filter(id == 1 | id == 2 | id==3 | id == 4) # all 2020 sharks
# 
# ## Check of the data I have for all sharks
# print(ssf_data%>%group_by(shark)%>%summarise(n=n()),n=22)

ssf_data<-ssf_data%>%ungroup()%>%mutate(st.WT=base::scale(watertemp_),st.WT2=base::scale(watertemp_^2),st.depth=base::scale(depth_m),st.depth2=base::scale(depth_m^2),st.distshore=base::scale(dist_to_shore_),year=year(t2_))
unique(ssf_data$id)


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
}

parameters {
  vector[K] beta;
  vector<lower=0>[K] sigma;
  vector[I] a_re; // RE for steps
  vector[J] i_re1; // RE effects for each individual (temp)
  vector[J] i_re2; // RE effects for each individual (temp^2)
  vector[J] i_re3; // RE effects for each individual (depth)
  vector[J] i_re4; // RE effects for each individual (depth^2)
  vector[J] i_re5; // RE effects for each individual (distance to shore)
}


model {
  vector[N] mu;
  
  // priors
  sigma ~ normal(0, 5);
  a_re ~ normal(0, 100); // The paper has this as a_re ~ normal(0, 1000000)
  i_re1 ~ normal(0, sqrt(sigma[1]));
  i_re2 ~ normal(0, sqrt(sigma[2]));
  i_re3 ~ normal(0, sqrt(sigma[3]));
  i_re4 ~ normal(0, sqrt(sigma[4]));
  i_re5 ~ normal(0, sqrt(sigma[5]));
  beta ~ normal(0,5);

  //likelihood
 for(i in 1:N){
  
  mu[i] = a_re[stepid[i]]+
    (beta[1] + i_re1[indid[i]]) * X[i,1] + //temp effect 
    (beta[2] + i_re2[indid[i]]) * X[i,2]+
    (beta[3] + i_re3[indid[i]]) * X[i,3]+
    (beta[4] + i_re4[indid[i]]) * X[i,4]+
    (beta[5] + i_re5[indid[i]]) * X[i,5];
    
  //+ (beta[1] + i_re4[indid[i]]) * X[i,5] + //dummy month June
   //(beta[2] + i_re5[indid[i]]) * X[i,6] + //dummy month July
//   (beta[3] + i_re6[indid[i]]) * X[i,7] + //dummy month Aug
//   (beta[4] + i_re7[indid[i]]) * X[i,8] + //dummy month Sept
//   (beta[5] + i_re8[indid[i]]) * X[i,9] + //dummy month Oct
//   (beta[6] + i_re9[indid[i]]) * X[i,10] + //dummy month Nov
//   (beta[7] + i_re10[indid[i]]) * X[i,11]; //dummy month Dec
    
    
 }
  y ~ poisson_log(mu);
}
"

X <- model.matrix(y ~ st.WT + st.depth + st.WT2 + st.depth2+st.distshore, data = ssf_data)
X

stan_dat <- list(N = nrow(ssf_data), I = length(unique(ssf_data$step_id)), 
                 J = length(unique(ssf_data$id)), 
                 K = 5,
                 y = ssf_data$y, 
                 X=X[,2:6],
                 stepid = as.numeric(factor(ssf_data$step_id)), 
                 indid = ssf_data$id,
                 year=ssf_data$year)




mod1<-stan(model_code = model,data=stan_dat,iter=2500,chains=4,cores=4)
saveRDS(mod1,"fullmodel_4Nov.RDS")
mod1<-readRDS("fullmodel_4Nov.RDS")

traceplot(mod1,pars=c("sigma"))
traceplot(mod1,pars=c("beta"))
print(mod1,pars=c("sigma"))
print(mod1,pars=c("beta"))
pairs(mod1,pars=c("beta"))
pairs(mod1,pars=c("sigma"))

print(mod1)




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



## Plotting RSS
beta_temp<-rstan::extract(mod1,pars="beta[1]")
beta_temp2<-rstan::extract(mod1,pars="beta[2]")

#Relative selection strength for Temperature
range(vps$temp_c,na.rm=T) # 11.9 - 24.4
reftemp<-(12-mean(ssf_data$watertemp_))/sd(ssf_data$watertemp_)
reftemp2<-(12^2-mean(ssf_data$watertemp_^2))/sd(ssf_data$watertemp_^2)
temp<-seq(11.9,23,length=25)
st.temp<-(temp-mean(ssf_data$watertemp_))/sd(ssf_data$watertemp_)
st.temp2<-((temp^2)-mean((ssf_data$watertemp_^2)))/sd((ssf_data$watertemp_^2))

rss.pred<-matrix(ncol=25,nrow=2000)
for(i in 1:2000){
for(t in 1:25){
  rss.pred[i,t]<-beta_temp$`beta[1]`[i]*st.temp[t]+beta_temp2$`beta[2]`[i]*st.temp2[t]
}}

rss.ref<-vector()
for(i in 1:2000){
    rss.ref[i]<-beta_temp$`beta[1]`[i]*reftemp+beta_temp2$`beta[2]`[i]*reftemp2
}

rss<-matrix(ncol=25,nrow=2000)
for(i in 1:2000){
  for(t in 1:25){
    rss[i,t]<-exp(rss.pred[i,t]-rss.ref[i])  
  }
}

ggplot()+
  geom_line(aes(y=colMeans(rss),x=st.temp))+
  geom_ribbon(aes(x=st.temp,ymin=colQuantiles(rss,prob=0.025),ymax=colQuantiles(rss,prob=0.975)),alpha=0.5)+
  geom_hline(aes(yintercept=1))+
  theme_bw()+
  xlab("Temperature")+
  ylab("Relative Selection Strength
       (Reference: 12degrees)")+
  scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+
  theme(text=element_text(size=16))

#Relative selection strength for depth
range(vps$depth_m,na.rm=T) # 0- 10
refdepth<-(10-mean(ssf_data$depth_m))/sd(ssf_data$depth_m)
refdepth2<-((10^2)-mean((ssf_data$depth_m^2)))/sd((ssf_data$depth_m^2))
depth<-seq(0,10,length=25)
st.depth<-(depth-mean(ssf_data$depth_m))/sd(ssf_data$depth_m)
st.depth2<-((depth^2)-mean((ssf_data$depth_m^2)))/sd((ssf_data$depth_m^2))

beta_depth<-rstan::extract(mod1,pars=c("beta[3]","beta[4]"))

rss.pred<-matrix(ncol=25,nrow=2000)
for(i in 1:2000){
  for(t in 1:25){
    rss.pred[i,t]<-beta_depth$`beta[3]`[i]*st.depth[t]+beta_depth$`beta[4]`[i]*st.depth2[t]
  }}

rss.ref<-vector()
for(i in 1:2000){
  rss.ref[i]<-beta_depth$`beta[3]`[i]*refdepth+beta_depth$`beta[4]`[i]*refdepth2
}

rss<-matrix(ncol=25,nrow=2000)
for(i in 1:2000){
  for(t in 1:25){
    rss[i,t]<-exp(rss.pred[i,t]-rss.ref[i])  
  }
}

ggplot()+
  geom_line(aes(y=colMeans(rss),x=st.depth))+
  geom_ribbon(aes(x=st.depth,ymin=colQuantiles(rss,prob=0.025),ymax=colQuantiles(rss,prob=0.975)),alpha=0.5)+
  geom_hline(aes(yintercept=1))+
  theme_bw()+
  xlab("Depth (m)")+
  ylab("Relative Selection Strength
       (Reference: 10m)")+
  scale_x_continuous(breaks=c(-1,0,1,2),labels=c(1.45,4.15,6.85,9.55))+
  theme(text=element_text(size=16))

#Relative selection strength for distance to shore
range(ssf_data$dist_to_shore_,na.rm=T) # 0.09 - 1450
refdist<-(100-mean(ssf_data$dist_to_shore_))/sd(ssf_data$dist_to_shore_)
dist<-seq(0,1450,length=25)
st.dist<-(dist-mean(ssf_data$dist_to_shore_))/sd(ssf_data$dist_to_shore_)


beta_dist<-rstan::extract(mod1,pars=c("beta[5]"))

rss.pred<-matrix(ncol=25,nrow=2000)
for(i in 1:2000){
  for(t in 1:25){
    rss.pred[i,t]<-beta_dist$`beta[5]`[i]*st.dist[t]
  }}

rss.ref<-vector()
for(i in 1:2000){
  rss.ref[i]<-beta_dist$`beta[5]`[i]*refdist
}

rss<-matrix(ncol=25,nrow=2000)
for(i in 1:2000){
  for(t in 1:25){
    rss[i,t]<-exp(rss.pred[i,t]-rss.ref[i])  
  }
}

ggplot()+
  geom_line(aes(y=colMeans(rss),x=st.dist))+
  geom_ribbon(aes(x=st.dist,ymin=colQuantiles(rss,prob=0.025),ymax=colQuantiles(rss,prob=0.975)),alpha=0.5)+
  geom_hline(aes(yintercept=1))+
  theme_bw()+
  xlab("Distance to Shore (m)")+
  ylab("Relative Selection Strength
       (Reference: 100m)")+
  scale_x_continuous(breaks=c(-1,0,1,2,3,4),labels=c(148,416,683,950,1218,1485))+
  theme(text=element_text(size=16))


# Plot individual lines
#Relative selection strength for Temperature
beta_temp<-rstan::extract(mod1,pars="beta[1]")
beta_temp2<-rstan::extract(mod1,pars="beta[2]")
ind<-rstan::extract(mod1,pars=c("i_re1","i_re2"))

range(vps$temp_c,na.rm=T) # 11.9 - 24.4
reftemp<-(12-mean(ssf_data$watertemp_))/sd(ssf_data$watertemp_)
reftemp2<-(12^2-mean(ssf_data$watertemp_^2))/sd(ssf_data$watertemp_^2)
temp<-seq(11.9,23,length=25)
st.temp<-(temp-mean(ssf_data$watertemp_))/sd(ssf_data$watertemp_)
st.temp2<-((temp^2)-mean((ssf_data$watertemp_^2)))/sd((ssf_data$watertemp_^2))

temp_short<-seq(11.9,19,length=25)
st.temp<-(temp-mean(ssf_data$watertemp_))/sd(ssf_data$watertemp_)
st.temp2<-((temp^2)-mean((ssf_data$watertemp_^2)))/sd((ssf_data$watertemp_^2))

rss.pred<-array(dim = c(5000,25,21))
for(i in 1:5000){
  for(t in 1:25){
    for(s in 1:21){
    rss.pred[i,t,s]<-(beta_temp$`beta[1]`[i]+ind$i_re1[i,s])*st.temp[t]+(beta_temp2$`beta[2]`[i]+ind$i_re2[i,s])*st.temp2[t]
  }}}

rss.ref<-matrix(ncol=21,nrow=5000)
for(i in 1:5000){
  for(s in 1:21){
  rss.ref[i,s]<-(beta_temp$`beta[1]`[i]+ind$i_re1[i,s])*reftemp+(beta_temp2$`beta[2]`[i]+ind$i_re2[i,s])*reftemp2
}}

rss<-array(dim = c(5000,25,21))
for(i in 1:5000){
  for(t in 1:25){
    for(s in 1:21){
    rss[i,t,s]<-exp(rss.pred[i,t,s]-rss.ref[i,s])  
  }
}}

ggplot()+
  geom_line(aes(y=colMeans(rss[,,1]),x=st.temp))+
  geom_line(aes(y=colMeans(rss[,,2]),x=st.temp))+
  geom_line(aes(y=colMeans(rss[,,3]),x=st.temp))+ 
  geom_line(aes(y=colMeans(rss[,,4]),x=st.temp))+ 
  geom_line(aes(y=colMeans(rss[,,5]),x=st.temp))+ 
  geom_line(aes(y=colMeans(rss[,,6]),x=st.temp))+ 
  geom_line(aes(y=colMeans(rss[,,7]),x=st.temp))+ 
  geom_line(aes(y=colMeans(rss[,,8]),x=st.temp))+ 
  geom_line(aes(y=colMeans(rss[,,9]),x=st.temp))+ 
  geom_line(aes(y=colMeans(rss[,,10]),x=st.temp))+
  geom_line(aes(y=colMeans(rss[,,11]),x=st.temp))+ 
 # geom_line(aes(y=colMeans(rss[,1:16,12]),x=st.temp[1:16]))+ ## 12 is odd
 # geom_line(aes(y=colMeans(rss[,1:19,13]),x=st.temp[1:19]))+ ## 13 is odd
  geom_line(aes(y=colMeans(rss[,,14]),x=st.temp))+ 
  geom_line(aes(y=colMeans(rss[,,15]),x=st.temp))+
  geom_line(aes(y=colMeans(rss[,,16]),x=st.temp))+
  #geom_line(aes(y=colMeans(rss[,,17]),x=st.temp)) ## 17 is odd
  geom_line(aes(y=colMeans(rss[,,18]),x=st.temp))+
  geom_line(aes(y=colMeans(rss[,,19]),x=st.temp))+
  geom_line(aes(y=colMeans(rss[,,20]),x=st.temp))+
  geom_line(aes(y=colMeans(rss[,,21]),x=st.temp))+
  geom_hline(aes(yintercept=1))+
  theme_bw()+
  xlab("Temperature")+
  ylab("Relative Selection Strength
       (Reference: 12degrees)")+
  scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+
  theme(text=element_text(size=16))

ggplot()+
  geom_line(aes(y=colMeans(rss[,1:16,12]),x=st.temp[1:16]))+
  geom_line(aes(y=colMeans(rss[,1:16,13]),x=st.temp[1:16]))+
  geom_line(aes(y=colMeans(rss[,1:16,17]),x=st.temp[1:16]))+  
  geom_hline(aes(yintercept=1))+
  theme_bw()+
  xlab("Temperature")+
  ylab("Relative Selection Strength
       (Reference: 12degrees)")+
  scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+
  theme(text=element_text(size=16))


unique(ssf_data$shark)[12]
ssf_data%>%filter(shark=="2021-23" & case_==TRUE)%>%ungroup()%>%summarise(mintemp=min(watertemp_),maxtemp=max(watertemp_))
# only 9 steps - temps only experienced up to 19.7 degrees

unique(ssf_data$shark)[13]
ssf_data%>%filter(shark=="2021-31" & case_==TRUE)%>%ungroup()%>%summarise(mintemp=min(watertemp_),maxtemp=max(watertemp_))
# only 35 steps - temps only experienced up to 20.3 degrees

unique(ssf_data$shark)[17]
ssf_data%>%filter(shark=="2021-50" & case_==TRUE)%>%ungroup()%>%summarise(mintemp=min(watertemp_),maxtemp=max(watertemp_))
# only 20 steps - temps only experienced up to 17.6 degrees


#Relative selection strength for depth
range(vps$depth_m,na.rm=T) # 0- 10
refdepth<-(10-mean(ssf_data$depth_m))/sd(ssf_data$depth_m)
refdepth2<-((10^2)-mean((ssf_data$depth_m^2)))/sd((ssf_data$depth_m^2))
depth<-seq(0,10,length=25)
st.depth<-(depth-mean(ssf_data$depth_m))/sd(ssf_data$depth_m)
st.depth2<-((depth^2)-mean((ssf_data$depth_m^2)))/sd((ssf_data$depth_m^2))

beta_depth<-rstan::extract(mod1,pars=c("beta[3]","beta[4]"))
ind<-rstan::extract(mod1,pars=c("i_re3","i_re4"))

rss.pred<-array(dim = c(5000,25,21))
for(i in 1:5000){
  for(t in 1:25){
    for(s in 1:21){
    rss.pred[i,t,s]<-(beta_depth$`beta[3]`[i]+ind$i_re3[i,s])*st.depth[t]+(beta_depth$`beta[4]`[i]+ind$i_re4[i,s])*st.depth2[t]
  }}}

rss.ref<-matrix(ncol=21,nrow=5000)
for(i in 1:5000){
  for(s in 1:21){
  rss.ref[i,s]<-(beta_depth$`beta[3]`[i]+ind$i_re3[i,s])*refdepth+(beta_depth$`beta[4]`[i]+ind$i_re4[i,s])*refdepth2
}}

rss<-array(dim = c(5000,25,21))
for(i in 1:5000){
  for(t in 1:25){
    for(s in 1:21){
    rss[i,t,s]<-exp(rss.pred[i,t,s]-rss.ref[i,s])  
  }
}}

ggplot()+
  # geom_line(aes(y=colMeans(rss[,,1]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,2]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,3]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,4]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,5]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,6]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,7]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,8]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,9]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,10]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,11]),x=st.depth))+
   geom_line(aes(y=colMeans(rss[,,12]),x=st.depth))+ ## 12 is odd
  geom_line(aes(y=colMeans(rss[,,13]),x=st.depth))+ ## 13 is odd
  # geom_line(aes(y=colMeans(rss[,,14]),x=st.depth))+ 
  # geom_line(aes(y=colMeans(rss[,,15]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,16]),x=st.depth))+
  geom_line(aes(y=colMeans(rss[,,17]),x=st.depth))+ ## 17 is odd
  geom_line(aes(y=colMeans(rss[,,18]),x=st.depth))+
  # geom_line(aes(y=colMeans(rss[,,19]),x=st.depth))+
  geom_line(aes(y=colMeans(rss[,,20]),x=st.depth))+
  geom_line(aes(y=colMeans(rss[,,21]),x=st.depth))+
  geom_hline(aes(yintercept=1))+
  theme_bw()+
  xlab("Depth (m)")+
  ylab("Relative Selection Strength
       (Reference: 10m)")+
  scale_x_continuous(breaks=c(-1,0,1,2),labels=c(1.45,4.15,6.85,9.55))+
  theme(text=element_text(size=16))

#Relative selection strength for distance to shore
range(ssf_data$dist_to_shore_,na.rm=T) # 0.09 - 1450
refdist<-(100-mean(ssf_data$dist_to_shore_))/sd(ssf_data$dist_to_shore_)
dist<-seq(0,1450,length=25)
st.dist<-(dist-mean(ssf_data$dist_to_shore_))/sd(ssf_data$dist_to_shore_)


beta_dist<-rstan::extract(mod1,pars=c("beta[5]"))
ind<-rstan::extract(mod1,pars=c("i_re5"))

rss.pred<-array(dim = c(5000,25,21))
for(i in 1:5000){
  for(t in 1:25){
    for(s in 1:21){
    rss.pred[i,t,s]<-(beta_dist$`beta[5]`[i]+ind$i_re5[i,s])*st.dist[t]
  }}}

rss.ref<-matrix(ncol=21,nrow=5000)
for(i in 1:5000){
  for(s in 1:21){
  rss.ref[i,s]<-(beta_dist$`beta[5]`[i]+ind$i_re5[i,s])*refdist
}}

rss<-array(dim = c(5000,25,21))
for(i in 1:5000){
  for(t in 1:25){
    for(s in 1:21){
    rss[i,t,s]<-exp(rss.pred[i,t,s]-rss.ref[i,s])  
  }
}}

ggplot()+
  geom_line(aes(y=colMeans(rss[,,1]),x=st.dist))+
  geom_line(aes(y=colMeans(rss[,,2]),x=st.dist))+
  geom_line(aes(y=colMeans(rss[,,3]),x=st.dist))+ 
  geom_line(aes(y=colMeans(rss[,,4]),x=st.dist))+ 
  geom_line(aes(y=colMeans(rss[,,5]),x=st.dist))+ 
  geom_line(aes(y=colMeans(rss[,,6]),x=st.dist))+ 
  geom_line(aes(y=colMeans(rss[,,7]),x=st.dist))+ 
  geom_line(aes(y=colMeans(rss[,,8]),x=st.dist))+ 
  geom_line(aes(y=colMeans(rss[,,9]),x=st.dist))+ 
  geom_line(aes(y=colMeans(rss[,,10]),x=st.dist))+
  geom_line(aes(y=colMeans(rss[,,11]),x=st.dist))+ 
  geom_line(aes(y=colMeans(rss[,,12]),x=st.dist))+ ## 12 is odd
  geom_line(aes(y=colMeans(rss[,,13]),x=st.dist))+ ## 13 is odd
  geom_line(aes(y=colMeans(rss[,,14]),x=st.dist))+ 
  geom_line(aes(y=colMeans(rss[,,15]),x=st.dist))+
  geom_line(aes(y=colMeans(rss[,,16]),x=st.dist))+
  geom_line(aes(y=colMeans(rss[,,17]),x=st.dist))+ ## 17 is odd
  geom_line(aes(y=colMeans(rss[,,18]),x=st.dist))+
  geom_line(aes(y=colMeans(rss[,,19]),x=st.dist))+
  geom_line(aes(y=colMeans(rss[,,20]),x=st.dist))+
  geom_line(aes(y=colMeans(rss[,,21]),x=st.dist))+
  geom_hline(aes(yintercept=1))+
  theme_bw()+
  xlab("Distance to Shore (m)")+
  ylab("Relative Selection Strength
       (Reference: 100m)")+
  scale_x_continuous(breaks=c(-1,0,1,2,3,4),labels=c(148,416,683,950,1218,1485))+
  theme(text=element_text(size=16))






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
