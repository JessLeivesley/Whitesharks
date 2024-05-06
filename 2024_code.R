## NEW METHOD FOR SAMPLING TEMPERATURE AND DEPTH - MUCH QUICKER AND ALLOWS MORE RANDOM STEPS TO BE SAMPLED.

## ---- Libraries ----
library(ggplot2)
library(lubridate)
library(raster)
library(rgdal)
library(amt)
library(sf)
library(dplyr)
library(geosphere)
library(tidyr)
library(tmap)
library(nplyr)
library(rstan)
library(fitdistrplus)
library(glmmTMB)

## ---- Create step object for all individuals ----
# read in all shark data
vps <- readRDS("vps_2020-2021_full.RDS")

# chose one shark and put dataset into format we need - ID, long, lat, and timestamp
alldata <- vps %>% 
  filter(is.na(temp_c)==F & is.na(depth_m)==F) %>% 
  dplyr::select(x = "lon",y="lat", t = "DateTimeUTC", id="shark",depth=depth_m,temp="temp_c",depth_bin2m="depth_bin2m")

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

# now we are going to make a track from this data, while also adjusting the GPS coordinates (needs to be the same projection as the temperature data which we can play around with) 
alldata_list<-alldata_list%>%mutate(trk=lapply(data, function(d){
  amt::make_track(d,.x=x, .y=y, .t=t, depth_bin2m=depth_bin2m,depth=depth, crs=4326)
}))

# look at the sampling resolution
print(alldata_list%>%mutate(sr=lapply(trk,summarize_sampling_rate))%>% dplyr::select(id,sr)%>%unnest(cols=c(sr)),n=22)

# Resample at the selected resolution
dat1<-alldata_list%>%mutate(dat_clean=map(trk, ~ {
  .x %>% track_resample(rate = minutes(20), tolerance = seconds(120))
}))

# Explore the number of movements I have per individual
print(dat1%>% dplyr::select(id,data)%>%unnest(cols=c(data))%>%group_by(id)%>%summarise(mindepth=min(depth),maxdepth=max(depth)),n=22)

# Now we are going to sample the potential moves an individual could make. Sample ~300
dat_ssf <- dat1 %>%
  mutate(
    stps = map(dat_clean, ~ .x %>%
                 amt::filter_min_n_burst(min_n = 3) %>%
                 amt::steps_by_burst()))

# Shark 2021-20 only has 9 points. 
dat_ssf<-dat_ssf[1:21,]

dat_ssf <- dat_ssf %>%
  mutate(
    stps = map(stps, ~ mutate(.x, sl_ = replace(sl_, sl_ == 0, 1e-5)))
  )

dat_ssf2 <- dat_ssf %>%
  mutate(dist = map(stps, ~fitdist((.x$sl_), "gamma")))
# Warning messages are okay.

set.seed(123)
for(i in 1:21){
  print(i)
  
  dat_ssf2$rando[[i]]<-random_steps(dat_ssf2$stps[[i]],
                                    n_control = 50, ## Update this to 300 eventually
                                    sl_distr=fitdist(dat_ssf2$stps[[i]]$sl_,"gamma"),
                                    rand_sl = rgamma(n = 1e05,
                                                     shape = dat_ssf2$dist[[i]]$estimate[1],
                                                     rate = dat_ssf2$dist[[i]]$estimate[2]))
 
}


## ---- Add polygon ID to the shark dataset ----
# this allows us to add temperature at the selected resolution

# read in polygonID
polygon = read_sf("PolygonID/Grid_250m.shp")
polygon
polygon<-st_transform(polygon, "+proj=longlat +ellps=WGS84 +datum=WGS84")

# create a raster object and add ID
r <- raster()
extent(r) <- extent(polygon)
rp <- rasterize(polygon, r, 'ID')
plot(rp)

# extract the polygonID for each movement
for(i in 1:21){
  print(i)
  dat_ssf2$rando[[i]]<-extract_covariates(dat_ssf2$rando[[i]],rp)
  #remove moves outside of our array
  dat_ssf2$rando[[i]]<-dat_ssf2$rando[[i]]%>%dplyr::filter(is.na(dat_ssf2$rando[[i]]$layer)==F)
}



## ---- Now add depth to true and simulated movement ----

# Get the max depth per polygon
# this will be used to simulate the potential depth the individual could be at
temp<-readRDS("available_temp_2020-2021_clean.RDS")
Polygon_depth<-temp%>%group_by(PolygonID_250m)%>%summarise(MaxDepth=mean(Total_Water_Colum_m))

# Combine total water column depth with all moves
for(i in 1:21){
  print(i)
  dat_ssf2$rando[[i]]<-left_join(dat_ssf2$rando[[i]],Polygon_depth,by=c("layer"="PolygonID_250m"))
  dat_ssf2$rando[[i]]<-left_join(dat_ssf2$rando[[i]],dat_ssf2$trk[[i]],by=c("x2_"="x_","y2_"="y_","t2_"="t_"))
  #Polygon 54 doesnt exist, and Polygon 55 has negative water column depth??
  dat_ssf2$rando[[i]]<-dat_ssf2$rando[[i]]%>%filter(layer!=55)
  dat_ssf2$rando[[i]]<-dat_ssf2$rando[[i]]%>%filter(layer!=54)
}


set.seed(12)
for(i in 1:21){
  # Simulate depth
  dat_ssf2$rando[[i]]$depth<-ifelse(dat_ssf2$rando[[i]]$case_==TRUE,dat_ssf2$rando[[i]]$depth,dunif(1,0,dat_ssf2$rando[[i]]$MaxDepth))
  # Now create depth bin
  dat_ssf2$rando[[i]]$depthbin<-ifelse(dat_ssf2$rando[[i]]$depth<1,"0-1m",
                                       ifelse(dat_ssf2$rando[[i]]$depth<3,"2-3m",
                                              ifelse(dat_ssf2$rando[[i]]$depth<5,"4-5m",
                                                     ifelse(dat_ssf2$rando[[i]]$depth<7,"6-7m",
                                                            ifelse(dat_ssf2$rando[[i]]$depth<9,"8-9m","10-11m")))))
}



## ---- Now add water temperature to true and simulated movement ----

for(i in 1:21){
  print(i)
  # create rounded hour
  dat_ssf2$rando[[i]]$t2_h<-round(dat_ssf2$rando[[i]]$t2_,units="hours")
  # now combine temperature based on hour
  dat_ssf2$rando[[i]]<-left_join(dat_ssf2$rando[[i]],temp,by=c("layer"="PolygonID_250m","depthbin"="depth_bin2m","t2_h"="DateTimeUTC"))
}

# There are 404 missing temperature records - where the temperature array wasnt measure temp yet!
# remove them
# also remove missing turn angle measures
for(i in 1:21){
  dat_ssf2$rando[[i]]<-dat_ssf2$rando[[i]]%>%dplyr::filter(is.na(dat_ssf2$rando[[i]]$temp_c)==F)
  dat_ssf2$rando[[i]]<-dat_ssf2$rando[[i]]%>%dplyr::filter(is.na(dat_ssf2$rando[[i]]$ta_)==F)
  dat_ssf2$rando[[i]]$sharkid<-rep(dat_ssf2$id[[i]],nrow(dat_ssf2$rando[[i]]))
}


# ---- Unlist the data and format ----

# Now unlist the dataset
ssf_data<-bind_rows(dat_ssf2$rando)
as.data.frame(ssf_data)%>%group_by(sharkid)%>%count()
glimpse(ssf_data)

# format shark id column, and create a unique id for each step across sharks 
ssf_data<-ssf_data%>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(sharkid)), 
    step_id = paste(id, step_id_, sep = "-"))

# look at the number of rows per individual
as.data.frame(ssf_data)%>%group_by(id)%>%count()

# how many random moves per move am I left with per step? 
as.data.frame(ssf_data)%>%filter(case_==FALSE)%>%group_by(id,step_id)%>%count()

# Keep the 20 random moves per step 
## Update this once its confirmed it works
ssf_data<-as.data.frame(ssf_data)%>%group_by(step_id)%>%dplyr::filter(row_number()<21)%>%ungroup()

# ---- Add in distance to shore ----
`# Read in shoreline shapefile
ca = read_sf("Shoreline/Shoreline.shp")

# Cast the spatial lines into a multi-linestring and sample
samp = st_sample(st_cast(ca$geometry[1], "MULTILINESTRING"), 1000)
tm_shape(samp) + 
  tm_dots()
print(samp)

# transform to correct projection
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
    dists$dist[i] = distm(c(ca_border_coord$X[i], ca_border_coord$Y[i]), c(s_lon, s_lat), fun = distHaversine)
  }
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
  }
  return(shore_distance)
}

# the function need a column called lat and long
ssf_data$long<-ssf_data$x2_
ssf_data$lat<-ssf_data$y2_

# run this - but takes a LONG TIME 
ssf_data$dist_to_shore_<-shark_shore_dist(ssf_data)
`
# ---- Create season variable ----

# first use summarised data to distinguish when thermocline is present
# filter to surface and deepest temperature
temp_range<-temp%>%filter(depth_m ==0 | depth_m == 10)

# summarise the mean daily temperature at each depth
temp_range_daily<-temp_range%>%group_by(DateTimeUTC,depth_m)%>%summarise(temp_c=mean(temp_c))

# transform from long to wide format
temp_range_daily_w<-temp_range_daily%>%pivot_wider(names_from = depth_m, values_from = temp_c,id_cols=c(DateTimeUTC))

# calculate the difference in temperature between surface and 10m
deltatemp_data<-temp_range_daily_w%>%mutate(delta_temp=`0`-`10`)

# create month and year variables
deltatemp_data$month_num<-month(deltatemp_data$DateTimeUTC)
deltatemp_data$year<-year(deltatemp_data$DateTimeUTC)

## Creating the homogenous and heterogenous layers
deltatemp_data_sum<-deltatemp_data%>%group_by(month_num,year)%>%summarise(mean(`0`),mean(`10`,na.rm=T),meanDT=mean(delta_temp,na.rm=T))

# If the difference is greater than 2.5 degrees then there is a thermocline, else its homogenous
deltatemp_data_sum$season<-NA
deltatemp_data_sum$season<-ifelse(deltatemp_data_sum$meanDT>2.5, "hetero", "cold")

# Combine this new variable with the previous SSF dataset
ssf_data$month_num <-month(ssf_data$t2_)
ssf_data$year<-year(ssf_data$t2_)

ssf_dataset<-left_join(ssf_data,deltatemp_data_sum,by=c("month_num", "year"))

# ---- Clean the dataset ----
# selecting relevant columns
ssf_dataset<-ssf_dataset%>%dplyr::select(sl_,ta_,case_,step_id,depth,temp_c,season,id, dist_to_shore_)

# Standardize all the variables
ssf_dataset<-as.data.frame(ssf_dataset)%>%ungroup()%>%mutate(st.WT=scale(temp_c),st.WT2=base::scale(temp_c^2),st.depth=base::scale(depth),st.depth2=base::scale(depth^2), logsl_=log(sl_),cos_ta_=cos(ta_),season=as.factor(season), st.distshore=base::scale(dist_to_shore_))

# add in numeric season variable
ssf_dataset$season_num<-as.numeric(as.factor(ssf_dataset$season))

# ---- Save dataset so don't have to re run the time expensive steps ----
saveRDS(ssf_dataset,"steps20.RDS")

# Now we can just load in ssf_dataset and start from here. 
ssf_dataset<-readRDS("steps20.RDS")

# ---- Fitting models ----

# ---- Stan model ----
model <- "
data {
  int<lower=1> N; // no data points
  int<lower=1> I; // no steps (over all individuals)
  int<lower=1> J; // no individuals
  int<lower=1> K; // number of predictors

  
  int<lower=0> y[N]; // response
  matrix[N,K] X; 

  int<lower=1, upper=N> stepid[N]; // step id
  int<lower=1, upper=J> indid[N]; // individual id
  int<lower=1,upper=2> season[N];
}

parameters {
  vector[2] beta_temp;
  vector[2] beta_depth;
  vector[2] beta_dist;
  vector[2] beta_temp_dist;
  vector[2] beta_depth_dist;
  vector[2] beta_temp_depth;
  vector[2] beta_three;
  vector[2] beta_sl;
  vector[2] beta_ta;
  vector<lower=0>[9] sigma;
  vector[I] a_re; // RE for steps
  vector[J] i_re1; // RE effects for each individual (temp)
  vector[J] i_re2;
  vector[J] i_re3;
  vector[J] i_re4;
  vector[J] i_re5;
  vector[J] i_re6;
  vector[J] i_re7;
  vector[J] i_re8;
  vector[J] i_re9;
}

model {
  vector[N] mu;
  
  // priors
  sigma ~ normal(0, 1);
  a_re ~ normal(0, 1000); // The paper has this as a_re ~ normal(0, 1000000)
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
  beta_temp_dist ~ normal(0,2);
  beta_depth_dist ~ normal(0,2);
  beta_temp_depth ~ normal(0,2);
  beta_three ~ normal(0,2);
  beta_sl ~ normal(0,2);
  beta_ta ~ normal(0,2);

  //likelihood
  mu = a_re[stepid]+
    (beta_temp[season] + i_re1[indid])  .* X[,1] + //temp effect 
    (beta_depth[season] + i_re2[indid]) .* X[,2] + // depth effect
    (beta_dist[season] + i_re3[indid]) .* X[,3] + // distance to shore effect
    (beta_temp_dist[season] + i_re4[indid]) .* X[,1] .* X[,3] + //temp.*distance to shore
    (beta_depth_dist[season] + i_re5[indid]) .* X[,2] .* X[,3]+ //depth.*distance to shore
    (beta_temp_depth[season] + i_re6[indid]) .* X[,1] .* X[,2] + //temp.*depth
    (beta_three[season] + i_re7[indid]) .* X[,1] .* X[,2] .* X[,3] + //temp.*depth.*distance to shore
    (beta_sl[season]+i_re8[indid]) .* X[,4] +
    (beta_ta[season]+i_re9[indid]) .* X[,5];
    
     y ~ poisson_log(mu);
}
"

# Create dependent variable matrix
names(ssf_dataset)
ssf_dataset$y<-as.numeric(ssf_dataset$case_)

X <- model.matrix(y ~ st.WT +st.depth + st.distshore + logsl_ + cos_ta_, data = ssf_dataset)
head(X)

# Package the data
stan_dat <- list(N = nrow(ssf_dataset), 
                 I = length(unique(ssf_dataset$step_id)), 
                 J = length(unique(ssf_dataset$id)), 
                 K = 5,
                 y = ssf_dataset$y, 
                 X = X[,2:6],
                 stepid = as.numeric(factor(ssf_dataset$step_id)), 
                 indid = as.numeric(factor(ssf_dataset$id)),
                 year=ssf_dataset$year,
                 season=as.numeric(as.factor(ssf_dataset$season)))


# Run the model, (this is SLOW, multiple days)
mod1<-stan(model_code = model,data=stan_dat,warmup=1000,iter=5000,chains=4,cores = 4,init_r = 0.1,pars=c("beta_temp","beta_depth","beta_dist","beta_temp_dist","beta_depth_dist","beta_temp_depth","beta_three","beta_sl","beta_ta","sigma","a_re"))

# View output
print(mod1)

# View traceplots
traceplot(mod1,pars=c("beta_temp","beta_depth","beta_dist","beta_temp_dist","beta_depth_dist","beta_temp_depth","beta_three","beta_sl","beta_ta","sigma"))
traceplot(mod1,pars=c("a_re[1]"))

# ---- GLMTMB (1) ----
# Analysis in Chatterjee et al. 2024

# Create ID variable as factor
ssf_dataset$id_f<-as.factor(ssf_dataset$id)

# This analysis needs to be done separately for each season
season1<-ssf_dataset%>%filter(season=="hetero")
season2<-ssf_dataset%>%filter(season=="cold")

# Heterogeneous season
shark.model.het <- glmmTMB(case_ ~ logsl_ + cos_ta_ + 
                          st.depth + st.WT + st.distshore + 
                          st.WT:st.distshore + st.depth:st.distshore + st.WT:st.depth + 
                          st.WT:st.distshore:st.depth +
                          (1|step_id) + 
                          (0+ st.depth|id) + (0+ st.WT|id) + (0+ st.distshore|id) +
                          (0+st.WT:st.distshore|id) + (0+st.depth:st.distshore|id) + (0 + st.WT:st.depth|id) +
                          (0+st.WT:st.distshore:st.depth|id),
                        REML = TRUE, family=poisson(), data = season1)
summary(shark.model.het)

shark.model.cold <- glmmTMB(case_ ~ logsl_ + cos_ta_ + 
                             st.depth + st.WT + st.distshore + 
                             st.WT:st.distshore + st.depth:st.distshore + st.WT:st.depth + 
                             st.WT:st.distshore:st.depth +
                             (1|step_id) + 
                             (0+ st.depth|id) + (0+ st.WT|id) + (0+ st.distshore|id) +
                             (0+st.WT:st.distshore|id) + (0+st.depth:st.distshore|id) + (0 + st.WT:st.depth|id) +
                             (0+st.WT:st.distshore:st.depth|id),
                           REML = TRUE, family=poisson(), data = season2)
summary(shark.model.cold)

# ---- GLMTMB (2) ----
# Analysis in Muff et al. 2020

MBStruc.fix.het = glmmTMB(case_ ~ logsl_ + cos_ta_ + 
                        st.depth + st.WT + st.distshore + 
                        st.WT:st.distshore + st.depth:st.distshore + st.WT:st.depth + 
                        st.WT:st.distshore:st.depth +
                        (1|step_id) + 
                        (0+ st.depth|id) + (0+ st.WT|id) + (0+ st.distshore|id) +
                        (0+st.WT:st.distshore|id) + (0+st.depth:st.distshore|id) + (0 + st.WT:st.depth|id) +
                        (0+st.WT:st.distshore:st.depth|id),
                      REML = TRUE, family=poisson(), data = season1, doFit=FALSE) 

#' Then fix the standard deviation of the first random term, which is the `(1|str_ID)` component  in the above model equation:
MBStruc.fix.het$parameters$theta[1] = log(5e5) # convergence problems start at 1e6

#' We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
MBStruc.fix.het$mapArg = list(theta=factor(c(NA,1,2,3,4,5,6,7)))

#' Then fit the model and look at the results:
glmm.TMB.fixed.het = glmmTMB:::fitTMB(MBStruc.fix.het) 
summary(glmm.TMB.fixed.het)


MBStruc.fix.cold = glmmTMB(case_ ~ logsl_ + cos_ta_ + 
                            st.depth + st.WT + st.distshore + 
                            st.WT:st.distshore + st.depth:st.distshore + st.WT:st.depth + 
                            st.WT:st.distshore:st.depth +
                            (1|step_id) + 
                            (0+ st.depth|id) + (0+ st.WT|id) + (0+ st.distshore|id) +
                            (0+st.WT:st.distshore|id) + (0+st.depth:st.distshore|id) + (0 + st.WT:st.depth|id) +
                            (0+st.WT:st.distshore:st.depth|id),
                          REML = TRUE, family=poisson(), data = season2, doFit=FALSE) 

#' Then fix the standard deviation of the first random term, which is the `(1|str_ID)` component  in the above model equation:
MBStruc.fix.cold$parameters$theta[1] = log(5e5) # convergence problems start at 1e6

#' We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
MBStruc.fix.cold$mapArg = list(theta=factor(c(NA,1,2,3,4,5,6,7)))

#' Then fit the model and look at the results:
glmm.TMB.fixed.cold = glmmTMB:::fitTMB(MBStruc.fix.cold) 
summary(glmm.TMB.fixed.cold)


## To do with whiteshark
# - make distance to shore calculation more efficient
# - determine if variance should be fixed at a large value or not (what value?)
# - stan vs glmmTMB (do they give similar answers?)
