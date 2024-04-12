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
print(dat1%>% dplyr::select(id,data)%>%unnest(cols=c(data))%>%group_by(id)%>%summarise(mindepth=min(depth),maxdepth=max(depth)),n=21)

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

set.seed(123)
for(i in 1:21){
  print(i)
  
  dat_ssf2$rando[[i]]<-random_steps(dat_ssf2$stps[[i]],
                                    n_control = 300,
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
ssf<-extract_covariates(obs_avail,rp)

#remove moves outside of our array
ssf_filt<-ssf%>%dplyr::filter(is.na(ssf$layer)==F)

# how many random moves per move am I left with? 
sum_ssf<-as.data.frame(ssf_filt)%>%filter(case_==FALSE)%>%group_by(step_id_)%>%count()
range(sum_ssf$n) 

## ---- Now add depth to true and simulated movement ----

# Get the max depth per polygon
# this will be used to simulate the potential depth the individual could be at
temp<-readRDS("available_temp_2020-2021_clean.RDS")
Polygon_depth<-temp%>%group_by(PolygonID_250m)%>%summarise(MaxDepth=mean(Total_Water_Colum_m))

# Combine total water column depth with all moves
ssf_totaldepth<-left_join(ssf_filt,Polygon_depth,by=c("layer"="PolygonID_250m"))

# get depth of all movements
ssf_depth<-(left_join(ssf_totaldepth,shark2020_20,by=c("x2_"="x_","y2_"="y_","t2_"="t_")))

#Polygon 54 doesnt exist, and Polygon 55 has negative water column depth??
ssf_depth<-ssf_depth%>%filter(layer!=55)
ssf_depth<-ssf_depth%>%filter(layer!=54)

# Simulate depth
set.seed(12)
ssf_depth$depth<-ifelse(ssf_depth$case_==TRUE,ssf_depth$depth,dunif(1,0,ssf_depth$MaxDepth))

# Now create depth bin
ssf_depth$depth_bin2m

# label to match temperature depth categories
ssf_depth$depthbin<-ifelse(ssf_depth$depth<1,"0-1m",
       ifelse(ssf_depth$depth<3,"2-3m",
              ifelse(ssf_depth$depth<5,"4-5m",
                     ifelse(ssf_depth$depth<7,"6-7m",
                            ifelse(ssf_depth$depth<9,"8-9m","10-11m")))))

# create rounded hour
ssf_depth$t2_h<-round(ssf_depth$t2_,units="hours")

# now combine temperature based on hour
ssf_full<-left_join(ssf_depth,temp,by=c("layer"="PolygonID_250m","depthbin"="depth_bin2m","t2_h"="DateTimeUTC"))

# There are 404 missing temperature records - where the temperature array wasnt measure temp yet!
ssf_full<-ssf_full%>%filter(is.na(temp_c)==F)

# how many random moves per move am I left with? 
sum_ssf<-as.data.frame(ssf_full)%>%filter(case_==FALSE)%>%group_by(step_id_)%>%count()
range(sum_ssf$n)

# Keep the 100 random moves per step
summary(ssf_full$step_id_)
ssf_data<-as.data.frame(ssf_full)%>%group_by(step_id_)%>%dplyr::filter(row_number()<100)%>%ungroup()

##---- Create distance to shore variable for each true and proposed movement

# Read in shoreline shapefile
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

# function need a column called lat and long
ssf_data$long<-ssf_data$x2_
ssf_data$lat<-ssf_data$y2_

# takes a LONG TIME 
ssf_data$dist_to_shore_<-shark_shore_dist(ssf_data)


##---- Create season variable
# filter to surface and deepest temperature
temp_range<-temp%>%filter(depth_m ==0 | depth_m == 10)

# summarise the mean daily temperature at each depth
temp_range_daily<-temp_range%>%group_by(DateTimeUTC,depth_m)%>%summarise(temp_c=mean(temp_c))

# transform from long to wide format
temp_range_daily_w<-temp_range_daily%>%pivot_wider(names_from = depth_m, values_from = temp_c,id_cols=c(DateTimeUTC))

# calculate the difference in temperature between surface adn 10m
deltatemp_data<-temp_range_daily_w%>%mutate(delta_temp=`0`-`10`)

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

#select relevant columns
ssf_dataset<-ssf_dataset%>%dplyr::select(sl_,ta_,case_,step_id_,depth,temp_c,season,dist_to_shore_)

ssf_dataset<-as.data.frame(ssf_dataset)%>%ungroup()%>%mutate(st.WT=scale(temp_c),st.WT2=base::scale(temp_c^2),st.depth=base::scale(depth),st.depth2=base::scale(depth^2),st.distshore=base::scale(dist_to_shore_), logsl_=log(sl_),cos_ta_=cos(ta_),season=as.factor(season))
levels(ssf_dataset$season)
summary(ssf_dataset$season)
ssf_dataset$season_num<-as.numeric(as.factor(ssf_dataset$season))

##---- Fit SSF model with whatever distribution you want ----
m1<- fit_issf(ssf_dataset, 
              # Response
              case_ ~ 
                st.WT + st.depth + st.distshore+season+
                # Stratum
                strata(step_id_),
              # Need this later for model predictions
              model = TRUE)
summary(m1)
