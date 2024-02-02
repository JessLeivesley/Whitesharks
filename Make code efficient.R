##### 20th Decemeber #####

# Building up the model so that I know it converges and runs in a reasonable amount of time.
# so far I have vectorized the code and narrowed the a_re and sigma priors. 
# once I've built this up to the most complex model I want then im going to increase the a_re sigma. 

# so far with main effects of temp, depth and distance to shore, temp*dist and dist*depth with normal(0,1)'s the model takes ~6637 seconds to run 1000 iterations. (ESS low warning though)


model <- "
data {
  int<lower=1> N; // no data points
  int<lower=1> I; // no steps (over all individuals)
  int<lower=1> J; // no individuals
  int<lower=1> K; // number of predictors

  
  int<lower=0> y[N]; // response
  matrix[N,K] X; // temperature
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
    (beta_temp_dist[season] + i_re4[indid]) .* X[,1] .* X[,3] + //temp*distance to shore
    (beta_depth_dist[season] + i_re5[indid]) .* X[,2] .* X[,3]+ //depth*distance to shore
    (beta_temp_depth[season] + i_re6[indid]) .* X[,1] .* X[,2] + //temp*depth
    (beta_three[season] + i_re7[indid]) .* X[,1] .* X[,2] .* X[,3]+ //temp*depth*distance to shore
    (beta_sl[season]+i_re8[indid]) .* X[,4] +
    (beta_ta[season]+i_re9[indid]) .* X[,5];
    
     y ~ poisson_log(mu);
}
"
X <- model.matrix(y ~ st.WT +st.depth + st.distshore + logsl_ + cos_ta_, data = ssf_data)
head(X)

stan_dat <- list(N = nrow(ssf_data), I = length(unique(ssf_data$step_id)), 
                 J = length(unique(ssf_data$id)), 
                 K = 5,
                 y = ssf_data$y, 
                 X=as.matrix(X[,2:6]),
                stepid = as.numeric(factor(ssf_data$step_id)), 
                 indid = ssf_data$id,
                 year=ssf_data$year,
                 season=as.numeric(as.factor(ssf_data$season)))

initf1 <- function() {
  list(beta_temp = rnorm(2,0,0.01),beta_depth = rnorm(2,0,0.01),beta_dist = rnorm(2,0,0.01), beta_temp_dist = rnorm(2,0,0.01), beta_depth_dist = rnorm(2,0,0.01), sigma = runif(9,0,1), i_re1 = rnorm(21,0,0.01),i_re2 = rnorm(21,0,0.01),i_re3 = rnorm(21,0,0.01),i_re4=rnorm(21,0,0.01), i_re5 = rnorm(21,0,0.01),i_re6 = rnorm(21,0,0.01),i_re7 = rnorm(21,0,0.01),i_re8 = rnorm(21,0,0.01),i_re9 = rnorm(21,0,0.01), a_re = rnorm(6495,0,0.01))
}


mod1<-stan(model_code = model,data=stan_dat,warmup=500,iter=8000,chains=2,cores = 2,init = initf1,pars=c("beta_temp","beta_depth","beta_dist","beta_temp_dist","beta_depth_dist","beta_temp_depth","beta_three","beta_sl","beta_ta","sigma","a_re"))
print(mod1)
pairs(mod1,pars=c("sigma","beta_temp"))

saveRDS(mod1,"fullmodel_19thJan.RDS")

which(summary(mod1)$summary[,10]>1.05)

traceplot(mod1,pars=c("sigma"))
traceplot(mod1,pars=c("beta_temp","beta_depth","beta_dist"))
traceplot(mod1,pars=c("beta_temp_dist","beta_depth_dist","beta_temp_depth"))
traceplot(mod1,pars=c("beta_three"))


#Relative selection strength for Temperature across seasons, depths and distances to shore
# 1 = cold homo; 2 = hetero, 3 = hot homo
beta_temp<-rstan::extract(mod1,pars="beta_temp")
beta_depth<-rstan::extract(mod1,pars="beta_depth")
beta_dist<-rstan::extract(mod1,pars="beta_dist")
beta_tdist_int<-rstan::extract(mod1,pars="beta_temp_dist")
beta_ddist_int<-rstan::extract(mod1,pars="beta_depth_dist")



reftemp<-(20-mean(ssf_data$watertemp_))/sd(ssf_data$watertemp_)
temp<-seq(11.9,23,length=25)
st.temp<-(temp-mean(ssf_data$watertemp_))/sd(ssf_data$watertemp_)
depth<-seq(0,10,by=2)
st.depth<-(depth-mean(ssf_data$depth_m))/sd(ssf_data$depth_m)
dist<-c(10,210,410,610,810,1010)
st.dist<-(dist-mean(ssf_data$dist_to_shore_))/sd(ssf_data$dist_to_shore_)

dim(beta_temp$beta_temp)
dim(beta_temp$beta_temp)

rss.pred<-array(dim=c(3000,25,3,6,6))
for(i in 1:3000){
  for(t in 1:25){
    for(s in 1:3){
      for(d in 1:6){
        for(dist in 1:6){
          rss.pred[i,t,s,d,dist]<-beta_temp$beta_temp[i,s]*st.temp[t]+beta_depth$beta_depth[i,s]*st.depth[d]+beta_dist$beta_dist[i,s]*st.dist[dist]+
            beta_ddist_int$beta_depth_dist[i,s]*st.temp[t]*st.dist[dist]+beta_tdist_int$beta_temp_dist[i,s]*st.temp[t]*st.dist[dist]
        }}}}}

rss.ref<-array(dim=c(3000,3,6,6))
for(i in 1:3000){
  for(s in 1:3){
    for(d in 1:6){
      for(dist in 1:6){
        rss.ref[i,s,d,dist]<-beta_temp$beta_temp[i,s]*reftemp+beta_depth$beta_depth[i,s]*st.depth[d]+
          beta_dist$beta_dist[i,s]*st.dist[dist]+
          beta_ddist_int$beta_depth_dist[i,s]*st.depth[d]*st.dist[dist]+
          beta_tdist_int$beta_temp_dist[i,s]*reftemp*st.dist[dist]
      }}}}

rss<-array(dim=c(3000,25,3,6,6))
for(i in 1:3000){
  for(t in 1:25){
    for(s in 1:3){
      for(d in 1:6){
        for(dist in 1:6){
          rss[i,t,s,d,dist]<-exp(rss.pred[i,t,s,d,dist]-rss.ref[i,s,d,dist])  
        }}
    }}}


plot_list<-list()


depthindex<-rep(1:6,6)
distindex<-rep(1:6,each=6)

var_list<-list()

for(i in 1:36){
  var_list[[i]]<-c(depthindex[i],distindex[i])
}


for(i in 1:36){
  p = ggplot()+
   geom_line(aes_string(y=colMeans(rss[,,1,var_list[[i]][1],var_list[[i]][2]]),x=st.temp),col="blue")+
    geom_line(aes_string(y=colMeans(rss[,,2,var_list[[i]][1],var_list[[i]][2]]),x=st.temp),col="purple")+
      geom_line(aes_string(y=colMeans(rss[,,3,var_list[[i]][1],var_list[[i]][2]]),x=st.temp),col="red")+
   geom_ribbon(aes_string(x=st.temp,ymin=colQuantiles(rss[,,1,var_list[[i]][1],var_list[[i]][2]],prob=0.025),ymax=colQuantiles(rss[,,1,var_list[[i]][1],var_list[[i]][2]],prob=0.975)),alpha=0.25,fill="blue")+
  geom_ribbon(aes_string(x=st.temp,ymin=colQuantiles(rss[,,2,var_list[[i]][1],var_list[[i]][2]],prob=0.025),ymax=colQuantiles(rss[,,2,var_list[[i]][1],var_list[[i]][2]],prob=0.975)),alpha=0.25,fill="purple")+
      geom_ribbon(aes_string(x=st.temp,ymin=colQuantiles(rss[,,3,var_list[[i]][1],var_list[[i]][2]],prob=0.025),ymax=colQuantiles(rss[,,3,var_list[[i]][1],var_list[[i]][2]],prob=0.975)),alpha=0.25,fill="red")+
    theme_classic()+
    geom_hline(aes(yintercept=1),lty=2)+
    xlab("Water Temperature")
  
  
  plot_list[[i]] = p
}



(plot_list[[1]]+ggtitle("0m Depth")+ylab("10m
Offshore
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
  (plot_list[[7]]+ylab("210m
Offshore
RSS")+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())|plot_list[[8]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[9]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[10]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[11]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[12]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL))/
  (plot_list[[13]]+ylab("410m
Offshore
RSS")+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())|plot_list[[14]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[15]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[16]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[17]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[18]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL))/
  (plot_list[[19]]+ylab("610m
Offshore
RSS")+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())|plot_list[[20]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[21]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[22]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[23]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[24]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL))/
  (plot_list[[25]]+ylab("810m
Offshore
RSS")+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())|plot_list[[26]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[27]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[28]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[29]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL)|plot_list[[30]]+
     theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank())+ylab(NULL))/
  (plot_list[[31]]+ylab("1010m
Offshore
RSS")+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))|plot_list[[32]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL)
   |plot_list[[33]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL)|plot_list[[34]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL)|plot_list[[35]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL)|plot_list[[36]]+scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c(12.6,14.8,17,19.2,21.4))+ylab(NULL))


