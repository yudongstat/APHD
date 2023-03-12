rm(list=ls())
source("G:/My Drive/Topic3_phase type/paper/JRSSC Desktop/general purpose code/phase_type_estimation.R")
library(FAdist)

####################################Toy example 1: warranty data with sales lags
warr<-3;#warranty limit
censor<-9#end of study
sample_size<-100;#sample size
weib_par<-c(1.1,0.8)
llogis_par<-c(1.5,5.2)
set.seed(2023)
lags<-rweibull(sample_size,shape = weib_par[1],scale = weib_par[2])###generate sales lags
set.seed(2023)
ts<-rllog(sample_size,shape = 1/llogis_par[1],scale = log(llogis_par[2]))###generate product lifetime
ids<-which((lags+ts<censor)&(ts<warr))###determine which subject will be returned as warranty claims
lags<-lags[ids]###only keep the lags of those returned products
ts<-ts[ids]###only keep the lifetimes of those returned products
ncensor<-sample_size-length(ts)###number of unreturned products

####implement the proposed phase-type method
warranty_example<-phase_type_estimation(ts = ts,lags = lags, ncensor = ncensor, warranty_length = 1.5, end_of_study = 4.5)




####################################Toy example 2: Survival data with reporting delays
sample_size=100;#sample size
end_of_study<-7
weib_par<-c(2,6)
lnorm_par<-c(-0.4,0.6)
set.seed(2023)
ts<-rweibull(sample_size,shape = weib_par[1],scale = weib_par[2])###generate product lifetime
set.seed(2023)
lags<-rlnorm(sample_size,meanlog = lnorm_par[1],sdlog = lnorm_par[2])###generate reporting delays
ids=which(lags+ts<=end_of_study)###identify those patients who have been reported to the database before the end-of-study date
lags=lags[ids]###only keep the reporting delays of those reported patients
ts=ts[ids]###only keep the lifetime of those reported patients
ncensor=sample_size-length(ts)###number of unreported patients

####implement the proposed phase-type method
survival_example<-phase_type_estimation(ts = ts,lags = lags, ncensor = ncensor, warranty_length = Inf, end_of_study = end_of_study)
