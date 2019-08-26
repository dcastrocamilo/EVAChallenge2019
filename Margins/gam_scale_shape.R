########################################################################################################################
########################################################################################################################
library(QRM)
library(mgcv)
library(parallel)
source("~/Dropbox/DanielaLindaThomas/gamGPD/game.R")

### 1) Fitting the GPD model #################################################
load("~/Dropbox/DanielaLindaThomas/gamGPD/data/all.excess.Rdata")
#Note: we are feeding the gam with the exceedances => the threshold will be set at 0

## simplest model (constant GPD parameters)
mod0 <- gamGPDfit(all.excess, threshold=0, datvar="datvar", xiFrhs=~1, nuFrhs=~1, eps.xi = 2e-5)

### model with trend and effect of months
mod_time <- gamGPDfit(all.excess, threshold=0, datvar="datvar", xiFrhs=~s(year)+s(month,bs ="cc",k= 12), 
                      nuFrhs=~s(year)+s(month,bs ="cc",k= 12))

### full model (both scale and shape vary with covariates)
###### Model for the scale and shape parameters
xiFrhs= ~s(year)+s(month,bs ="cc",k= 12)+ti(lat)+ti(dist)+ti(lat,dist)
nuFrhs= ~s(year)+s(month,bs="cc",k = 12)+ti(lat)+ti(dist)+ti(lat,dist)
mod_full <- gamGPDfit(all.excess, threshold=0, datvar="datvar", xiFrhs=xiFrhs, nuFrhs=nuFrhs)
#this model did not converge due to the non-convergence of the model for the shape parameter 

### full model for the scale & model for the shape as a function of years, months, and latitude
###### Model for the scale and shape parameters
xiFrhs= ~s(year)+s(month,bs ="cc",k= 12)+s(lat)
nuFrhs= ~s(year)+s(month,bs="cc",k = 12)+ti(lat)+ti(dist)+ti(lat,dist)
mod_xi_time_lat <- gamGPDfit(all.excess, threshold=0, datvar="datvar", xiFrhs=xiFrhs, nuFrhs=nuFrhs)
save(mod_xi_time_lat,file="mod_xi_time_lat.Rdata")
# all covariates had significant effect on the shape and the scale parameters

## final model
modGPD <- mod_xi_time_lat

### 2) Conduct the bootstrap  #################################################
## conduct the bootstrap and save the object
B=300 #number of replicates

xiFrhs= ~s(year)+s(month,bs ="cc",k= 12)+s(lat)
nuFrhs= ~s(year)+s(month,bs ="cc",k= 12)+ti(lat)+ti(dist)+ti(lat,dist)

sfile <- "modGPDboot.rds"
if(file.exists(sfile)){
  modGPDboot <- readRDS(sfile)
} else {
  set.seed(1) # set seed to be reproducible
  system.time(modGPDboot <- gamGPDboot(all.excess, B=B, threshold=0, datvar="datvar", xiFrhs=xiFrhs, nuFrhs=nuFrhs,
                                       niter=150, include.updates=FALSE, eps.xi=2e-5, eps.nu=1e-5,
                                       boot.progress=TRUE, verbose=FALSE))     
  ## save the bootstrapped object
  saveRDS(modGPDboot, file=sfile)
}


### 3) Fitting the frequency of exceedance (standard gam() application) #########
load("~/Dropbox/DanielaLindaThomas/gamGPD/data/nbr_excess.Rdata")

## year as a covariate
cl <- makeCluster(detectCores())
lam1 <- bam(n.u~s(year), data=nbr_excess, family=poisson,cluster=cl) # fit gam model
stopCluster(cl)
save(lam1,file="lam1.Rdata")

## year and month as covariates
cl <- makeCluster(detectCores())
lam2 <- bam(n.u~s(year)+s(month,bs ="cc",k= 12), data=nbr_excess, family=poisson,cluster=cl) # fit gam model
stopCluster(cl)
save(lam2,file="lam2.Rdata")

## check the effect of month (additionally to year)
lr <- as.numeric(-2*(logLik(lam1)-logLik(lam2))) # likelihood ratio test statistic
1-pchisq(lr, df=1) # p-value of the likelihood ratio test
## => month is also significant

## year, month, and latitude as covariates
cl <- makeCluster(detectCores())
lam3 <- bam(n.u~s(year)+s(month,bs ="cc",k= 12)+s(lat), data=nbr_excess, family=poisson,cluster=cl) # fit gam model
stopCluster(cl)
save(lam3,file="lam3.Rdata")

## check the effect of latitude (additionally to year and month)
lr <- as.numeric(-2*(logLik(lam2)-logLik(lam3))) # likelihood ratio test statistic
1-pchisq(lr, df=1) # p-value of the likelihood ratio test
## => latitude is also significant

## full model (year, month, latitude, distance, and interaction of latitude and distance)
cl <- makeCluster(detectCores())
lam_full <- bam(n.u~s(year)+s(month,bs ="cc",k= 12)+ti(lat)+ti(dist)+ti(lat,dist), data=nbr_excess, select=TRUE, family=poisson,cluster=cl) # fit gam model
stopCluster(cl)
save(lam_full,file="lam_full.Rdata")
## we keep this model as the gam is asked to reduce the degrees of freedom of non significant covariates to null (select=TRUE)

### Predict the occurence of exceedance at site 1 in January 2010 (for example) #################

lambda_site1_jan_2010 <- predict(lam_full,newdata = data.frame("year"=1985,"month"=1,"dist"=9147.422,"lat"=29.48),se.fit=TRUE)
exp(lambda_site1_jan_2010$fit)

### 4) Fitting the rate of the number of exceedances using logistic regression (standard gam() application) #########
# this is for prediction purposes: for a combination of the covariates where no excess has been observed, we could still predict of probability of exceedance
load("~/Dropbox/DanielaLindaThomas/gamGPD/data/nbr_excess.Rdata")

## full model (year, month, latitude, distance, and interaction of latitude and distance)
cl <- makeCluster(detectCores())
lam_full_logit <- bam(n.u~s(year)+s(month,bs ="cc",k= 12)+ti(lat)+ti(dist)+ti(lat,dist), data=nbr_excess, select=TRUE, link=logit,cluster=cl) # fit gam model
stopCluster(cl)
save(lam_full_logit,file="lam_full_logit.Rdata")

