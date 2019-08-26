wd = paste0(getwd(), '/Linda/MarginalFit/gamGPD/')
setwd(wd)
#### The code below displays the maps for the scale and shape parameters of the GPD for fixed month and year
load("data/all.excess.Rdata")
#this dataframe contains:
#datvar: exceedances with respect to a month and location-specific threshold (set as the 95% empirical quantile of the observations at a given location and a specific month)
#year: the year of the excess
#month: the month of the excess
#long: the longitude of the site at which the excess was observed (this is not used in the model but is needed to draw the plots)
#lat: the latitude of the site at which the excess was observed
#dist: the distance to the coast of the site at which the excess was observed. This distance is computed using the marmap package based on the bathymetric data from the NOAA server

load("data/mod_xi_time_lat.Rdata")
#this object contains a list with the following elements:
# xi: estimated shape
# beta: estimated beta (this is not needed in the gpd as it is a function of the scale and shape parameters)
# nu: estimated scale
# se.xi: standard error for xi
# se.nu: standard error for nu
# xi.covar: (unique) covariates for xi
# nu.covar: (unique) covariates for nu
# covar: *available* (not necessarily all) covariate combinations used for fitting beta (= xi *and* nu)
# y: the exceedances
# res: the residuals
# MRD: mean relative distances between old/new (xi, nu) for all iterations
# logL: log-likelihood at the estimated parameters
# xiObj: gamObject for estimated xi (return object of mgcv::gam())
# nuObj: gamObject for estimated nu (return object of mgcv::gam())
### The model for xi contains the effect of the time (year), the effect of the months (modelled with a cyclic cubic spline), and the effect of the latitude
### The model for nu contains the effect of the time (year), the effect of the months (modelled with a cyclic cubic spline), the effect of the latitude, the effect of the distance to the coast, and the effect of an interaction with the latitude and the distance to the coast
### All modelled effects are significant in both gams.

load("data/quant.mont.95.Rdata")
#this is a 12*16703 matrix
#each column represents the empirical 95% quantile of the monthly observations at a site 


require("ggplot2")
require("gridExtra")

### Plot for January at three different years

cond1 <- which((all.excess[,"month"]==1)&(all.excess[,"year"] == 1997))
dat <- data.frame("Longitude"=all.excess[cond1,"long"],"Latitude"=all.excess[cond1,"lat"],
                  "quant"=mod_xi_time_lat$xi[cond1])
plot_xi_1 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="shape")

dat <- data.frame("Longitude"=all.excess[cond1,"long"],"Latitude"=all.excess[cond1,"lat"],
                  "quant"=mod_xi_time_lat$beta[cond1])
plot_sigma_1 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="scale")

cond2 <- which((all.excess[,"month"]==1)&(all.excess[,"year"] == 2005))
dat <- data.frame("Longitude"=all.excess[cond2,"long"],"Latitude"=all.excess[cond2,"lat"],
                  "quant"=mod_xi_time_lat$xi[cond2])
plot_xi_2 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="shape")

dat <- data.frame("Longitude"=all.excess[cond2,"long"],"Latitude"=all.excess[cond2,"lat"],
                  "quant"=mod_xi_time_lat$beta[cond2])
plot_sigma_2 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="scale")

cond3 <- which((all.excess[,"month"]==1)&(all.excess[,"year"] == 2011))
dat <- data.frame("Longitude"=all.excess[cond3,"long"],"Latitude"=all.excess[cond3,"lat"],
                  "quant"=mod_xi_time_lat$xi[cond3])
plot_xi_3 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="shape")

dat <- data.frame("Longitude"=all.excess[cond3,"long"],"Latitude"=all.excess[cond3,"lat"],
                  "quant"=mod_xi_time_lat$beta[cond3])
plot_sigma_3 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="scale")


grid.arrange(plot_sigma_1,plot_xi_1,
             plot_sigma_2,plot_xi_2,
             plot_sigma_3,plot_xi_3, ncol=2)

### Plot for July at three different years

cond1 <- which((all.excess[,"month"]==7)&(all.excess[,"year"] == 1997))
dat <- data.frame("Longitude"=all.excess[cond1,"long"],"Latitude"=all.excess[cond1,"lat"],
                  "quant"=mod_xi_time_lat$xi[cond1])
plot_xi_1 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="shape")

dat <- data.frame("Longitude"=all.excess[cond1,"long"],"Latitude"=all.excess[cond1,"lat"],
                  "quant"=mod_xi_time_lat$beta[cond1])
plot_sigma_1 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="scale")

cond2 <- which((all.excess[,"month"]==7)&(all.excess[,"year"] == 2005))
dat <- data.frame("Longitude"=all.excess[cond2,"long"],"Latitude"=all.excess[cond2,"lat"],
                  "quant"=mod_xi_time_lat$xi[cond2])
plot_xi_2 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="shape")

dat <- data.frame("Longitude"=all.excess[cond2,"long"],"Latitude"=all.excess[cond2,"lat"],
                  "quant"=mod_xi_time_lat$beta[cond2])
plot_sigma_2 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="scale")

cond3 <- which((all.excess[,"month"]==7)&(all.excess[,"year"] == 2011))
dat <- data.frame("Longitude"=all.excess[cond3,"long"],"Latitude"=all.excess[cond3,"lat"],
                  "quant"=mod_xi_time_lat$xi[cond3])
plot_xi_3 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="shape")

dat <- data.frame("Longitude"=all.excess[cond3,"long"],"Latitude"=all.excess[cond3,"lat"],
                  "quant"=mod_xi_time_lat$beta[cond3])
plot_sigma_3 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="scale")


grid.arrange(plot_sigma_1,plot_xi_1,
             plot_sigma_2,plot_xi_2,
             plot_sigma_3,plot_xi_3, ncol=2)

### Plot for October at three different years

cond1 <- which((all.excess[,"month"]==10)&(all.excess[,"year"] == 1997))
dat <- data.frame("Longitude"=all.excess[cond1,"long"],"Latitude"=all.excess[cond1,"lat"],
                  "quant"=mod_xi_time_lat$xi[cond1])
plot_xi_1 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="shape")

dat <- data.frame("Longitude"=all.excess[cond1,"long"],"Latitude"=all.excess[cond1,"lat"],
                  "quant"=mod_xi_time_lat$beta[cond1])
plot_sigma_1 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="scale")

cond2 <- which((all.excess[,"month"]==10)&(all.excess[,"year"] == 2005))
dat <- data.frame("Longitude"=all.excess[cond2,"long"],"Latitude"=all.excess[cond2,"lat"],
                  "quant"=mod_xi_time_lat$xi[cond2])
plot_xi_2 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="shape")

dat <- data.frame("Longitude"=all.excess[cond2,"long"],"Latitude"=all.excess[cond2,"lat"],
                  "quant"=mod_xi_time_lat$beta[cond2])
plot_sigma_2 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="scale")

cond3 <- which((all.excess[,"month"]==10)&(all.excess[,"year"] == 2011))
dat <- data.frame("Longitude"=all.excess[cond3,"long"],"Latitude"=all.excess[cond3,"lat"],
                  "quant"=mod_xi_time_lat$xi[cond3])
plot_xi_3 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="shape")

dat <- data.frame("Longitude"=all.excess[cond3,"long"],"Latitude"=all.excess[cond3,"lat"],
                  "quant"=mod_xi_time_lat$beta[cond3])
plot_sigma_3 <- ggplot() + 
  geom_polygon(data=dat, aes(x=Longitude, y=Latitude),colour="white", fill="white")+
  geom_point(data=dat,aes(x=Longitude, y=Latitude,colour=quant))+
  scale_colour_gradientn(colours = rev(heat.colors(10)),name="scale")


grid.arrange(plot_sigma_1,plot_xi_1,
             plot_sigma_2,plot_xi_2,
             plot_sigma_3,plot_xi_3, ncol=2)

#### Example of how to get (or predict if there were 0 observations for the combination of covariates) the values of the scale and shape for a fixed combination of covariates
# this give the shape for the site 1 (lat=29.48 and distance to coast=9147.422) in January 1995
predict(mod_xi_time_lat$xiObj,newdata=data.frame("month"=1,"year"=1995,"lat"=29.48,"dist"=9147.422))
# this give the scale for the site 1 (lat=29.48 and distance to coast=9147.422) in January 1995
exp(predict(mod_xi_time_lat$nuObj,newdata=data.frame("month"=1,"year"=1995,"lat"=29.48,"dist"=9147.422)))/
  (1 + predict(mod_xi_time_lat$xiObj,newdata=data.frame("month"=1,"year"=1995,"lat"=29.48,"dist"=9147.422)))

### The code below displays the probability of exceedance for the same fixed month and years as above

load("data/nbr_excess.Rdata")
#this dataframe contains:
#n.u: the number of exceedances with respect to a month and location-specific threshold (set as the 95% empirical quantile of the observations at a given location and a specific month)
#n: the number of observations for a specific month, year and location
#month: the month of the excess
#year: the year of the excess
#long: the longitude of the site at which the excess was observed (this is not used in the model but is needed to draw the plots)
#lat: the latitude of the site at which the excess was observed
#dist: the distance to the coast of the site at which the excess was observed. This distance is computed using the marmap package based on the bathymetric data from the NOAA server

load("data/lam_full_logit.Rdata")
#this is a standard gam object where the rate of exceedances is modelled using a logistic regression with all the covariates at hand, i.e., the year, the month, the latitude, the distance and the interaction between the latitude and the distance,
# require(mgcv); summary(lam_full_logit) displays the significance of all the considered effects
# the gam was asked to select only the significant covariates and to reduce to (almost) 0 the degrees of freedom for the non-significant ones
# the effect of the distance alone (not interacting with latitude) was reduced to almost 0 in the final model

#### Example of how to get (or predict if there were 0 observations for the combination of covariates) the rate of exceedance for a fixed combination of covariates

predict(lam_full_logit,newdata=data.frame("month"=1,"year"=1995,"lat"=29.48,"dist"=9147.422),type="response")

### get fitted values for all available combinations of covariates
prob.fit <- cbind(nbr_excess, prob=lam_full_logit$fitted.values)
#save(prob.fit,file='prob.fit.Rdata')

### The code below displays the rate of the non homogeneous Poisson process followed by the number of exceedances for the same fixed month and years as above

load("data/lam_full.Rdata")
#this is a standard gam object where the number of exceedances is modelled using a Poisson regression with all the covariates at hand, i.e., the year, the month, the latitude, the distance and the interaction between the latitude and the distance,
# require(mgcv); summary(lam_full) displays the significance of all the considered effects
# the gam was asked to select only the significant covariates and to reduce to (almost) 0 the degrees of freedom for the non-significant ones
# the effects of all the covariates were significant

#### Example of how to get the number of exceedance for a fixed combination of covariates

predict(lam_full,newdata=data.frame("month"=1,"year"=1995,"lat"=29.48,"dist"=9147.422),type="response")

### get fitted values for all available combinations of covariates
Lambda.fit <- cbind(nbr_excess, Lambda=lam_full$fitted.values)
#save(Lambda.fit,file='Lambda.fit.Rdata')
