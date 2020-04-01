rm(list=ls())
library(mgcv)
library(parallel)
library(doParallel)
library(marmap)

### Read data
load("DATA_TRAINING.RData")

#################################
### Compute distance to coast   #
#################################
tmp <- getNOAA.bathy(lon1 = 32, lon2 = 44,
                     lat1 = 12, lat2 = 30.5, resolution = 10)


# Create vectors of latitude and longitude
lon <- loc[,1]
lat <- loc[,2]
# Compute distances (in meters) between each point and the coastline
# this takes a while...
dist2coast <- dist2isobath(tmp, lon, lat, isobath = 0,locator=FALSE)

###############################################################################
### 1) Fitting the Gaussian model (1) to capture the long-term trend          #
###############################################################################

### Data preparation
all.obs   = as.vector(anom.training)
all.site  = rep(1:ncol(anom.training), each= nrow(anom.training))
all.dist  = rep(dist2coast[,1], each= nrow(anom.training))
all.lat   = rep(loc[,2], each= nrow(anom.training))
all.month = rep(month, ncol(anom.training))
all.year  = rep(year, ncol(anom.training))

dat= data.frame("obs"=all.obs,
                "site"=all.site,
                "dist"=all.dist,
                "lat"=all.lat,
                "month"=all.month,
                "year"=all.year)

###### we use only the long-term trend in the Gaussian model
###### we keep one site out of 10 and model both the mean and the scale of the Gaussian model

# Keep one site out of 10: selection made randomly
set.seed(22)
site2keep        <- sample(1:ncol(anom.training), round(ncol(anom.training)/50))

mod_Gaussls_year <- gam(list(obs~-1+as.factor(year),~-1+as.factor(year)), 
                        data=dat, subset= which(dat$site %in% site2keep),
                        family=gaulss())
###############################################################################
#fitted mean and sd for each of the 31 years

fitted.Gauss <- predict(mod_Gaussls_year, type="response")

fitted.mu <- rep(unique(fitted.Gauss[,1]), each=365)
fitted.sd <- rep(1/unique(fitted.Gauss[,2]), each=365)

mu_mat <- matrix(NA, ncol=ncol(anom.training), nrow=nrow(anom.training))
for(j in 1:ncol(mu_mat)){
  mu_mat[,j] <- fitted.mu
}

sd_mat <- matrix(NA, ncol=ncol(anom.training), nrow=nrow(anom.training))
for(j in 1:ncol(mu_mat)){
  sd_mat[,j] <- fitted.sd
}

#get the normalized dataset \tilde{Z}(s,t)
anom.training.gauss <- (anom.training-mu_mat)/sd_mat

#############################################################################################
### In what follows, we model the upper right tail using either a GAM specification   #######
### or a local approach to capture the spatial non-stationarity in the GP parameters  #######
### First, fix the threshold at 0.75 (around the 78.3% empirical quantile)            #######
#############################################################################################

### Threshold exceedances
exc.gauss     <- c(anom.training.gauss-0.75) #stacked site by site, i.e., exc of site 1, exc of site 2, etc...
dat.exc.gauss <- data.frame("exc"=exc.gauss,
                            "lat"=rep(loc$lat, each=nrow(anom.training)),
                            "lon"=rep(loc$lon, each=nrow(anom.training)),
                            "dist"=rep(dist2coast$distance, each=nrow(anom.training)))

########################################################################
############                The GAM approach                  ##########
########################################################################

############           subsample from exceedances 
set.seed(22)
subsample_tail.exc <- sample(which(dat.exc.gauss$exc>0), 
                             round(length(which(dat.exc.gauss$exc>0))/10))

dat.exc.gauss.subsample.exc <- dat.exc.gauss[subsample_tail.exc,]

library(evgam)
dat.exc.gauss_gpd <- dat.exc.gauss.subsample.exc
dat.exc.gauss_gpd$exc[dat.exc.gauss_gpd$exc <= 0] <- NA

#############################################################################################################
### 2) The GAM specification with a latitude-longitude interaction and distance to coast in both parameters #
#############################################################################################################
fmla_gpd                   <- list(exc ~ s(lat, lon)+ s(dist), ~ s(lat, lon)+ s(dist))
m_gpd_lat_lon_dist_sub_exc <- evgam(fmla_gpd, dat.exc.gauss_gpd, family = "gpd")

#Estimated scale and shape parameters for the 16703 locations in the dataset
dat2pred          <- data.frame("lat"=loc$lat, "lon"=loc$lon, "dist"=dist2coast$distance)
gpd_lat_dist_pred <- predict(m_gpd_lat_lon_dist_sub_exc, dat2pred, type = "response", se.fit = TRUE)
gpd.fitted        <- gpd_lat_dist_pred$fitted
gpd.fitted.se     <- gpd_lat_dist_pred$se.fit
gpd_estimates     <- list("fitted_gpd"=gpd.fitted, "sd_gpd"=gpd.fitted.se)

########################################################################################################################
### 3) Fitting the spatial pattern in the probability of exceedance
########################################################################################################################

### dataset containing the number of exceedances and the total number of observations per site
### as well as the latitude, longitude, and distance to coast of each site

exc.gauss.mat <- anom.training.gauss-0.75
fct_extract <- function(s){
  nbr.excess  <- data.frame("n.u"=length(which(exc.gauss.mat[,s]>0)),
                            "n"=length(exc.gauss.mat[!is.na(exc.gauss.mat[,s]),s]),
                            "lon"=loc[s,1],
                            "lat"=loc[s,2],
                            "dist"=dist2coast[s,"distance"])
  return(nbr.excess)
}

registerDoParallel(cores = detectCores())
nbr_excess_s <- foreach(s=1:ncol(anom.training.gauss)) %dopar% {
  fct_extract(s)
}
stopImplicitCluster()

nbr_excess <- do.call(rbind, nbr_excess_s)


cl <- makeCluster(detectCores())
prob.exc.model <- bam((n.u)~s(lat,lon)+s(dist), data=nbr_excess,
                      family=poisson,select=TRUE,cluster=cl)
stopCluster(cl)

### get fitted values for all available combinations of covariates, i.e., all locations
prob.fit <- cbind(nbr_excess, prob=prob.exc.model$fitted.values/nbr_excess$n)

# Finally, to obtain vectors containing one estimate for each location in loc:
scale = gpd_estimates$fitted_gpd$scale
xi = gpd_estimates$fitted_gpd$shape
pexc = prob.fit$prob

########################################################################
############                The local approach                ##########
########################################################################

# get 40 nearest neighbors of each pixel ####
library(FNN)
k.nn = 40
u = 0.75
nn = get.knn(loc, k = k.nn)$nn.index
dim(nn)

# functions for maximum likelihood estimation of the GPD ####

# GPD density
function(x,sigma,xi){
  if(abs(xi)<10^{-4}){
    exp(-x/sigma)/sigma
  }else{
    ifelse(1+xi*x/sigma<=0,0,sigma^{-1}*(1+xi*x/sigma)^{-1/xi-1})
  }
}
# negative GP log-likelihood ####
fun2opt=function(par, exc){
  if(par[1] <= 0) return(Inf)
  -sum(log(dgp(x=exc, sigma = par[1], xi = par[2])))
}

# estimate marginal parameters using anom.training.gauss defined above ####
xi = scale = pexc = rep(NA, ncol(anom.training.gauss))
for(i in 1:ncol(anom.training.gauss)){
  sample.i = as.numeric(anom.training.gauss[, unique(c(i, nn[i,]))])
  pexc[i] = mean(sample.i > u, na.rm = TRUE)
  exc.i = na.omit(sample.i[sample.i > u] - u + 10^{-6})
  tmp = optim(par = c(.5, 0), fun2opt, exc = exc.i, method = "Nelder")$par
  scale[i] = tmp[1]
  xi[i] = tmp[2]  
}
# now, xi, pexc and scale are vectors containing one estimate for each location in loc
