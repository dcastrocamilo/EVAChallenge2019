#############################################################################################################################################
##                                                                                                                                         ##
## This code, suited to run in a Cluster with mc.cores cores, transform the data into Gaussian scale using the GEV mean to define extremes ##
##                                                                                                                                         ##
#############################################################################################################################################
#!/usr/bin/Rscript
args <- commandArgs(TRUE)
for (arg in args) eval(parse(text = arg))
rm(arg, args)

library(parallel)
library(mgcv)

print('Loading data')
load("~/EVAChallenge2019/IBEXcluster/Data/DATA_TRAINING.RData")
load("~/EVAChallenge2019/IBEXcluster/Data/dist2coast.Rdata")
load("~/EVAChallenge2019/IBEXcluster/Data/mod_xi_time_lat.Rdata")
load("~/EVAChallenge2019/IBEXcluster/Data/prob.fit.Rdata")
load("~/EVAChallenge2019/IBEXcluster/Data/quant.mont.95.Rdata")
loc = data.frame(loc)

## -------------------------------------------- ##
## Transformation to exponential scale          ##
Tx <- function(x, xi, mu, sigma){               ##        
  if(is.na(x))                                  ##
    return(NA)                                  ##
  else{                                         ##
    if(xi == 0)                                 ##
      return((x - mu)/sigma)                    ##
    else                                        ##
      return((1/xi)*log(1 + xi*(x - mu)/sigma)) ##
  }                                             ##
                                                ##
}                                               ##
## Transformation to Gaussian scale             ##
Sx <- function(x){ # x is an exponential r.v.   ##
  if(is.na(x))                                  ##
    return(NA)                                  ##
  else{                                         ##
    u <- pexp(x)                                ##
    return(qnorm(u))                            ##
  }                                             ##
}                                               ##
## -------------------------------------------- ##


# This function transform the original data to Gaussian scale for a given location i ----
get.data.i <- function(i){
  lat.i = loc$lat[i]
  dist.i = dist2coast$distance[i]
  
  w = ecdf(anom.training[, i])(anom.training[, i])
  
  anom.training.exp = anom.training.unif = anom.training.gauss = NULL
  xis = sigmas.gp = mus.gev = sigmas.gev = NULL
  
  for(j in unique(year)){
    
    for(k in 1:12){
      # Shape GP
      xi <- as.numeric(predict(mod_xi_time_lat$xiObj, newdata = data.frame("month" = k, "year" = j, "lat" = lat.i, "dist" = dist.i)))
      xis = c(xis, xi)
      #print(xi)
      # Scale GP
      sigma.gp <- as.numeric(exp(predict(mod_xi_time_lat$nuObj, newdata = data.frame("month"= k,"year" = j, "lat" = lat.i, "dist" = dist.i)))/
                               (1 + predict(mod_xi_time_lat$xiObj, newdata = data.frame("month" = k, "year" = j, "lat" = lat.i, "dist" = dist.i))))
      sigmas.gp = c(sigmas.gp, sigma.gp)
      #print(sigma.gp)
      # Exceedance probability
      pr = as.numeric(prob.fit$prob[prob.fit$month == k & prob.fit$year == j & prob.fit$dist == dist.i & prob.fit$lat == lat.i][1])
      # Threshold
      thresh = quant.mont.95[k, i]
      # Location and Scale GEV
      if(xi == 0){
        mu.gev = thresh + sigma.gp*log(pr)
        sigma.gev = sigma.gp
      }else{
        mu.gev = thresh - sigma.gp*(pr^(-xi) - 1)/(xi*pr^(-xi))
        sigma.gev = sigma.gp - sigma.gp*(pr^(-xi) - 1)/pr^(-xi)
      }
      mus.gev = c(mus.gev, mu.gev)
      sigmas.gev = c(sigmas.gev, sigma.gev)
      # Observations in original scale (all days within month {k} of year {j})
      x = anom.training[month == k & year == j, i]
      # Observations in uniform scale (all days within month {k} of year {j}). eCDF computed for each month {k} at each site {i}
      # w.x = ecdf(anom.training[month == k, i])(anom.training[month == k & year == j, i])
      # Exceedance w.r.t. location mu
      is.exc = x > mu.gev 
      ###################
      # For exceedances #
      ###################
      y = u = z = numeric(length(x))
      if(sum(is.exc, na.rm = T) > 0){
        x.exc = x[which(is.exc)]
        # Observations in exponential scale
        y[which(is.exc)] <- sapply(x.exc, Tx, xi = xi, mu = mu.gev, sigma = sigma.gev)
        # Observations in Gaussian scale
        z[which(is.exc)] <- sapply(y[which(is.exc)], Sx)
      }
      #######################
      # For non-exceedances #
      #######################
      if(sum(!is.exc, na.rm = T) > 0){
        x.nonexc = x[which(!is.exc)]
        w.x = w[month == k & year == j]
        # Observations in uniform scale
        y[which(!is.exc)] = w.x[which(!is.exc)]
        # Observations in Gaussian scale
        z[which(!is.exc)] <- sapply(y[which(!is.exc)], qnorm)
      }
      ###########
      # For NAs #
      ###########
      if(sum(is.na(x)) > 0){
        x.na = which(is.na(x))
        y[is.na(x)] = u[is.na(x)] = z[is.na(x)] = NA
      }
      
      anom.training.gauss = c(anom.training.gauss, z)
    }
  }
  
  return(list(anom.training.gauss = anom.training.gauss, xis = xis, sigmas.gp = sigmas.gp, mus.gev = mus.gev, sigmas.gev = sigmas.gev))
}

print('Running code')
out <- mclapply(Rs, get.data.i, mc.cores = mc.cores)
save(out, file = paste0("~/EVAChallenge2019/Margins/OutputsByLocMu/out_Rs=", min(Rs), "-", max(Rs), ".Rdata"))
print('Outputs saved')


