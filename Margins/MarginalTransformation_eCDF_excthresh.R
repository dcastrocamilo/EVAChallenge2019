##################################################################################################################################################
##                                                                                                                                              ##
## This code, suited to run in a Cluster with mc.cores cores, transform the data into Gaussian scale using the GPD threshold to define extremes ##
##                                                                                                                                              ##
##################################################################################################################################################
#!/usr/bin/Rscript
args <- commandArgs(TRUE)
for (arg in args) eval(parse(text = arg))
rm(arg, args)

library(parallel)
library(mgcv)
library(edfun)

# rm(list=setdiff(ls(), c("anom.training", "dist2coast", "mod_xi_time_lat", "prob.fit", "quant.mont.95", "loc", "year", "month", "time")))

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
## -------------------------------------------- ##

# i = 1; j = 1994; k = 1

# This function transform the original data to Gaussian scale for a given location i ----
get.data.i <- function(i, years = NULL, months = NULL){
  lat.i = loc$lat[i]
  dist.i = dist2coast$distance[i]
  
  x = anom.training[, i]
  # w = edfun(x[!is.na(x)])$pfun(x)*(length(x)/(length(x) + 1))
  
  anom.training.gauss = anom.training.unif = NULL
  # unif.nexc = unif.exc = gauss.nexc = gauss.exc = NULL # for testing
  xis = sigmas.gp = mus.gev = sigmas.gev = NULL
  
  if(is.null(years)) years = unique(year)
  if(is.null(months)) months = 1:12
  
  for(j in years){
    
    for(k in months){
      # Observations in original scale for month {k} and year {j}
      x.kj = x[month == k & year == j]
      
      # Shape GP
      xi <- as.numeric(predict(mod_xi_time_lat$xiObj, newdata = data.frame("month" = k, "year" = j, "lat" = lat.i, "dist" = dist.i)))
      xis = rbind(xis, c(k, j, i, xi))
      # xis = c(xis, xi)
      
      # Scale GP
      sigma.gp <- as.numeric(exp(predict(mod_xi_time_lat$nuObj, newdata = data.frame("month"= k,"year" = j, "lat" = lat.i, "dist" = dist.i)))/
                               (1 + predict(mod_xi_time_lat$xiObj, newdata = data.frame("month" = k, "year" = j, "lat" = lat.i, "dist" = dist.i))))
      sigmas.gp = rbind(sigmas.gp, c(k, j, i, sigma.gp))
      # sigmas.gp = c(sigmas.gp, sigma.gp)
      
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
      mus.gev = rbind(mus.gev, c(k, j, i, mu.gev))
      sigmas.gev = rbind(sigmas.gev, c(k, j, i, sigma.gev))
      # mus.gev = c(mus.gev, mu.gev)
      # sigmas.gev = c(sigmas.gev, sigma.gev)
      
      if(all(is.na(x.kj))){
        z = u = x.kj
        
      }else{
        # Observations in uniform scale for month {k} and year {j}
        w.kj = edfun(x.kj[!is.na(x.kj)])$pfun(x.kj)*(length(x.kj)/(length(x.kj) + 1))
        # w.kj = w[month == k & year == j]
        
        # Exceedance w.r.t. threshold thresh
        is.exc = x.kj > thresh
        
        # Transformation for non-exceedances
        z = qnorm(w.kj) # Gaussian scale for exceedances and non-exceedances
        u = w.kj
        # u.nexc = w.kj[!is.exc] # uniform scale for non-exceedances # for testing
        # z.nexc = qnorm(u.nexc) # Gaussian scale for non-exceedances # for testing
        
        
        # Transformation for exceedances
        if(sum(is.exc, na.rm = T) > 0){
          x.exc = x.kj[is.exc] # Exceedances for month {k} and year {j}
          tmpE <- sapply(x.exc, Tx, xi = xi, mu = mu.gev, sigma = sigma.gev) # Observations in exponential scale
          u.exc <- sapply(tmpE, pexp) 
          u[is.exc] <- u.exc
          # z.exc <- sapply(u.exc, qnorm) # Gaussian scale for non-exceedances # for testing
          z[is.exc] <- sapply(u.exc, qnorm) # replacing z values in the presence of exceedances
        }
        
      }
      anom.training.gauss = c(anom.training.gauss, z)
      anom.training.unif = c(anom.training.unif, u)
      
      # qqnorm(z)
      # qqline(z)
      # hist(u, freq = F)
      # unif.nexc = c(unif.nexc, u.nexc)
      # unif.exc = c(unif.exc, u.exc)
      # gauss.nexc = c(gauss.nexc, z.nexc)
      # gauss.exc = c(gauss.exc, z.exc)
      
    }
  }
  
  xis = data.frame(xis); colnames(xis) = c('month', 'year', 'location', 'xi')
  sigmas.gp = data.frame(sigmas.gp); colnames(sigmas.gp) = c('month', 'year', 'location', 'sigma.gp')
  mus.gev = data.frame(mus.gev); colnames(mus.gev) = c('month', 'year', 'location', 'mu.gev')
  sigmas.gev = data.frame(sigmas.gev); colnames(sigmas.gev) = c('month', 'year', 'location', 'sigma.gev')
  
  
  return(list(anom.training.gauss = anom.training.gauss, anom.training.unif = anom.training.unif, 
              xis = xis, sigmas.gp = sigmas.gp, mus.gev = mus.gev, sigmas.gev = sigmas.gev))
              # unif.nexc = unif.nexc, unif.exc = unif.exc, gauss.nexc = gauss.nexc, gauss.exc = gauss.exc, w = w))
}

print('Running code')
out <- mclapply(Rs, get.data.i, mc.cores = mc.cores)
save(out, file = paste0("~/EVAChallenge2019/MarginsThresh/OutputsByLocTh/Locs=", min(Rs), "-", max(Rs), ".Rdata"))
print('Outputs saved')


