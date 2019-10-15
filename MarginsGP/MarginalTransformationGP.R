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

# Data ----
print('Loading data')
load(paste0(getwd(), "/Data/DATA_TRAINING.RData"))
load(paste0(getwd(), "/Data/dist2coast.Rdata"))
load(paste0(getwd(), "/Data/mod_xi_time_lat.Rdata"))
loc = data.frame(loc)

# GEV GP parameters ----
load(file = paste0(getwd(), "/Data/param.Rdata"))
load(paste0(getwd(), "/Data/prob.fit.Rdata"))
load(paste0(getwd(), "/Data/quant.mont.95.Rdata"))


## Transformation from "raw" data to uniform scale
Tx <- function(x, xi, thresh, sigma, pr){             
  if(is.na(x))                                 
    return(NA)                                 
  else{                                        
    1-(pr*(1-pgpd(x, thresh, sigma, xi)))      
  }                                            
}

# i = 1; j = 1994; k = 1
# This function transform the original data to Gaussian scale for a given location i ----
get.data.i <- function(i, months = NULL, years = NULL){
  lat.i = loc$lat[i]
  dist.i = dist2coast$distance[i]
  
  x = anom.training[, i]
  
  anom.training.gauss = anom.training.unif = NULL
  
  if(is.null(years)) years = unique(year)
  if(is.null(months)) months = 1:12
  
  for(j in years){
    
    for(k in months){
      # Observations in original scale for month {k} and year {j}
      x.kj = x[month == k & year == j]
    
      which.is = which(param$month == k & param$year == j & param$location == i)
      xi = param$xi[which.is]
      sigmagp = param$sigma.gp[which.is]
      # sigmagev = param$sigma.gev[which.is]
      # mugev = param$mu.gev[which.is]
      
      pr = as.numeric(prob.fit$prob[prob.fit$month == k & prob.fit$year == j & prob.fit$dist == dist.i & prob.fit$lat == lat.i][1])
      thresh = quant.mont.95[k, i]
      
      if(all(is.na(x.kj))){
        z = u = x.kj
        
      }else{
        # Observations in uniform and Gaussian scale for month {k} and year {j}
        u = edfun(x.kj[!is.na(x.kj)])$pfun(x.kj)*(length(x.kj)/(length(x.kj) + 1))
        z = qnorm(u)
        
        is.exc = x.kj > thresh
        
        if(sum(is.exc, na.rm = T) > 0){
          x.exc = x.kj[is.exc] # Exceedances for month {k} and year {j}
          u.exc <- sapply(x.exc, Tx, xi = xi, thresh = thresh, sigma = sigmagp, pr = pr)
          u[is.exc] <- u.exc
          z[is.exc] <- sapply(u.exc, qnorm) # replacing z values in the presence of exceedances
        }
        
      }
      anom.training.gauss = c(anom.training.gauss, z)
      anom.training.unif = c(anom.training.unif, u)
      
      # qqnorm(z)
      # qqline(z)
      # hist(u, freq = F)
      
    }
  }
  return(list(anom.training.gauss = anom.training.gauss, anom.training.unif = anom.training.unif))
}

print('Running code')
out <- mclapply(Rs, get.data.i, mc.cores = mc.cores)
save(out, file = paste0("~/EVAChallenge2019/MarginsThresh/OutputsByLocTh/Locs=", min(Rs), "-", max(Rs), ".Rdata"))
print('Outputs saved')


