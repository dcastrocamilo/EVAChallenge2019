######################################################################################
##                                                                                  ##
## Back-transform Gaussian data and compare to original values (should be the same) ##
##                                                                                  ##
######################################################################################
# load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/mu.gev_mu.Rdata")
# load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/sigma.gp_mu.Rdata")
# load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/xi.gp_mu.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/IBEXcluster/Data/quant.mont.95.Rdata")
load("~//Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/IBEXcluster/Data/dist2coast.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/IBEXcluster/Data/mod_xi_time_lat.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/IBEXcluster/Data/prob.fit.Rdata")
# load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/anom.training.gauss_mu.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/DATA_TRAINING.RData")
load("~/Dropbox/Projects/EVAChallenge2019/Margins/Checks/out_Rs=1-1_excthres.Rdata")

# library(devtools)
# devtools::install_github('talgalili/edfun')
library(edfun)
library(rje)

Tx.inv <- function(y, xi, mu, sigma){
  if(is.na(y))            
    return(NA)            
  else{                   
    if(xi == 0)           
      return(sigma*y + mu)
    else                  
      return((exp(xi*y)-1)*sigma/xi + mu)
  }
}

i = 1 # location
x = anom.training[, i] # Observations in original scale
z = out[[1]]$anom.training.gauss # Observations in Gaussian scale

qqnorm(z)
qqline(z, col = 2)

u = pnorm(z) # Observations in Uniform scale

# GP and GEV parameters
xis = out[[1]]$xis
sigmas.gp = out[[1]]$sigmas.gp
sigmas.gev = out[[1]]$sigmas.gev
mus.gev = out[[1]]$mus.gev

x.tilde = exc = rep(NA, length(x))

for(j in unique(year)){
  
  for(k in 1:12){
    # Observations in original scale for month {k} and year {j}
    x.kj = x[month == k & year == j]
    
    if(!all(is.na(x.kj))){
      # Shape GP
      xi = xis[xis$month == k, xis$year == j, xis$location == i, ]
      # Scale GP
      sigma.gp = sigmas.gp[sigmas.gp$month == k, sigmas.gp$year == j, sigmas.gp$location == i, ]
      # Location GEV
      mus.gp = mus.gev[mus.gev$month == k, mus.gev$year == j, mus.gev$location == i, ]
      # Scale GEV
      sigmas.gev = sigmas.gev[sigmas.gev$month == k, sigmas.gev$year == j, sigmas.gev$location == i, ]
      # Exceedance probability
      pr = as.numeric(prob.fit$prob[prob.fit$month == k & prob.fit$year == j & prob.fit$dist == dist.i & prob.fit$lat == lat.i][1])
      # Threshold
      thresh = quant.mont.95[k, i]
      
      # Observations in Uniform scale for month {k} and year {j}
      u.kj = u[month == k & year == j]
      
      is.exc = x.kj > thresh
      
      # Back-transformation for non-exceedances
      tmp = edfun(x.kj[!is.na(x.kj)])$qfun(x.kj)*(length(x.kj)/(length(x.kj) + 1))
      
      # Back-transformation for exceedances
      if(sum(is.exc) > 0)
        tmp[is.exc] = sapply(qexp(u.kj[is.exc]), Tx.inv, xi = xi, mu = mu, sigma = sigma)
      
      plot(x.kj, tmp)
      
      x.tilde[month == k & year == j] = tmp
      exc[month == k & year == j] = is.exc
    }
    
  }
}

png('~/Dropbox/Projects/EVAChallenge2019/Margins/Checks/Backtransform_excthres_loc_1_all.png', width = 960, height = 960)
plot(x, x.tilde)
abline(0, 1, col = 2, lwd = 2)
dev.off()

png('~/Dropbox/Projects/EVAChallenge2019/Margins/Checks/Backtransform_excthres_loc_1_exc&noexc.png', width = 960, height = 960)
par(mfrow = c(1,2), mar = c(1,4,5,1), pty = 's')
plot(x.exc, y.exc, main = 'x that exceed the threshold')
abline(0, 1, col = 2, lwd = 2)

plot(x.nexc, y.nexc, main = 'x that exceed do not the threshold')
abline(0, 1, col = 2, lwd = 2)
dev.off()







