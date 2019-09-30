load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/mu.gev_mu.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/sigma.gp_mu.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/xi.gp_mu.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/Daniela/Margins/anom.training.gauss_mu.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019_myversion/DanielaLindaThomas/DATA_TRAINING.RData")

# library(devtools)
# devtools::install_github('talgalili/edfun')
library("edfun")

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
z = anom.training.gauss[, i] # Obs in Gaussian scale
x = anom.training[, i] # Obs in original scale
xedfun = edfun(x[!is.na(x)])
ecdf.inv = xedfun$qfun # this is the inverve function

mu.x = mu.gev.hat[, (i+2)]
xi.x = xi.hat[, (i+2)]
sigma.x = sigma.gp.hat[, (i+2)]
y = rep(NA, length(x)) # inverse of anom.training.gauss_mu
y.nexc = y.exc = NULL # here I save the backtramsformed data for non-exceedances and exceedances, respectively.
x.nexc = x.exc = NULL # here I save the original non-exceedances and exceedances data, respectively.
for(j in 1:length(x)){
  mu.j = mu.x[mu.gev.hat$month == month[j] & mu.gev.hat$year == year[j]]
  xi.j = xi.x[xi.hat$month == month[j] & xi.hat$year == year[j]]
  sigma.j = sigma.x[sigma.gp.hat$month == month[j] & sigma.gp.hat$year == year[j]]
  
  if(!is.na(x[j])){
    if(x[j] > mu.j){
      y[j] = Tx.inv(qexp(pnorm(z[j])), xi = xi.j, mu = mu.j, sigma = sigma.j)
      y.nexc = c(y.nexc, y[j])
      x.nexc = c(x.nexc, x[j])
    }else{
      y[j] = ecdf.inv(pnorm(z[j]))
      y.exc = c(y.exc, y[j])
      x.exc = c(x.exc, x[j])
    }
  }
  
}
plot(x, y)
abline(0, 1, col = 2, led = 2)








