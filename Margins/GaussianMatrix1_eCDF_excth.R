#############################################################################################
##                                                                                         ##
## 1/4 codes to put the Gaussian data into a matrix, similar to the original anom.training ##
##                                                                                         ##
#############################################################################################
library(rje)
ncpus = 16
Replic = 16703
njobs = ceiling(Replic/ncpus)

anom.training.gauss = xi.hat = sigma.gp.hat = mu.gev.hat = sigma.gev.hat = iter.error = NULL

for(i in 1:(njobs/4)){ # i is location
  print(i)
  resto = Replic%%ncpus
  tmp = ncpus
  Rs <- (ncpus)*(i-1) + c(1:(tmp))
  ncpus = tmp
  tryload = tryCatch(load(paste0("~/EVAChallenge2019/MarginsThresh/OutputsByLocThx/Locs=", min(Rs), "-", max(Rs), ".Rdata")),
                     error = function(e) e)
  print(tryload)
  if(!inherits(tryload, 'error')){
    for(j in 1:ncpus){
      res = out[[j]]
      anom.training.gauss = cbind(anom.training.gauss, res$anom.training.gauss)
      
      xi.hat = rbind(xi.hat, res$xis)
      sigma.gp.hat = rbind(sigma.gp.hat, res$sigmas.gp)
      mu.gev.hat = rbind(mu.gev.hat, res$mus.gev)
      sigma.gev.hat = rbind(sigma.gev.hat, res$sigmas.gev)
      
    }
  }
  else
    iter.error = c(iter.error, i)
  
  save(anom.training.gauss, file = "~/EVAChallenge2019/MarginsThresh/OutputsAllTh/anom.training.gauss_excth_1.Rdata")
  save(iter.error, file = "~/EVAChallenge2019/MarginsThresh/OutputsAllTh/iter.error_excth_1.Rdata")
  save(xi.hat, file = "~/EVAChallenge2019/MarginsThresh/OutputsAllTh/xi.hat_excth_1.Rdata")
  save(sigma.gp.hat, file = "~/EVAChallenge2019/MarginsThresh/OutputsAllTh/sigma.gp.hat_excth_1.Rdata")
  save(sigma.gev.hat, file = "~/EVAChallenge2019/MarginsThresh/OutputsAllTh/sigma.gev.hat_excth_1.Rdata")
  save(mu.gev.hat, file = "~/EVAChallenge2019/MarginsThresh/OutputsAllTh/mu.gev.hat_excmu_1.Rdata")
  
}





