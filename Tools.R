
#function to generate data samples based on the posterior model ####
sim.from.posterior=function(x,idx2keep) {
  mu = x$latent[idx2keep, 1]
  s.d = sqrt(1/x$hyperpar[1])
  err = rnorm(length(mu), 0, s.d)
  u =  mu + err
  return(u) 
}

# function to get minimum over N(s_i,t_i) (radius 50km, 3-day buffer before and after) from inla posterior sample ####
# samples.inla : the output of inla.posterior.sample
# loc : pixel coordinate matrix (2 columns, latitude-longitude)
# idx.space.sub : pixel indices which include the neighborhood N(s_i,t_i)
# id.station : index of station in the vector of pixel locations 
# deltat : temporal buffer (number of days before/after validation date) used for fitting the model
#nS.sub : number of spatial locations in INLA model
get.min.sample = function(samples.inla, loc, id.station, idx.space.sub, deltat, nS.sub, doBackT = TRUE){
  idx.spat.pred=which(rdist.earth(loc[id.station,,drop=FALSE],loc[idx.space.sub,],miles=FALSE)<50) # get indices of locations in 50km radius
  idx2extract = c()
  for(i in 1:7){
    idx2extract=c(idx2extract, (deltat-3+i-1)*nS.sub + idx.spat.pred)
  }
  # generate matrix N times number of ST sample points
  sample = sapply(samples.inla, sim.from.posterior, idx2keep = idx2extract)
  if(doBackT){
    sample = mybackT(sample, xis = inla.data$xis[idx2extract], scales = inla.data$scales[idx2extract], pexcs=inla.data$pexcs[idx2extract], years = inla.data$year[idx2extract], stations = inla.data$stationID[idx2extract])
  }
  apply(sample,2,min)
}

# function to transform simulations of local models on Gaussian scale back to original scale ####
mybackT = function(pred.norm, xis, scales, pexcs, years, stations){ 
  #rows of pred.norm correspond to location-time
  #columns of pred.norm correspond to replicated samples
  pred.unif = pnorm(pred.norm) # predictions in uniform scale
  pred.orig = matrix(NA, nrow = nrow(pred.norm), ncol = ncol(pred.norm))
  #first, transform back using ranks:
  for(i in unique(stations)){
    pred.orig[stations == i,] = quantile(anom.training.gauss.Ztilde[,i], probs = pred.unif[stations == i,], na.rm = TRUE)
  }
  #then, replace tails by GPD back-transformed simulations:
  for(i in 1:ncol(pred.norm)){
    idx.above = which(pred.unif[,i] > 1 - pexcs)
    if(length(idx.above)>0){
      pred.orig[idx.above,i]= thresh + 
        qgp((pnorm(pred.norm[idx.above,i])-(1-pexcs[idx.above]))/pexcs[idx.above], scales[idx.above], xis[idx.above])
    }
    pred.orig[,i] = pred.orig[,i]*sd.by.year[years - 1985 + 1] + mean.by.year[years - 1985 + 1]
  }
  pred.orig
}

#density of generalized Pareto distribution ####
dgp=function(x,sigma,xi){
  if(abs(xi)<10^{-4}){
    exp(-x/sigma)/sigma
  }else{
    ifelse(1+xi*x/sigma<=0,0,sigma^{-1}*(1+xi*x/sigma)^{-1/xi-1})
  }
}

# cdf of generalized Pareto distribution ####
pgp = function(x, sigma, xi){
  if(abs(xi)<10^{-4}){
    1-exp(-x/sigma)
  }else{
    1- ifelse(1+xi*x/sigma<=0, 0, (1+xi*x/sigma)^{-1/xi})
  }
}

# quantile function of generalized Pareto distribution ####
qgp = function(x, sigma, xi){
  if(abs(xi)<10^{-4}){
    -sigma * log(1-x)
  }else{
    sigma/xi * ((1-x)^{-xi}-1)
  }
}

# negative log-likelihood of generalized Pareto distribution ####
nll.gp=function(par, exc){
  if(par[1] <= 0) return(Inf)
  -sum(log(dgp(x=exc, sigma = par[1], xi = par[2])))
}

# function to transform exceedances in data to standard normal scale ####
exc2normal = function(exc, pexc, scale, xi){
  qnorm((1-pexc) + pexc * pgp(exc, scale, xi))
}

# twCRPS function ####
# INPUTS:
# prediction: matrix (dimension 162000x400, size 494.4 Mb). The i-th row contains the predicted distribution for the i-th point in the validation set, evaluated at each of the 400 design points
# true.observations: vector (length 162000, size 1.2 Mb). Observation vector (unknown to the teams), containing the true spatio-temporal minimum X(s,t)=min A(s,t) for each point in the validation set
twCRPS <- function(prediction, true.observations){
  n.validation <- length(true.observations) # number of observations in the validation set (here, equal to 162000)
  xk <- -1+c(1:400)/100 # 'design points' used to evaluate the predicted distribution
  n.xk <- length(xk) # number of design points (here, equal to 400)
  weight <- function(x){ # Define weight function for computing the threshold-weighted CRPS
    return( pnorm((x-1.5)/0.4) )
  }
  wxk <- weight(xk) # weights for each of the 'design points'
  
  true.observations.mat <- matrix(true.observations,nrow=n.validation,ncol=n.xk) # true observations put in matrix form
  xk.mat <- matrix(xk,nrow=n.validation,ncol=n.xk,byrow=TRUE) # design points put in matrix form
  wxk.mat <- matrix(wxk,nrow=n.validation,ncol=n.xk,byrow=TRUE) # weights put in matrix form
  
  twCRPS.res <- 4*mean((prediction-(true.observations.mat<=xk.mat))^2*wxk.mat)
  return(twCRPS.res)
}
