# Code file to run the local Gaussian dependence model and obtain posterior predictive samples of minimum over space-time cylinders.

library(fields)
library(INLA)

only.gaussian = FALSE # use Model 1 (no marginal transformations)
nthread = 3 #number of threads to use for running INLA 

thresh = 0.75 # threshold used for tail modeling of Z.tilde-process
# In the case of tail modeling (models 2 (GAM) and 3 (NN)), we assume that  vectors xi, scale, pexc, sd.by.year and mean.by.year have been generated precedingly (see file GAM_local_approaches.R).

# use PARDISO matrix computation library if installed 
# (not required but inla runs faster)
#inla.setOption(pardiso.license="/home/topitz/INLA/pardiso.lic")
#inla.pardiso.check()

# you may want to set a working directory:
#setwd(...)

iT = 1 # validation day (between 1 and 324)

source(paste0(WORK,"Tools-final.R"))

#number of INLA posterior samples to draw, and size of each sample:
nsamples = 1
N=500 

deltat=4 #temporal buffer: number of days before/after validation day (must be >= 3 for the 7-day hotspots to be predicted)

load('DATA_TRAINING.RData') #training data from EVA challenge
metric=readRDS('cartesianLoc.Rda') #loc in metric coordinates

if(!only.gaussian){
  # transform tails of data ####
  # generate anom.training.gauss with corrected tails
  # keep the initial data (Z.tilde process in the paper)
  means.vec = mean.by.year[year-min(year)+1]
  sds.vec = sd.by.year[year-min(year)+1]
  anom.training.gauss = (anom.training - means.vec)/sds.vec
  u.new = qnorm(1 - pexc)
  anom.training.Ztilde = anom.training.gauss
  for(i in 1:ncol(anom.training.gauss)){
    is.exc = which(anom.training.gauss[, i] > thresh)
    is.not.exc = which(anom.training.gauss[, i] <= thresh)
    anom.training.gauss[is.exc, i] = exc2normal(anom.training.Ztilde[is.exc, i] - thresh, pexc[i], scale[i], xi[i])
    n.notNA = sum(!is.na(anom.training.gauss[, i]))
    anom.training.gauss[is.not.exc, i] = qnorm(rank(anom.training.gauss.Ztilde[, i], na.last = "keep")/(n.notNA+1))[is.not.exc]
  }
}else{
  anom.training.gauss = anom.training 
}

#extract location index and time index of validation points (s_i,t_i):
nT=nrow(anom.training.gauss) #number of observation times
nS=ncol(anom.training.gauss) #number of observation sites
nvalid=length(index.validation) #number of validation points
#extract time indices of the validation points:
t.valid=index.validation %% nT 
t.valid[t.valid==0]=nT 
#store also the unique validation times:
t.valid.unique=unique(t.valid)
#length(unique(t.valid)) 324
#extract space indices of the validation points:
s.valid=index.validation %/% nT+1
s.valid[t.valid==nT]=s.valid[t.valid==nT]-1
#combine space and time indices:
st.valid=cbind(t.valid,s.valid) 
#save(st.valid,file=paste0("st.valid.RData"))


#fit model for the iT-th validation day ####
idx.time.sub = (t.valid.unique[iT]-deltat):(t.valid.unique[iT]+deltat)
idx.space.sub = 1:ncol(anom.training.gauss) 
data.sub = anom.training.gauss[idx.time.sub, idx.space.sub] 
dim(data.sub)
mean(is.na(data.sub))
loc.valid=loc[s.valid[t.valid==t.valid.unique[iT]],] 
metric.valid=metric[s.valid[t.valid==t.valid.unique[iT]],]
n.sub=prod(dim(data.sub))
nT.sub=dim(data.sub)[1]
nS.sub=dim(data.sub)[2]
metric.sub=metric[idx.space.sub,]

#prepare INLA ####
station = rep(1:ncol(data.sub), nrow(data.sub))
stationID = rep(idx.space.sub, nrow(data.sub))
months = month[idx.time.sub]
months = rep(months, each=ncol(data.sub))
years = year[idx.time.sub]
years = rep(years, each=ncol(data.sub))
days = day[idx.time.sub]
days = rep(days, each=ncol(data.sub))

if(!only.gaussian){
  # Add marginal tail parameters that will be needed to backtransform data from gaussian to original scale: ####
  xis = rep(xi, nrow(data.sub))
  scales = rep(scale, nrow(data.sub))
  pexcs = rep(pexc, nrow(data.sub))
  inla.data = data.frame(intercept = 1, 
                         y=c(t(data.sub)), 
                         station = station, 
                         stationID = stationID, 
                         day = days, 
                         month = months, 
                         year = years, 
                         xis = xis, 
                         scales = scales, 
                         pexcs = pexcs)
}else{
  inla.data = data.frame(intercept = 1,
                         y=c(t(data.sub)), 
                         station = station, 
                         stationID = stationID, 
                         day = days, 
                         month = months, 
                         year = years)
}

# remove big data files that we do not need any more
rm(anom.training.gauss); gc()

#set up SPDE model ####
bnd.int = inla.nonconvex.hull(metric.sub, convex=-.015)
bnd.ext = inla.nonconvex.hull(metric.sub, convex=-.4)
mesh = inla.mesh.2d(loc=metric.sub, cutoff=15, boundary=list(bnd.int, bnd.ext),max.edge = c(20,300),offset=c(-0.025, -.4))
#plot(mesh);mesh$n
#points(metric.sub,pch=19,cex=.15,col="blue")
rep.id = rep(1:nT.sub, each=nS.sub) 
A = inla.spde.make.A(mesh = mesh, loc=metric.sub, index=rep(1:nS.sub,nT.sub),group = rep.id)
#dim(A)

spde = inla.spde2.pcmatern(
  mesh= mesh, alpha = 2, 
  prior.range = c(500, .5), # P(range < 500) = 0.5
  prior.sigma = c(0.5, .5)) # (P(sd > 0.5) = 0.5)

mesh.index = inla.spde.make.index(name = "field", n.spde = spde$n.spde, n.group = nT.sub) 

stack = inla.stack(tag = "est", 
                   data = list(y = inla.data$y), 
                   A = list(A, 1), 
                   effects = list(mesh.index, intercept = rep(1, length(inla.data$y))))

hyper.ar1=list(theta=list(prior="pccor1",param=c(.85,.5))) # P(AR(1)-coefficient>0.85) = 0.5

# INLA formula ####
form.spde = y ~ -1 + intercept +
  f(field, 
    model = spde, 
    group = field.group, 
    control.group = list(model = "ar1", hyper=hyper.ar1))
  
# run INLA ####
fit = inla(form.spde,
                family="gaussian",
                quantiles=NULL,
                control.family=list(hyper=list(theta = list(prior="pc.prec", param=c(0.1, 0.5)))), # P(measurement error sd > 0.1) = 0.5
                data = inla.stack.data(stack),
                control.predictor = list(link = 1, compute = FALSE, A = inla.stack.A(stack)),
                control.compute = list(config = TRUE),
                control.inla = list(strategy = "simplified.laplace", int.strategy = "eb"),
                verbose = TRUE,
                num.threads = nthread)
summary(fit)

# save summary of fit ####
fithyper = fit$summary.hyperpar
fixed = fit$summary.fixed
runtime = fit$cpu.used
resid = inla.data$y - fit$summary.fitted.values$mean[1:length(inla.data$y)]
rmse = mean(resid^2, na.rm = TRUE)
save(fithyper, fixed, runtime, resid, rmse, iT, file = paste0("summary", iT, ".RData"))

# save the full fit (big file, around 1 GB) ####
#save(fit, mesh, stack, file = paste0("fit", iT, ".RData"))

#simulate from fitted model ####
# (to obtain the predictive distributions of minimum over N(s_i,t_i))

# Posterior samples (linear predictor, hyperparameters)
index_pred = inla.stack.index(stack, "est")$data

#get the validation points for the current validation date ####
idx2valid = which(t.valid == t.valid.unique[iT])
len = length(idx2valid)

if(len > 0){ #with EVA challenge data, len > 0 always
  idxS2valid = s.valid[idx2valid]
  idxT2valid = t.valid[idx2valid[1]]
  obs.true = obs2pred[idx2valid]
  samples.min = vector(mode = "list", length = len)
  
  for(id.sample in 1:nsamples){
    set.seed(id.sample)
    samples.all = inla.posterior.sample(N, fit, num.threads=nthread)
    # save full sample if you want to (but generates big file) ####
    #save(samples.all, file = paste0("samples", iT, "-nsample", id.sample, ".RData"))
    for(m in 1:len){
      #cat("m is ", m, "\n")
      samples.min[[m]] = c(samples.min[[m]],
                           drop(get.min.sample(samples.all,
                                               loc,
                                               id.station = idxS2valid[m],
                                               idx.space.sub,
                                               deltat = deltat,
                                               nS.sub = nS.sub,
                                               doBackT = !only.gaussian)))
    }
  }

  for(m in 1:len){
    # get (empirical) cdf values of predictive distribution for space-time minimum
      sample.pred = samples.min[[m]]
      ecdf.m = ecdf(samples.min[[m]])
      pred.m = ecdf.m(-1+1:400/100)
      save(pred.m, m, iT, idx2valid, sample.pred, file = paste0("pred-", idx2valid[m], ".RData"))
  }
}
