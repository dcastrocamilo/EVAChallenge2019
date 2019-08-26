###############################################################################################################
##                                                                                                           ##
## INLA fit for a neighborhood N(s,t) at the North of the Red sea using the generic0 model and a AR(3) model ##
##                                                                                                           ##  
###############################################################################################################

library(fields)
library(INLA)
library(rje)
# MAC
# source('~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Dependence/Tools.R')
# load('~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/anom.training.gauss_mu.Rdata')
# load('~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/sigma.gp_mu.Rdata')
# load('~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/mu.gev_mu.Rdata')
# load('~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/xi.gp_mu.Rdata')
# load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/Data/DATA_TRAINING.RData")
# metric = readRDS('~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/cartesianLoc.Rda')
# LINUX
source('~/EVAChallenge2019/Dependence/Tools.R')
load('~/EVAChallenge2019/Margins/sigma.gp_mu.Rdata')
load('~/EVAChallenge2019/Margins/mu.gev_mu.Rdata')
load('~/EVAChallenge2019/Margins/xi.gp_mu.Rdata')
load('~/EVAChallenge2019/Dependence/anom.training.gauss_mu.Rdata')
load('~/EVAChallenge2019/Dependence/DATA_TRAINING.RData')
metric = readRDS('~/EVAChallenge2019/Dependence/cartesianLoc.Rda')

############################################
## Three locations (north, centre, south) ##
############################################
loc.id = c(1000, 8500, 15800)
# plot(loc)
# points(loc[1:4000,], col = 3, lwd = 5)
# points(loc[loc.id, ], col = 2, lwd = 5, pch = 16)

# plot(metric)
# points(metric[loc.id, ], col = 2, lwd = 5, pch = 16)

##################################################
## Cilynder of 50km radius around each location ##
##################################################
r = 50
tmp = rdist.earth(loc[loc.id, ], loc, miles = F)
which.neigh.north = which(tmp[1, ] < 50)
# which.neigh.centre = which(tmp[2, ] < 50)
# which.neigh.south = which(tmp[3, ] < 50)

##############################
##                          ##
## INLA fit North in N(s,t) ##
##                          ##
##############################
sub.loc.metric = metric[which.neigh.north, ]
sub.loc.degree = loc[which.neigh.north, ]
# Subset data 1: N(s,t) as defined in the challenge description ----
sub.data = anom.training.gauss[ ,which.neigh.north]
time.period = day[11000:11006]
sub.data = sub.data[time.period, ] 

# Subset data 2 ----
# sub.data = anom.training.gauss[ ,which.neigh.north]
# time.period = year >= 2010 & year <= 2015
# sub.data = sub.data[time.period, ]

# Model components ----
station = rep(1:ncol(sub.data), each = nrow(sub.data))
stationID = rep(which.neigh.north, each = nrow(sub.data))
months = month[time.period]
months = rep(months, each = ncol(sub.data))
years = year[time.period]
years = rep(years, each = ncol(sub.data))
days = day[time.period]
days = rep(days, each = ncol(sub.data))

# Extra components to back-transform the data ----
xi.h = xi.hat[, -(1:2)]
sigma.h = sigma.gp.hat[, -(1:2)]
mu.h = mu.gev.hat[, -(1:2)]
xihat = sigmahat = muhat = rep(NA, length(months))
for(i in 1:length(months)){
  printPercentage(i, length(months))
  xihat[i] = xi.h[xi.hat$month == months[i] & xi.hat$year == years[i], stationID[i]]
  sigmahat[i] = sigma.h[sigma.gp.hat$month == months[i] & sigma.gp.hat$year == years[i], stationID[i]]
  muhat[i] = mu.h[mu.gev.hat$month == months[i] & mu.gev.hat$year == years[i], stationID[i]]
}

inla.data = data.frame(intercept = 1, y = c(sub.data), station = station, stationID = stationID, day = days, month = months, year = years, mu.hat = muhat, xi.hat = xihat, sigma.hat = sigmahat)

# Introduce missing values ----
INLA.DATA = inla.data
my.missing = (1:9)*5
inla.data$y[my.missing] = NA

# Correlation structure ----
matrange = 40 # Fixed range and smoothness param
nu = 1
kappa = sqrt(2*nu)/matrange
dist = as.matrix(dist(metric[which.neigh.north, ]))
cormat = as.matrix(2^(1-nu)*(kappa*dist)^nu*besselK(dist*kappa,nu)/gamma(nu)) # Matérn correlation matrix
diag(cormat) = 1
prec = solve(cormat) # Compute precision matrix

# INLA formula ----
form1 = y ~ -1 + intercept + f(station, model = "generic0", Cmatrix = prec, constr = TRUE) + 
  f(day, model="ar", order = 3)

# form2 = y ~ -1 + intercept + f(station, model = "generic0", Cmatrix = prec, constr = TRUE) + 
#   f(month, model = "rw2", cyclic = TRUE, hyper=list(prec=list(initial=log(1/.01^2), fixed = TRUE)), constr = TRUE) +
#   f(day, model="ar", order = 3)

# INLA fit ----
fit1 = inla(form1,
            data = inla.data,
            control.predictor = list(link = 1),
            control.compute=list(config = TRUE),
            control.inla = list(strategy = "simplified.laplace",int.strategy = "eb"),
            verbose = T,
            num.threads = 2)
save(fit1, file = '~/EVAChallenge2019/Dependence/Out/fit1.Rdata')

#######################
## INLA output fit 1 ##
#######################
load('~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Dependence/Out/fit1.Rdata')
fit1$summary.fixed

# Spatial effect
sub.loc.degree[which(which.neigh.north == loc.id[1]), ]
quilt.plot(x = sub.loc.degree[,1], y = sub.loc.degree[,2], z = fit1$summary.random$station$mean)

# AR effect
ry = c(min(fit1$summary.random$day$`0.025quant`), max(fit1$summary.random$day$`0.975quant`))
par(mar = c(4,4,2,2))
plot(time.period, fit1$summary.random$day$mean, ylim = ry, type = 'n')
polygon(x = c(rev(time.period), time.period), y = c(rev(fit1$summary.random$day$`0.025quant`), fit1$summary.random$day$`0.975quant`), border = NA, col = 'lightgray')
lines(time.period, fit1$summary.random$day$mean, type = 'b', pch = 16)

# Posterior samples (linear predictor, hyperparameters)
N = 100
samples.all = inla.posterior.sample(N, fit1)

# Generate samples from the data
psam = sapply(samples.all, function(x) {
  index_pred = startsWith(rownames(x$latent), 'Predictor')
  mu = x$latent[index_pred, 1]
  s.d = sqrt(1/x$hyperpar[1])
  err = rnorm(length(mu), 0, s.d)
  u =  mu + err
  return(u) 
})

# Back-transform samples to original scale
psamp.orig.scale = matrix(NA, nrow(psam), ncol(psam))
for(i in 1:nrow(psam))
psamp.orig.scale[i, ] = sapply(psam[i, ], backT, xi = inla.data$xi.hat[i], mu = inla.data$mu.hat[i], sigma = inla.data$sigma.hat[i])

# Comparison between posterior samples and my missing values
par(mfrow = c(3,3), mar = c(2,2,2,2))
for(i in 1:length(my.missing)){
  plot(density(psamp.orig.scale[my.missing[i], ]), main = '')
  abline(v = INLA.DATA$y[my.missing[i]], col = 2)
}

# Compute X
X.min = apply(psamp.orig.scale, 2, min) # posterior sample of X in N(s,t)


