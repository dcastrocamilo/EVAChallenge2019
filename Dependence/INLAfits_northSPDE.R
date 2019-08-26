############################################################################################
##                                                                                        ##
## INLA fit for a neighborhood N(s,t) at the North of the Red sea using the SPDe approach ##
##                                                                                        ##  
############################################################################################

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

# SPDE ----
mesh = inla.mesh.2d(loc = sub.loc.degree, cutoff = 0.05, max.edge = c(0.3,10))
m = nrow(sub.data)
n = ncol(sub.data)
rep.id = rep(1:m, each = n)
x.loc = rep(sub.loc.degree[,1], m)
y.loc = rep(sub.loc.degree[,2], m)
A = inla.spde.make.A(mesh = mesh, loc = cbind(x.loc, y.loc), group = rep.id)

sdy = ceiling(sd(inla.data$y, na.rm = T))
dist = rdist.earth(sub.loc.degree, miles = F)
range.median = round(median(dist[lower.tri(dist)]), 1)
spde = inla.spde2.pcmatern(
  mesh= mesh, alpha = 2, 
  prior.range = c(range.median, 0.5), # P(range < range.median) = 0.5
  prior.sigma = c(2*sdy, 0.01)) # P(sigma > 2*sdy)= 0.01

mesh.index = inla.spde.make.index(name = "field", n.spde = spde$n.spde, n.group = m)

stack = inla.stack(tag = 'est', 
                   data = list(y = inla.data$y), 
                   A = list(A, 1), 
                   effects = list(mesh.index, 
                                  intercept = rep(1, length(inla.data$y))))

# INLA formula ----
form1.spde = y ~ -1 + intercept +
  f(field, model = spde, group = field.group, control.group = list(model = 'ar', order = 3))

# INLA fit 1 spde ----
fit1spde = inla(form1.spde,
                data = inla.stack.data(stack),
                control.predictor = list(link = 1, compute = TRUE, A = inla.stack.A(stack)),
                control.compute=list(config = TRUE),
                control.inla = list(strategy = "simplified.laplace",int.strategy = "eb"),
                verbose = T,
                num.threads = 2)
save(fit1spde, file = '~/EVAChallenge2019/Dependence/Out/fit1spde.Rdata')

############################
## INLA output fit 1 spde ##
############################
load('~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Dependence/Out/fit1spde.Rdata')
fit1spde$summary.fixed

fit1spde$summary.hyperpar
# Posterior samples (linear predictor, hyperparameters)
N = 100
samples.all = inla.posterior.sample(N, fit1spde)
index_pred = inla.stack.index(stack,"est")$data

# Generate samples from the data
psam = sapply(samples.all, function(x) {
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


