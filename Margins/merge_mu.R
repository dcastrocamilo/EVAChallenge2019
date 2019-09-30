##########################################################################################################
## Input data were divided in four parts, to get results faster. Here we merge all the transformed data ##
##########################################################################################################

# Part 1 ----
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/anom.training.gauss_excmu_1.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/iter.error_excmu_1.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/sigma.gp.hat_excmu_1.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/mu.gev.hat_excmu_1.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/xi.hat_excmu_1.Rdata")

anom.training.gauss1 = anom.training.gauss
iter.error1 = iter.error
xi.hat1 = xi.hat
sigma.gp.hat1 = sigma.gp.hat
mu.gev.hat1 = mu.gev.hat

rm(anom.training.gauss, iter.error, xi.hat, sigma.gp.hat, mu.gev.hat)

# Part 2 ----
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/anom.training.gauss_excmu_2.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/iter.error_excmu_2.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/sigma.gp.hat_excmu_2.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/mu.gev.hat_excmu_2.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/xi.hat_excmu_2.Rdata")

anom.training.gauss2 = anom.training.gauss
iter.error2 = iter.error
xi.hat2 = xi.hat
sigma.gp.hat2 = sigma.gp.hat
mu.gev.hat2 = mu.gev.hat

rm(anom.training.gauss, iter.error, xi.hat, sigma.gp.hat, mu.gev.hat)

# Part 3 ----
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/anom.training.gauss_excmu_3.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/iter.error_excmu_3.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/sigma.gp.hat_excmu_3.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/mu.gev.hat_excmu_3.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/xi.hat_excmu_3.Rdata")

anom.training.gauss3 = anom.training.gauss
iter.error3 = iter.error
xi.hat3 = xi.hat
sigma.gp.hat3 = sigma.gp.hat
mu.gev.hat3 = mu.gev.hat

rm(anom.training.gauss, iter.error, xi.hat, sigma.gp.hat, mu.gev.hat)

# Part 4 ----
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/anom.training.gauss_excmu_4.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/iter.error_excmu_4.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/sigma.gp.hat_excmu_4.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/mu.gev.hat_excmu_4.Rdata")
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/OutputsAllmu/xi.hat_excmu_4.Rdata")

anom.training.gauss4 = anom.training.gauss
iter.error4 = iter.error
xi.hat4 = xi.hat
sigma.gp.hat4 = sigma.gp.hat
mu.gev.hat4 = mu.gev.hat

rm(anom.training.gauss, iter.error, xi.hat, sigma.gp.hat, mu.gev.hat)

# Merge Gaussian-tranformed data ----
load("/Users/danielacastrocamilo/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/Data/DATA_TRAINING.RData")
dim(anom.training)

ncol(anom.training.gauss1)+
  ncol(anom.training.gauss2)+
  ncol(anom.training.gauss3)+
  ncol(anom.training.gauss4)

anom.training.gauss = cbind(anom.training.gauss1,
                             anom.training.gauss2,
                             anom.training.gauss3,
                             anom.training.gauss4)

dim(anom.training.gauss)

# # Check normality assumption ----
# i = ncol(anom.training)
# j = 2015
# k = 7
# for(i in 1:5){
#   for(j in unique(year)){
#     for(k in 1:12){
#       x = anom.training[month == k & year == j, i]
#       y = anom.training.gauss[month == k & year == j, i]
#       plot(density(y[!is.na(y)]))
#       qqnorm(y, pch = 16)
#       qqline(y, col = 2)
#     }
#   }
# }

set.seed(25072019)
sel = sample(1:ncol(anom.training), 25, replace = F)
pdf('~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/some_qqplots_mu.pdf',width = 10, height = 10)
par(mfrow = c(5,5), mar = c(2,2,2,2))
for(i in sel){
  x = anom.training[, i]
  y = anom.training.gauss[, i]
  qqnorm(y, pch = 16, ylab = '', xlab = '', main = '')
  qqline(y, col = 2)
  legend('topleft', paste0('St ', i), cex = 0.7)
}
dev.off()

# Check there are no errors ----
is.null(iter.error1)
is.null(iter.error2)
is.null(iter.error3)
is.null(iter.error4)

# Get GP, GEV parameters ----
dim(xi.hat1)
dim(xi.hat2)
dim(xi.hat3)
dim(xi.hat4)

dim(sigma.gp.hat1)
dim(sigma.gp.hat2)
dim(sigma.gp.hat3)
dim(sigma.gp.hat4)

dim(mu.gev.hat1)
dim(mu.gev.hat2)
dim(mu.gev.hat3)
dim(mu.gev.hat4)

xi.hat = cbind(xi.hat1, xi.hat2, xi.hat3, xi.hat4)
sigma.gp.hat = cbind(sigma.gp.hat1, sigma.gp.hat2, sigma.gp.hat3, sigma.gp.hat4)
mu.gev.hat = cbind(mu.gev.hat1, mu.gev.hat2, mu.gev.hat3, mu.gev.hat4)

# Add month and year to xi.hat, sigma.gp.hat, mu.gev.hat ----
load("~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/IBEXcluster/Data/DATA_TRAINING.RData")

df.xi = data.frame(month = rep(1:12, times = length(unique(year))), year = rep(unique(year), each = 12), xi.hat)
colnames(df.xi)[-(1:2)] = paste0('loc', 1:ncol(xi.hat))
xi.hat = df.xi
dim(df.xi)

df.sigma = data.frame(month = rep(1:12, times = length(unique(year))), year = rep(unique(year), each = 12), sigma.gp.hat)
colnames(df.sigma)[-(1:2)] = paste0('loc', 1:ncol(sigma.gp.hat))
sigma.gp.hat = df.sigma
dim(df.sigma)

df.mu = data.frame(month = rep(1:12, times = length(unique(year))), year = rep(unique(year), each = 12), mu.gev.hat)
colnames(df.mu)[-(1:2)] = paste0('loc', 1:ncol(mu.gev.hat))
mu.gev.hat = df.mu
dim(df.mu)

save(anom.training.gauss, file = '~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/anom.training.gauss_mu.Rdata')
save(xi.hat, file = '~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/xi.gp_mu.Rdata')
save(sigma.gp.hat, file = '~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/sigma.gp_mu.Rdata')
save(mu.gev.hat, file = '~/Dropbox/Projects/EVAChallenge2019/DanielaLindaThomas/Daniela/Margins/mu.gev_mu.Rdata')


