####################################################################################################################################################
##                                                                                                                                                ##                     
## This code allows us to identify the pattern of missing values: For every month of every year, either we have complete data or entirely missing ##
##                                                                                                                                                ##
####################################################################################################################################################
load("~/EVAChallenge2019/DATA_TRAINING.RData")

check.missings <- function(i){
  print(paste('Site', i))
  x = anom.training[, i]
  out = matrix(NA, length(unique(year)), 12)
  
  for(j in 1:length(unique(year))){
    l = unique(year)[j]
    for(k in 1:12){
      # Observations in original scale for month {k} and year {j}
      x.kj = x[month == k & year == l]
      out[j,k] = mean(is.na(x.kj))
    }
  }
  # out
  sum(out > 0 & out < 1)
}

check.missings(1)

# out = sapply(1:16, check.missings)
t0 = Sys.time()
out = sapply(1:16703, check.missings)
t1 = Sys.time() - t0

range(unlist(out)) 
# [1] 0 0 # this means, either we have complete months, or we dont. We dont have partially missing months