#####################################
## Backtransform to original scale ##
#####################################
Tx.inv <- function(y, xi, mu, sigma){ 
  if(is.na(y))                        
    return(NA)                        
  else{                               
    if(xi == 0)                       
      return((mu + sigma*y))          
    else                              
      return(mu + (sigma/xi)*(exp(xi*y) - 1)) 
  }                                             
}
Sx.inv <- function(z){
  if(is.na(z))                                  
    return(NA)                                  
  else{                                         
    v <- pnorm(z)
    return(qexp(v))                          
  }                                           
}                                             

backT <- function(z, xi, mu, sigma){
  y = Sx.inv(z)
  return(Tx.inv(y, xi, mu, sigma))
}

##########################################################
## Function calculating the twCRPS (from Competition.R) ##
##########################################################
# INPUTS:
# prediction: matrix (dimension 162000x400, size 494.4 Mb). The i-th row contains the predicted distribution for the i-th point in the validation set, evaluated at each of the 400 design points
# true.observations: vector (length 162000, size 1.2 Mb). Observation vector (unknown to the teams), containing the true spatio-temporal minimum X(s,t)=min A(s,t) for each point in the validation set
twCRPS <- function(prediction,true.observations){
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