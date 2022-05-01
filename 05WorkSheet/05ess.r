####################################
# STA421 FBM FS22: 05Worksheet Exercise 6
####################################


#################################################################
## R-code to compute the effective sample size for a MCMC chain
## According to Geyer (1992)
## Adapted from STA422 in FS12
#################################################################

library(coda)

# Function to test if a vector is monoton decreasing,
# a boolean value is returned
monotone <- function(vec){

  a <- TRUE
  
  if(length(vec) == 1){
    return(a)
  }

  for(i in 2:length(vec)){
    if(vec[i] > vec[i-1]){
      a <- FALSE
      break;
    }
  }
  
  return(a)
}

# Function to find the lag to stop the ESS
# calculation compare: 
#
#      Geyer (1992),
#     "Practical Markov Chain Monte Carlo".
#      Statistical Science, 7: 473- 511
#
# Gamma_i = gamma(2*i) + gamma(2*i + 1)
#
# m is the greatest integer at which Gamma_i > 0
# and Gamma_i is montone for i = 1, ..., m.
# Thereby gamma(i) is the sample autocorrelation 
# at lag i. 
#
# Parameter:
# vec - sample vector (mcmc object)
#
# Output
# m <- greatest integer where both criteria are
#     fulfilled
geyer <- function(vec){

  g <- c()
  res <- 1

  for(i in 1:(length(vec)/2 - 1)){
    g <- c(g, 
      autocorr(vec, lags = 2*i) + autocorr(vec, lags = 2*i + 1))
    if(monotone(g) == FALSE || g[i] < 0){
      break
    }
  }
  if(i==1){
    res <- 1
  }
  else{
    res <- i-1
  }
  return(res)
}

# Function to calculate the effective sample
# size for one MCMC chain.
#
# Parameter:
# mcmc - mcmc object 
# M - number of sampled values
#
# Output:
# effective sample size
ess <- function(mcmc, M){

  m <- geyer(mcmc)
  y <- M / (1 + 2 * sum(autocorr(mcmc, lag = 1:(2*m +1))))
  return(y)
}


## Example
# M <- 10000
# a <- rnorm(M)
# my.mcmc <- as.mcmc(a)
# ess(my.mcmc, 10000)
