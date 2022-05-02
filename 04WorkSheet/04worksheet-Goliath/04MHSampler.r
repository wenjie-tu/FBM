####################################
# STA421 FBM FS22: 04Worksheet Exercise 3
####################################


## Metropolis-Hastings for logistic model
## Two independent normal proposals

rm(list=ls())

## Data from Collett, D. (2003) Modelling Binary Data 2nd Edition
## Table 1.6 p. 7 Number of deaths from pneumonia amongst batches of 40 mice exposed to different doses of an anti-pneumococcus serum (grouped)
## Table 3.1 p. 71 (original) Number of deaths from pneumonia in mice exposed to various doses of an anti-pneumococcus serum

# the covariate values (dose)
# x_original <- c(0.0028, 0.0056, 0.0112, 0.0225, 0.0450)
x_original <- c(0.0028, 0.0028, 0.0056, 0.0112, 0.0225, 0.0450)
# the centered covariate values (centered dose)
x <- x_original - mean(x_original)
# number of mice deaths
# y <- c(35, 21, 9, 6, 1)
y <- c(26, 9, 21, 9, 6, 1)
# total number of mice
# n <- c(40, 40, 40, 40, 40)
n <- c(28, 12, 40, 40, 40, 40)

# Assumption
# variance of normal priors
sigma2 <- 10^(4)



####################
## Bayesian analysis
####################

# inverse logit: logit^(-1)(alpha + beta*x)
mypi <- function(alpha, beta, x){
  tmp <- exp(alpha + beta*x)
  pi <- tmp/(1+tmp)
  return(pi)
}



#####################################
## Step 1: R: (univariate proposal) Metropolis MCMC settings 
#####################################

set.seed(44566)

# number of MCMC iterations
n.iter <- 10000 
# burnin length
n.burnin <- 4000
# thinning parameter
n.thin <- 1
#n.thin <- floor((n.iter-n.burnin)/500)


######################################
## univariate random walk proposals ##
######################################


alpha_samples <- c()
beta_samples <- c()
# number of accepted proposals
alpha_yes <- 0
beta_yes <- 0

# starting values
alpha <- 0
beta <- 0


# standard deviations for the
# normal proposal
s_alpha <- 1
s_beta <- 60


# counter
count <- 0

# start the MCMC algorithm (the first iteration after the burn-in is 1)
for(i in -n.burnin:(n.iter*n.thin)){
  count <- count +1
  
  ## update alpha
  # generate a new proposal for alpha
  alpha_star <- rnorm(1, alpha, sd=s_alpha)
  
  # NOTE: it is more stable to calculate everything on the log scale
  enum <- sum(dbinom(y, size=n, prob=mypi(alpha_star, beta, x), log=TRUE)) + 
    dnorm(alpha_star, mean=0, sd=sqrt(sigma2), log=TRUE)
  denom <- sum(dbinom(y, size=n, prob=mypi(alpha, beta, x), log=TRUE))  + 
    dnorm(alpha, mean=0, sd=sqrt(sigma2), log=TRUE)
  
  # log accetpance rate (since we use a random walk proposal there is no
  #	proposal ratio in the acceptance probability)
  logacc <- enum - denom
  if(log(runif(1)) <= logacc){
    # accept the proposed value
    alpha <- alpha_star
    alpha_yes <- alpha_yes + 1
  }
  
  ## update beta
  # generate a new proposal for beta
  beta_star <- rnorm(1, beta, sd=s_beta)
  
  enum <- sum(dbinom(y, size=n, prob=mypi(alpha, beta_star, x), log=TRUE)) + 
    dnorm(beta_star, mean=0, sd=sqrt(sigma2), log=TRUE)
  denom<- sum(dbinom(y, size=n, prob=mypi(alpha, beta, x), log=TRUE)) + 
    dnorm(beta, mean=0, sd=sqrt(sigma2), log=TRUE)
  # log accetpance rate
  logacc <- enum - denom
  
  if(log(runif(1)) <= logacc){
    # accept the proposed value
    beta <- beta_star
    beta_yes <- beta_yes + 1
  }
  
  # after the burnin save every kth sample
  if((i > 0) && (i%%n.thin == 0)){
    alpha_samples <- c(alpha_samples, alpha)
    beta_samples <- c(beta_samples, beta)
  }
  if(i%%1000 ==0){
    # print the acceptance rates on the fly
    cat(c(i, alpha_yes/count, beta_yes/count), "\n")
  }
}




