####################################
# STA421 FBM FS22: 05Worksheet Exercise 3
####################################

###################
# Normal Example (Gibbs MCMC)
###################

remove(list=ls())
set.seed(44566)

# set the path to the 05normal_exmple_JAGS.txt file
path <- ""

##########################################################
# rjags interface to JAGS
##########################################################

library(rjags)

list.factories(type = "rng")
list.factories(type = "monitor")

list.factories(type = "sampler")
set.factory(name = "base::Slice", type = "sampler", state = FALSE)
list.factories(type = "sampler")
set.factory(name = "base::Slice", type = "sampler", state = TRUE)
list.factories(type = "sampler")

list.modules()
load.module("glm")
list.modules()
unload.module("glm")
list.modules()


####################
## Introduction
####################


# generating data
#mu <- 4
#sigma2 <- 16
#n <- 30
#y <- rnorm(n=n, mean=mu, sd=sqrt(sigma2))


# Load the data
# round(y,3)
y <- c(3.048,2.980,2.029,7.249,-0.259,3.061,4.059,6.370,7.902,1.926,
       9.094,10.489,-0.384,-3.096,2.315,5.830,-1.542,-1.544,5.714,
       -5.182,3.828,-4.038,2.169,5.087,-0.201,4.880,3.302,3.859,
       11.144,5.564)

par(mfrow = c(1, 1))
boxplot(y)
summary(y)
sd(y)

# Define the parameters of the prior distributions 
mu0 <- -3
sigma2_0 <- 4
a0 <- 1.6
b0 <- 0.4


################################
# INLA exact result (motivation)
################################

library(INLA)

formula <- y ~ 1
inla.output <- inla(formula,data=data.frame(y=y),
                    control.family = list(hyper =
                                            list(prec = list(prior="loggamma",param=c(a0,b0)))),
                    control.fixed = list(mean.intercept=mu0, prec.intercept=1/sigma2_0))


#############################
# Step 2: JAGS model file as a string in rjags with coda
#############################
set.seed(44566)

library(rjags)
library(coda)


#sessionInfo()

wb_data <- list( N=30, 
                 y=c(3.048,2.980,2.029,7.249,-0.259,3.061,4.059,6.370,7.902,1.926,
                     9.094,10.489,-0.384,-3.096,2.315,5.830,-1.542,-1.544,5.714,
                     -5.182,3.828,-4.038,2.169,5.087,-0.201,4.880,3.302,3.859,
                     11.144,5.564) 
)

wb_inits <- list( mu=-0.2381084, inv_sigma2=0.3993192 )

modelString = " # open quote for modelString
model{
# likelihood
for (i in 1:N){ 
y[i] ~ dnorm( mu, inv_sigma2 )    
}
# Priors
mu ~ dnorm( -3, 0.25 ) # prior for mu N(mu0, prec=1/sigma2_0)
inv_sigma2 ~ dgamma( 1.6, 0.4 ) # prior for precision G(a0, b0)

# transformations
# deterministic definition of variance
sigma2 <- 1/inv_sigma2

# deterministic definition of standard deviation
sigma <- sqrt(sigma2)
}
" # close quote for modelString

writeLines(modelString, con="TempModel.txt") # write to a file

##########################################################
# JAGS only one chain
##########################################################

# model initiation
model.jags <- jags.model(
  file = "TempModel.txt", 
  data = wb_data,
  inits = wb_inits,
  n.chains = 1,
  n.adapt = 4000
)

str(model.jags)
class(model.jags)
attributes(model.jags)
list.samplers(model.jags)

# burn-in

update(model.jags, n.iter = 4000)

# sampling
fit.jags.coda <- coda.samples(
  model = model.jags, 
  variable.names = c("mu", "sigma2", "inv_sigma2"), 
  n.iter = 10000,
  thin = 1
)

str(fit.jags.coda)
class(fit.jags.coda)
attributes(fit.jags.coda)

summary(fit.jags.coda)
#print(fit.jags.coda)
plot(fit.jags.coda)


# store samples for each parameter from the chain into separate objects
m.fit.jags.coda <- as.matrix(fit.jags.coda)
mu.sim <- m.fit.jags.coda[,"mu"] 
sigma2.sim <- m.fit.jags.coda[,"sigma2"]
inv_sigma2.sim <- m.fit.jags.coda[,"inv_sigma2"] 

library(MASS)

par(mfrow=c(1,1))
# plot for mean
rg <- range(inla.output$marginals.fixed$"(Intercept)"[,2])
truehist(mu.sim, prob=TRUE, col="yellow", xlab=expression(mu),ylim=rg)
lines(density(mu.sim),lty=3,lwd=3, col=2)
lines(inla.output$marginals.fixed$"(Intercept)",lwd=2)
legend("topright",c("MCMC: JAGS","INLA"),lty=c(3,1),lwd=c(2,2),col=c(2,1),cex=1.0,bty="n")

# plot for variance
m_var <-inla.tmarginal(function(x) 1/x, inla.output$marginals.hyperpar[[1]])
rg <- range(m_var[,2])
truehist(sigma2.sim, prob=TRUE, col="yellow", xlab=expression(sigma^2),ylim=rg)
lines(density(sigma2.sim),lty=3,lwd=3, col=2)
lines(m_var,lwd=2)
legend("topright",c("MCMC: JAGS","INLA"),lty=c(3,1),lwd=c(2,2),col=c(2,1),cex=1.0,bty="n")

# plot for precision
truehist(inv_sigma2.sim, prob=TRUE, col="yellow", xlab=expression(1/sigma^2))
lines(density(inv_sigma2.sim),lty=3,lwd=3, col=2)
lines(inla.output$marginals.hyperpar[[1]],lwd=2)
legend("topright",c("MCMC: JAGS","INLA"),lty=c(3,1),lwd=c(2,2),col=c(2,1),cex=1.0,bty="n")














##########################################################
# JAGS several chains
##########################################################

wb_inits <- function() {
  list(mu = rnorm(1),
       inv_sigma2 = runif(1)
  )  
}

# model initialisation
model.jags <- jags.model(
  file = "TempModel.txt", 
  data = wb_data,
  inits = wb_inits,
  n.chains = 4,
  n.adapt = 4000
)

# burn-in

update(model.jags, n.iter = 4000)

# sampling/monitoring
fit.jags.coda <- coda.samples(
  model = model.jags, 
  variable.names = c("mu", "sigma2", "inv_sigma2"), 
  n.iter = 10000,
  thin = 10
)

#n.thin<-floor((n.iter-n.adapt)/500)
#floor((10000-4000)/500)=12

summary(fit.jags.coda)
plot(fit.jags.coda)



# store samples for each parameter from the chains into separate vectors
m.fit.jags.coda <-as.matrix(fit.jags.coda)
mu.sim <- m.fit.jags.coda[,"mu"] 
sigma2.sim <- m.fit.jags.coda[,"sigma2"]
inv_sigma2.sim <- m.fit.jags.coda[,"inv_sigma2"] 


par(mfrow=c(1,1))
# plot for mean
rg <- range(inla.output$marginals.fixed$"(Intercept)"[,2])
truehist(mu.sim, prob=TRUE, col="yellow", xlab=expression(mu),ylim=rg)
lines(density(mu.sim),lty=3,lwd=3, col=2)
lines(inla.output$marginals.fixed$"(Intercept)",lwd=2)
legend("topright",c("MCMC: JAGS","INLA"),lty=c(3,1),lwd=c(2,2),col=c(2,1),cex=1.0,bty="n")

# plot for variance
m_var <-inla.tmarginal(function(x) 1/x, inla.output$marginals.hyperpar[[1]])
rg <- range(m_var[,2])
truehist(sigma2.sim, prob=TRUE, col="yellow", xlab=expression(sigma^2),ylim=rg)
lines(density(sigma2.sim),lty=3,lwd=3, col=2)
lines(m_var,lwd=2)
legend("topright",c("MCMC: JAGS","INLA"),lty=c(3,1),lwd=c(2,2),col=c(2,1),cex=1.0,bty="n")

# plot for precision
truehist(inv_sigma2.sim, prob=TRUE, col="yellow", xlab=expression(1/sigma^2))
lines(density(inv_sigma2.sim),lty=3,lwd=3, col=2)
lines(inla.output$marginals.hyperpar[[1]],lwd=2)
legend("topright",c("MCMC: JAGS","INLA"),lty=c(3,1),lwd=c(2,2),col=c(2,1),cex=1.0,bty="n")


## CODA
#summary(fit.jags.coda)
#effectiveSize(fit.jags.coda)
#lapply(fit.jags.coda, effectiveSize)
#gelman.diag(fit.jags.coda,autoburnin=TRUE)
#gelman.plot(fit.jags.coda,autoburnin=TRUE)
#geweke.diag(fit.jags.coda)
#geweke.plot(fit.jags.coda)
#heidel.diag(fit.jags.coda)
#raftery.diag(fit.jags.coda)
#coda:::traceplot(fit.jags.coda)


# "DIC" penalised expected deviance computation
dic1<-dic.samples(model=model.jags, n.iter=1000, type="popt")
#dic2<-dic.samples(model=model.jags2, n.iter=1000, type="popt")

# "DIC" penalised expected deviance comparison
# There is no absolute scale for DIC comparison
# SE is very helpful

#diffdic(dic1,dic2)









######
# Additional sampling in several chains, preparation for BGR/Gelman
# with runjags
######

library(runjags)

###############
# runjags interface with a link to a file
###############


wb_data <- list( N=30, 
                 y=c(3.048,2.980,2.029,7.249,-0.259,3.061,4.059,6.370,7.902,1.926,
                     9.094,10.489,-0.384,-3.096,2.315,5.830,-1.542,-1.544,5.714,
                     -5.182,3.828,-4.038,2.169,5.087,-0.201,4.880,3.302,3.859,
                     11.144,5.564) 
)

wb_inits <- function() {
  list(mu = rnorm(1),
       inv_sigma2 = runif(1)
  )  
}

fit.runjags<-run.jags(model=paste(path,"05normal_exmple_JAGS.txt",sep=""),
                      monitor=c("mu", "sigma2", "inv_sigma2"),
                      data=wb_data,
                      inits=wb_inits,
                      n.chains=4,
                      burnin=4000,
                      sample=5000,
                      adapt=1000,
                      thin=2)


plot(fit.runjags)
print(fit.runjags)

# CODA
fit.runjags.coda<-as.mcmc.list(fit.runjags)
summary(fit.runjags.coda)
# conduct CODA

















##########################################################
# R2jags wrapper to rjags interface to JAGS several chains
##########################################################



library(R2jags)

#rm(list=ls())


###############
# R2jags wrapper with a link to a file
###############


wb_data <- list( N=30, 
                 y=c(3.048,2.980,2.029,7.249,-0.259,3.061,4.059,6.370,7.902,1.926,
                     9.094,10.489,-0.384,-3.096,2.315,5.830,-1.542,-1.544,5.714,
                     -5.182,3.828,-4.038,2.169,5.087,-0.201,4.880,3.302,3.859,
                     11.144,5.564) 
)

#define parameters
params<-c("mu", "sigma2", "inv_sigma2")

# define inits
inits1 <- list(mu=rnorm(1), inv_sigma2=runif(1),
               .RNG.name="base::Super-Duper", .RNG.seed=1)
inits2 <- list(mu=rnorm(1), inv_sigma2=runif(1),
               .RNG.name="base::Wichmann-Hill", .RNG.seed=2)
wb_inits <- list(inits1,inits2)


fit.R2jags<-jags(data=wb_data,
                 inits=wb_inits,
                 parameters.to.save=params,
                 model.file="05normal_exmple_JAGS.txt",
                 n.chains=2,
                 n.iter=50000,
                 n.burnin=4000,
                 n.thin=5,
                 DIC = TRUE,
                 jags.seed = 321,
                 refresh =100,
                 digits = 4,
                 jags.module = c("glm","dic"))


# Standard plots of the monitored variables
plot(fit.R2jags)

# Display summary statistics
print(fit.R2jags)

# traceplot
traceplot(fit.R2jags)


# CODA
fit.R2jags.coda<-as.mcmc(fit.R2jags)
summary(fit.R2jags.coda)
# conduct CODA

