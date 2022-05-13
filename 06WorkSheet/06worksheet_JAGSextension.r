####################################
# STA421 FBM FS22: 06Worksheet Exercise 4
####################################


library(rjags)
library(coda)

###############
# rjags interface with code in a model string
###############


pl1.data<-list(N = 16, 
               y = c(23., 12., 19.,  9.,  39.,  6.,  9., 10., 120., 18., 107., 26., 82., 16., 126., 23.),
               n = c(107., 44., 51., 39., 139., 20., 78., 35., 208., 38., 150., 45., 138., 20., 201., 34.),
               C1 = c(0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 1.))


pl1.params<-c("mu", "beta", "tau", "p1.star", "p2.star")


pl1_modelString <- "
model 
{

#	sampling model (likelihood)
for (j in 1:N)	{
y[j] ~ dbin(p[j],n[j])
logit(p[j]) <- mu + beta*C1[j] + eta[j]
eta[j] ~ dnorm(0, tau.prec)

#	prediction for posterior predictive checks
y.pred[j] ~ dbin(p[j],n[j])
PPC[j] <- step(y[j]-y.pred[j])-0.5*equals(y[j],y.pred[j])
}

#	priors
mu ~ dunif(-10,10)
beta ~ dunif(-10,10)
tau ~ dunif(0,10)
tau.prec <- 1/tau/tau

#	population effect
p1 <- 1/(1+exp(-mu)) 
p2 <- 1/(1+exp(-mu-beta))

#	predictive distribution for new study effect
eta.star ~ dnorm(0,tau.prec)
p1.star <- 1/(1+exp(-mu-eta.star))
p2.star <- 1/(1+exp(-mu-beta-eta.star))

}
"

writeLines(pl1_modelString, con="TempModel.txt") # write to a file

# model initiation
rjags.pl1 <- jags.model(
  file = "TempModel.txt", 
  data = pl1.data,
  n.chains = 4,
  n.adapt = 4000
)

# str(rjags.pl1)
# class(rjags.pl1)
# attributes(rjags.pl1)

# burn-in

update(rjags.pl1, n.iter = 4000)

# sampling/monitoring
fit.rjags.pl1.coda <- coda.samples(
  model = rjags.pl1, 
  variable.names = pl1.params, 
  n.iter = 10000,
  thin = 1
)


summary(fit.rjags.pl1.coda)
plot(fit.rjags.pl1.coda)


# CODA