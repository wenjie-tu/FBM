## load packages:
library(rstan)
library(bayesplot)


## Sample data from normal distribution with mu =3 and sigma^2 = 1000
Y <- rnorm(n = 100, mean = 3, sd = 10)


## make list object out of sample data:
lst_score_data <- list(y = Y, N = length(Y))

# Note please setwd(), the normal.stan file needs to be in the same place as this R file.
fit_score <- stan(file = "normal.stan", data = lst_score_data)

# Make Traceplot:
traceplot(fit_score, pars = c("mu", "sigma"))

# Print Summary result:
print(fit_score, pars = c("mu", "sigma"))



# generate data frame from model fit: 
df_fit_score <- as.data.frame(fit_score)

## color scheme:
color_scheme_set("red")

# plot posterior distributions: 
mcmc_hist(df_fit_score, pars = c("mu", "sigma"))

