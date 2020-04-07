### Script for implementing identifiability solutions for crossed effects 
### (as per Equation 14).

library(rjags)
load.module("dic")

### Load simulated (fake) data that represents observations of y (response variable) and
### x (covariate) for multiple (3) observations associated per each plot and date combination,
### where plot and date represent crossed effects:

data.eq14 = list(N = 120, Ndates = 10, Nplots = 4,
                 plot = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 
                          1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 
                          1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 
                          1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 
                          1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4), 
                 date = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                          3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                          5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
                          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
                          9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), 
                 x = c(-3.350, -3.193, 2.662, -4.679, -3.471, 2.115, -3.073, -0.632, -1.474, 4.537, 4.458, 
                       -1.491, 3.244, 2.334, 4.750, 4.618, 3.296, -1.454, -4.639, -2.613, -4.779, 4.133, -1.576, 
                       -3.061, -1.151, 2.288, 3.569, 1.473, -2.937, -2.273, -0.122, -1.307, -0.032, 0.439, 
                       -1.885, -4.275, -1.647, -4.979, 3.742, -0.366, 1.287, 0.586, -1.559, 2.982, 2.474, 
                       -2.850, -2.702, -4.275, 2.659, 3.267, 0.856, -1.034, 0.120, 3.111, -4.978, 0.891, 2.442, 
                       -0.242, 2.438, 3.180, -2.952, -3.623, 1.979, 4.600, 4.991, 1.199, -3.525, -2.043, -1.121, 
                       3.329, -3.917, 4.684, 0.313, 1.942, -0.528, 1.475, -2.777, -0.570, -0.253, 1.286, 1.325, 
                       -0.277, 4.660, -3.524, 1.623, -2.496, 4.272, 2.029, 1.669, -0.347, -0.315, -0.195, 3.570, 
                       4.169, -1.964, -3.184, -1.906, -0.590, -3.840, 1.165, 0.030, 0.788, -3.137, 3.002, 4.195, 
                       2.445, -3.157, -0.133, -4.354, 1.135, 4.819, -2.556, 2.244, 1.807, 2.303, 3.594, 3.421, 
                       2.793, 0.103, -4.272), 
                 y = c(-3.823, -4.313, 8.908, -5.854, -4.174, 5.956, -4.197, 1.623, -1.792, 11.878, 13.117, -0.847,
                       4.585, 3.371, 6.744, 6.907, 4.827, -4.012, -9.786, -5.785, -10.967, 7.356, -5.127, -6.277, 
                       -2.329, 7.268, 7.522, 3.536, -6.207, -2.836, 0.283, -2.259, 1.728, -0.902, -2.046, -6.591, 
                       -0.389, -7.209, 10.506, 4.537, 5.540, 4.587, -0.568, 9.887, 6.319, -1.607, -0.974, -7.014,
                       7.117, 9.580, 4.449, -0.030, 2.223, 8.054, -8.167, 4.948, 7.068, 1.086, 7.685, 8.979, -0.490, 
                       -2.455, 7.423, 13.535, 16.440, 4.937, -3.790, -1.693, 4.138, 9.554, -5.815, 11.434, -1.379, 
                       2.148, -2.836, -0.262, -6.659, -3.164, -2.828, 0.733, 2.105, -3.081, 8.056, -8.781, 7.205, 
                       -0.362, 12.744, 7.362, 5.597, 2.128, 3.556, 2.491, 11.016, 10.554, -0.465, -3.143, 2.484, 
                       4.099, -3.566, 6.175, 4.618, 4.723, -2.144, 10.899, 12.665, 10.948, -3.530, 4.026, -5.339, 
                       4.622, 13.212, -1.760, 6.782, 6.548, 7.115, 9.286, 7.429, 7.445, 3.029, -3.577))


##########################################################
### JAGS (OpenBUGS) code with sum-to-zero constraints (sz).
### Code uses the indexing variables of plot and date.
##########################################################

mod.string.sz <- "
  model{
    for(i in 1:N){
      y[i] ~ dnorm(mu[i],tau)
      # Mean model as per the original specification in Eq 14:
      mu[i] <- b0 + b1*x[i] + eps[plot[i]] + gamma[date[i]]
    }
    # Zero-centered priors for random effects with sum-to-zero contraints
    # Plot random effects:
    for(p in 1:(Nplots-1)){
      eps[p] ~ dnorm(0, tau.eps)
    }
    # Sum-to-zero constraint for plots:
    eps[Nplots] <- -sum(eps[1:(Nplots-1)])
    
    # Date random effects:
    for(d in 1:(Ndates-1)){
      gamma[d] ~ dnorm(0,tau.gamma)
    }
    # Sum-to-zero constraint for dates:
    gamma[Ndates] <- -sum(gamma[1:(Ndates-1)])
    
    # Relatively non-informative priors for root nodes:
    b0 ~ dnorm(0,1E-6)
    b1 ~ dnorm(0,1E-6)
    tau ~ dgamma(0.1, 0.1)
    tau.eps ~ dgamma(0.1, 0.1)
    tau.gamma ~ dgamma(0.1, 0.1)
    sig <- 1/sqrt(tau)
    sig.eps <- 1/sqrt(tau.eps)
    sig.gamma <- 1/sqrt(tau.gamma)
  }
"

##########################################################
### JAGS (OpenBUGS) code for post-sweeping of random 
### effects.
##########################################################

mod.string.ps <- " 
  model{
    for(i in 1:N){
      y[i] ~ dnorm(mu[i],tau)
      # Mean model as per the original specification in Eq 11:
      mu[i] <- b0 + b1*x[i] + eps[plot[i]] + gamma[date[i]]
    }
    # Zero-centered priors for non-identifiable random effects
    # Plot random effects:
    for(p in 1:Nplots){
      # Non-identifiable random effect (don't monitor)
      eps[p] ~ dnorm(0, tau.eps)
      # Identifiable random effect (monitor this)
      eps.star[p] <- eps[p] - ave.eps
    }
    # Mean plot-level random effect:
    ave.eps <- mean(eps[])
    
    # Date random effects:
    for(d in 1:Ndates){
      # Non-identifiable random effect
      gamma[d] ~ dnorm(0,tau.gamma)
      # Identifiable random effect (monitor this)
      gamma.star[d] <- gamma[d] - ave.gamma
    }
    # Mean date random effect:
    ave.gamma <- mean(gamma[])
    
    # Relatively non-informative priors for root nodes:
    # Non-identifiable intercept
    b0 ~ dnorm(0,1E-6)
    # Identifiable intercept (monitor this)
    b0.star <- b0 + ave.gamma + ave.eps
    b1 ~ dnorm(0,1E-6)
    tau ~ dgamma(0.1, 0.1)
    tau.eps ~ dgamma(0.1, 0.1)
    tau.gamma ~ dgamma(0.1, 0.1)
    sig <- 1/sqrt(tau)
    sig.eps <- 1/sqrt(tau.eps)
    sig.gamma <- 1/sqrt(tau.gamma)
  }
"


##########################################################
### Implement Bayesian models via JAGS

# Set arguments for jag.model function:
n.adapt = 5000
n.chains = 3
n.iter = 10000


##########################################################
### First: sum-to-zero constraint

mod.sz <- textConnection(mod.string.sz)
jm.sz =jags.model(mod.sz,
                  data=data.eq14,
                  n.chains=n.chains,
                  n.adapt=n.adapt)
coda.sz = coda.samples(jm.sz,variable.names=c("deviance","b0","b1","sig","sig.eps", "sig.gamma", "eps", "gamma"),
                       n.iter=n.iter)


##########################################################
### Final: post-sweeping of random effects

mod.ps <- textConnection(mod.string.ps)
jm.ps =jags.model(mod.ps,
                  data=data.eq14,
                  n.chains=n.chains,
                  n.adapt=n.adapt)
coda.ps = coda.samples(jm.ps,variable.names=c("deviance","b0.star","b1","sig","sig.eps", "sig.gamma", "eps.star", "gamma.star"),
                       n.iter=n.iter)


### MCMC plots
library(mcmcplots)
mcmcplot(window(coda.sz,thin=10))
mcmcplot(window(coda.ps,thin=10))

# Compute Raftery diagnostic to determine number of MCMC iterations required:
raft.sz<-raftery.diag(coda.sz)
raft.ps<-raftery.diag(coda.ps)


###################################################################
## Posterior statistics from each model:
## Compute posterior summary statistics, assuming a sufficient number
## of MCMC iterations have been run, otherwise, use coda.samples
## to update the MCMC iterations.


stats.sz <- summary(window(coda.sz), start=1000)
sum_sz <- cbind.data.frame(mean=stats.sz$statistics[,1],
                           lower2.5=stats.sz$quantiles[,1],
                           upper97.5=stats.sz$quantiles[,5])
sum_sz<-round(sum_sz,4)

stats.ps <- summary(window(coda.ps), start=1000)
sum_ps <- cbind.data.frame(mean=stats.ps$statistics[,1],
                           lower2.5=stats.ps$quantiles[,1],
                           upper97.5=stats.ps$quantiles[,5])
sum_ps<-round(sum_ps,4)