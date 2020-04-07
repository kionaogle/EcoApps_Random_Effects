### Script for implementing identifiability solutions for groups of nested random 
### effects (as per Equation 11).

library(rjags)
load.module("dic")

### Load simulated (fake) data that represents observations of y (response variable) and
### x (covariate) for multiple (3) observations associated with each plot, with plots 
### nested in watersheds:
data.eq11 = list(N = 105, Nwater = 7, Nplots = 5, 
                 plot = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 
                          1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 
                          1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 
                          1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5), 
                 wshed = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                               3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                               5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
                               7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7), 
                 x = c(0.030, 2.201, 4.777, 3.590, -1.892, 3.764, -2.600, -1.973, -3.739, 2.547, 1.888, 3.408, -3.231, 
                       0.930, 2.214, 4.051, -2.248, -0.144, -0.159, -1.048, -2.123, 1.342, 2.536, 0.241, -0.618, 1.678, 
                       2.576, -3.336, -3.280, 0.599, 2.286, -2.154, 3.014, 1.983, -2.686, -2.542, -4.263, 1.113, 0.133,
                       4.064, 0.128, 3.291, 4.773, 3.088, -2.181, 3.626, 3.895, 1.392, 0.517, -1.047, -1.941, 2.757, 4.375, 
                       -1.039, 0.171, 1.942, 0.707, -4.835, -1.762, -2.230, 3.346, -1.102, 2.854, -3.667, -0.825, -3.787,
                       -0.089, -1.860, 3.084, -4.438, 2.000, -0.869, -4.745, 4.307, 0.168, 4.757, -1.877, 4.731, 3.006, 3.992, 
                       -1.540, 0.985, -4.165, 4.883, 0.599, 4.914, 4.144, 3.160, 2.791, -2.406, 2.599, -0.139, 3.054, 2.459, 
                       -3.216, -2.750, 2.721, 4.258, -3.828, 4.390, 2.842, 4.650, 1.986, -4.539, -0.360), 
                 y = c(1.306, 6.968, 11.674, 9.677, -1.654, 9.273, -4.454, -4.087, -8.358, 7.448, 3.829, 8.724, -3.131, 5.386,
                       7.477, 12.766, 0.412, 4.492, 3.907, 0.672, -1.564, 7.118, 11.959, 5.734, 3.839, 8.457, 9.043, -4.891, 
                       -3.259, 3.876, 4.696, -3.645, 6.773, 5.058, -5.250, -5.224, -8.932, 1.449, 1.277, 10.621, 1.183, 8.481, 
                       9.286, 5.636, -6.161, 3.697, 5.455, 2.329, 2.576, -2.266, -3.676, 4.984, 6.979, -3.226, -1.312, 2.612, 
                       -0.584, -13.626, -6.124, -8.599, 6.349, -2.421, 3.217, -5.498, -0.018, -6.047, -0.291, -3.239, 7.361, 
                       -7.387, 4.761, -1.064, -10.253, 7.422, -0.743, 12.475, -1.469, 12.926, 7.778, 10.915, 0.077, 4.088, 
                       -5.555, 12.275, 2.214, 12.629, 9.342, 9.360, 10.634, -1.614, 4.205, -1.578, 5.332, 4.696, -7.734, 
                       -7.363, 6.005, 9.604, -5.809, 9.230, 6.116, 9.129, 5.805, -7.872, -2.073)) 

##########################################################
### JAGS (OpenBUGS) code for hierarchical centering (hc).
### Code uses the indexing variables of plot and wshed
##########################################################

mod.string.hc <- "
model{
  for(i in 1:N){
    y[i] ~ dnorm(mu[i],tau)
    # Mean model differs from the original specification in Eq 11 to reflect hierarhical centering:
    mu[i] <- b0.eps[plot[i],wshed[i]] + b1*x[i]
  }
  # Hierarchical priors for random effects, centered on overall intercept (b0):
  for(s in 1:Nwater){
    for(p in 1:Nplots){
      # Plot within watershed effects:
      b0.eps[p,s] ~ dnorm(b0.gamma[s], tau.eps)
    }
    # Watershed effects:
    b0.gamma[s] ~ dnorm(b0, tau.gamma)
    }
    
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
### JAGS (OpenBUGS) code with sum-to-zero constraints (sz).
##########################################################

mod.string.sz <- "
model{
  for(i in 1:N){
    y[i] ~ dnorm(mu[i],tau)
    # Mean model as per the original specification in Eq 11:
    mu[i] <- b0 + b1*x[i] + eps[plot[i], wshed[i]] + gamma[wshed[i]]
  }
  # Zero-centered priors for random effects with sum-to-zero contraint
  # Plot within watershed effects:
  for(s in 1:Nwater){
    for(p in 1:(Nplots-1)){
      eps[p,s] ~ dnorm(0, tau.eps)
    }
    # Sum-to-zero occurs for plots within each watershed:
    eps[Nplots,s] <- -sum(eps[1:(Nplots-1),s])
  }
  # Watershed random effects:
  for(s in 1:(Nwater-1)){
    gamma[s] ~ dnorm(0,tau.gamma)
  }
  # Sum-to-zero constraint:
  gamma[Nwater] <- -sum(gamma[1:(Nwater-1)])
  
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
      mu[i] <- b0 + b1*x[i] + eps[plot[i], wshed[i]] + gamma[wshed[i]]
    }
    # Zero-centered priors for non-identifiable random effects
    # Plot within watershed effects:
    for(s in 1:Nwater){
      for(p in 1:Nplots){
        # Non-identifiable random effect (don't monitor):
        eps[p,s] ~ dnorm(0, tau.eps)
        # Identifiable plot random effect (monitor this):
        eps.star[p,s] <- eps[p,s] - ave.eps[s]
      }
      # Mean plot-level random effects within each watershed:
      ave.eps[s] <- mean(eps[,s])
    }
    # Watershed random effects:
    for(s in 1:Nwater){
      # Non-identifiable random effect:
      gamma[s] ~ dnorm(0,tau.gamma)
      # Identifiable random effect (monitor this):
      gamma.star[s] <- gamma[s] + ave.eps[s] - ave.gamma - ave.ave.eps
    }
    # Mean watershed random effect:
    ave.gamma <- mean(gamma[])
    # Mean overall plot random effect:
    ave.ave.eps <- mean(ave.eps[])
    
    # Relatively non-informative priors for root nodes:
    # Non-identifiable intercept
    b0 ~ dnorm(0,1E-6)
    # Identifiable intercept (monitor this)
    b0.star <- b0 + ave.gamma + ave.ave.eps
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

### First: hierarhical centering:
mod.hc <- textConnection(mod.string.hc)
jm.hc =jags.model(mod.hc,
                  data=data.eq11,
                  n.chains=n.chains,
                  n.adapt=n.adapt)
coda.hc = coda.samples(jm.hc,variable.names=c("deviance","b0", "b1", "sig", "sig.eps", "sig.gamma", "b0.eps", "b0.gamma"),
                       n.iter=n.iter)

##########################################################
### Next: sum-to-zero constraint

mod.sz <- textConnection(mod.string.sz)
jm.sz =jags.model(mod.sz,
                  data=data.eq11,
                  n.chains=n.chains,
                  n.adapt=n.adapt)
coda.sz = coda.samples(jm.sz,variable.names=c("deviance","b0","b1","sig","sig.eps", "sig.gamma", "eps", "gamma"),
                       n.iter=n.iter)


##########################################################
### Final: post-sweeping of random effects

mod.ps <- textConnection(mod.string.ps)
jm.ps =jags.model(mod.ps,
                  data=data.eq11,
                  n.chains=n.chains,
                  n.adapt=n.adapt)
coda.ps = coda.samples(jm.ps,variable.names=c("deviance","b0.star","b1","sig","sig.eps", "sig.gamma", "eps.star", "gamma.star"),
                       n.iter=n.iter)


### MCMC plots
library(mcmcplots)
mcmcplot(window(coda.hc,thin=10))
mcmcplot(window(coda.sz,thin=10))
mcmcplot(window(coda.ps,thin=10))

# Compute Raftery diagnostic to determine number of MCMC iterations required:
raft.hc<-raftery.diag(coda.hc)
raft.sz<-raftery.diag(coda.sz)
raft.ps<-raftery.diag(coda.ps)


###################################################################
## Posterior statistics from each model:
## Compute posterior summary statistics, assuming a sufficient number
## of MCMC iterations have been run, otherwise, use coda.samples
## to update the MCMC iterations.

stats.hc <- summary(window(coda.hc), start=1000)
sum_hc <- cbind.data.frame(mean=stats.hc$statistics[,1],
                           lower2.5=stats.hc$quantiles[,1],
                           upper97.5=stats.hc$quantiles[,5])
sum_hc<-round(sum_hc,4)

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