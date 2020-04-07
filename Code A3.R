### Script for implementing hierarchical centering, sum-to-zero, and
### post-sweeping of random effects associated with the random effects
### regression in Equation 5 in the main text, and using the simulated
### data from Data S1. See the JAGS code in Code A2 for the original,
### non-identifiable model for additional annotation / comments.

# Load data: x = x variable, y1 = y variable with sig.eps = 0.1
# y2 = y with sig.eps = 1, and y3 = y with sig.eps = 10.
load("Data S1.RData")

# Create data lists for JAGS (or OpenBUGS) model:
# Small sig.eps:
data1 = list(y = y1, x = x)
# If using folded-Cauchy prior:
#data1 = list(y = y1, x = x, a = 0.5, b = 1/2)

# Medium sig.eps:
data2 = list(y = y2, x = x)
# If using folded-Cauchy prior:
#data2 = list(y = y2, x = x, a = 0.5, b = 1/2)

# Large sig.eps:
data3 = list(y = y3, x = x)
# If using folded-Cauchy prior:
#data3 = list(y = y3, x = x, a = 0.5, b = 1/2)

library(rjags)
load.module("dic")

##########################################################
### JAGS (OpenBUGS) code for hierarchical centering (hc).
##########################################################

mod.string.hc <- "
model{
  for(j in 1:10){
    for(r in 1:10){
      y[r,j] ~ dnorm(mu[r,j],tau)
      # Mean model differs from the original, non-identifiable model:
      mu[r,j] <- b0.eps[j] + b1*x[r,j]
    }
    # Hierarchical prior for random effects, centered on overall intercept (mean):
    b0.eps[j] ~ dnorm(b0,tau.eps)
    # If interested in monitoring the random effects, can
    # compute as a deviation term (but, not necessary)
    eps[j] <- b0.eps[j] - b0
  }
  # Relatively non-informative priors for root nodes:
  b0 ~ dnorm(0,1E-6)
  b1 ~ dnorm(0,1E-6)
  tau ~ dgamma(0.1, 0.1)
  sig <- 1/sqrt(tau)
  # Mean of random effects (not necessary):
  mean.eps <- mean(eps[])

  # Conjugate gamma prior for random effects precision:
  tau.eps ~ dgamma(0.1, 0.1)
  sig.eps <- 1/sqrt(tau.eps)

  # Folded Cauchy(0,A^2,df=1) prior for standard deviation:
  # alpha ~ dnorm(0,1)
  # tau.temp ~ dgamma(a,b)
  # sig.temp <- 1/sqrt(tau.temp)
  # sig.eps <- abs(alpha)*sig.temp
  # tau.eps <- pow(sig.eps,-2)

}
"

##########################################################
### JAGS (OpenBUGS) code with sum-to-zero constraints (sz).
##########################################################

mod.string.sz <-"
  model{
    for(j in 1:10){
      for(r in 1:10){
        y[r,j] ~ dnorm(mu[r,j],tau)
        # Mean model the same as the original, non-identifiable model:
        mu[r,j] <- b0 + b1*x[r,j] + eps[j]
      }
    }
    
    # Hierarchical prior for all but one random effect:
    for(j in 1:9){
      eps[j] ~ dnorm(0,tau.eps)
    }
    # Remaining random effect = minus sum of other random effects:
    eps[10] <- -sum(eps[1:9])
    
    # Relatively non-informative priors for root nodes:
    b0 ~ dnorm(0,1E-6)
    b1 ~ dnorm(0,1E-6)
    tau ~ dgamma(0.1, 0.1)
    sig <- 1/sqrt(tau)
    # Mean of the random effects (not necessary), which will be exactly zero:
    mean.eps <- mean(eps[])
    
    # Conjugate gamma prior for random effects precision:
    tau.eps ~ dgamma(0.1, 0.1)
    sig.eps <- 1/sqrt(tau.eps)
    
    # Folded Cauchy(0,A^2,df=1) prior for standard deviation:
    # alpha ~ dnorm(0,1)
    # tau.temp ~ dgamma(a,b)
    # sig.temp <- 1/sqrt(tau.temp)
    # sig.eps <- abs(alpha)*sig.temp
    # tau.eps <- pow(sig.eps,-2)
    
  }
"


##########################################################
### JAGS (OpenBUGS) code for post-sweeping of random 
### effects (ps).
##########################################################

mod.string.ps <- "
  model{
    for(j in 1:10){
      for(r in 1:10){
        y[r,j] ~ dnorm(mu[r,j],tau)
        # Mean model same as original, non-identifiable model:
        mu[r,j] <- b0 + b1*x[r,j] + eps[j]
      }
      # Zero-centered hierarchical prior for random effects, just
      # as used in the original model (do not monitor or report these):
      eps[j] ~ dnorm(0,tau.eps)
      # Compute identifiable random effects (monitor and report these,
      # if desired):
      eps.star[j] <- eps[j] - mean.eps
    }
    # Prior for non-identifiable intercept (don't monitor or report)
    b0 ~ dnorm(0,1E-6)
    # Compute identifiable intercept (monitor and report this):
    b0.star <- b0 + mean.eps
    # Relatively non-informative priors for other root nodes:
    b1 ~ dnorm(0,1E-6)
    tau ~ dgamma(0.1, 0.1)
    sig <- 1/sqrt(tau)
    # Mean of non-identifiable random effects (required)
    mean.eps <- mean(eps[])
    # Mean of identifiable random effects (not required), will be exactly zero.
    mean.eps.star <- mean(eps.star[])
    
    # Conjugate gamma prior for random effects precision:
    tau.eps ~ dgamma(0.1, 0.1)
    sig.eps <- 1/sqrt(tau.eps)
    
    # Folded Cauchy(0,A^2,df=1) prior for standard deviation:
    # alpha ~ dnorm(0,1)
    # tau.temp ~ dgamma(a,b)
    # sig.temp <- 1/sqrt(tau.temp)
    # sig.eps <- abs(alpha)*sig.temp
    # tau.eps <- pow(sig.eps,-2)
  }
"


##########################################################
### Implement Bayesian models 

# Pick the dataset to be analyzed:
# for "small" sig.eps
data.soln = data1
# for "medium" sig.eps
data.soln = data2
# for "large" sig.eps
data.soln = data3

# Set arguments for jag.model function:
n.adapt = 5000
n.chains = 3
n.iter = 10000

### First: hierarhical centering:
mod.hc <- textConnection(mod.string.hc)
jm.hc =jags.model(mod.hc,
                data=data.soln,
                n.chains=n.chains,
                n.adapt=n.adapt)
# If using folded-Cauchy prior, may need to provide jags.model with initials for related parameters, such as:
# inits = list(list(alpha = 0.5, tau.temp = 1), list(alpha = 1.5, tau.temp = 0.1), 
#             list(alpha = 2, tau.temp = 1.5))
coda.hc = coda.samples(jm.hc,variable.names=c("deviance","b0","b1","sig","sig.eps", "mean.eps",
                                              "eps"),
                     n.iter=n.iter)

##########################################################
### Next: sum-to-zero constraint

mod.sz <- textConnection(mod.string.sz)
jm.sz =jags.model(mod.sz,
                  data=data.soln,
                  n.chains=n.chains,
                  n.adapt=n.adapt)
# If using folded-Cauchy prior, may need to provide initials, as indicated above.
coda.sz = coda.samples(jm.sz,variable.names=c("deviance","b0","b1","sig","sig.eps", "mean.eps",
                                              "eps"),
                       n.iter=n.iter)


##########################################################
### Final: post-sweeping of random effects

mod.ps <- textConnection(mod.string.ps)
jm.ps =jags.model(mod.ps,
                  data=data.soln,
                  n.chains=n.chains,
                  n.adapt=n.adapt)
# If using folded-Cauchy prior, may need to provide initials, as indicated above.
coda.ps = coda.samples(jm.ps,variable.names=c("deviance","b0.star","b1","sig","sig.eps",
                                              "eps.star"),
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
## to update the MCMC iterations (e.g., see code in Code A2)

stats.hc <- summary(window(coda.hc), start=500)
sum_hc <- cbind.data.frame(mean=stats.hc$statistics[,1],
                             lower2.5=stats.hc$quantiles[,1],
                             upper97.5=stats.hc$quantiles[,5])
sum_hc<-round(sum_hc,4)

stats.sz <- summary(window(coda.sz), start=500)
sum_sz <- cbind.data.frame(mean=stats.sz$statistics[,1],
                           lower2.5=stats.sz$quantiles[,1],
                           upper97.5=stats.sz$quantiles[,5])
sum_sz<-round(sum_sz,4)

stats.ps <- summary(window(coda.ps), start=500)
sum_ps <- cbind.data.frame(mean=stats.ps$statistics[,1],
                           lower2.5=stats.ps$quantiles[,1],
                           upper97.5=stats.ps$quantiles[,5])
sum_ps<-round(sum_ps,4)