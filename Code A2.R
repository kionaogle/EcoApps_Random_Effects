### Script for running the JAGS (or OpenBUGS) code for implementing the
### non-identifiable random effects regression associated with 
### Equation 5 in the main text.

library(rjags)

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


###############################################################
### JAGS (or OpenBUGS) model code for implementing the 
### original, non-identifiable model, before reparameterizing
### or implementation of computational or coding tricks:
mod.string <- "
model{
  for(j in 1:10){
    for(r in 1:10){
      # Likelihood of (stochastic) data:
      y[r,j] ~ dnorm(mu[r,j],tau)
      # Mean model for a random effects regression:
      mu[r,j] <- b0 + b1*x[r,j] + eps[j]
    }
    # Zero-centered hierarchical prior for random effects:
    eps[j] ~ dnorm(0,tau.eps)
  }
  # Relatively non-informative priors for root nodes:
  b0 ~ dnorm(0,1E-6)
  b1 ~ dnorm(0,1E-6)
  tau ~ dgamma(0.1, 0.1)
  sig <- 1/sqrt(tau)
  
  # Conjugate gamma prior for precision:
  tau.eps ~ dgamma(0.1, 0.1)
  sig.eps <- 1/sqrt(tau.eps)
  
  # Folded Cauchy(0,A^2,df=1) prior for standard deviation:
  # alpha ~ dnorm(0,1)
  # tau.temp ~ dgamma(a,b)
  # sig.temp <- 1/sqrt(tau.temp)
  # sig.eps <- abs(alpha)*sig.temp
  # tau.eps <- pow(sig.eps,-2)

  # Compute mean of random effects (to illustrate non-identifiability)
  mean.eps <- mean(eps[])
}
"

#####################################################
### Fit the model to the three different datasets
#####################################################

load.module("dic")


### First model, sig.eps = sig/10
mod <- textConnection(mod.string)
n.adapt = 1000
n.chains = 3
n.iter = 10000
jm1 =jags.model(mod,
                 data=data1,
                 n.chains=n.chains,
                 n.adapt=n.adapt)
# If using folded-Cauchy prior, may need to provide jags.model with initials for related parameters, such as:
# inits = list(list(alpha = 0.5, tau.temp = 1), list(alpha = 1.5, tau.temp = 0.1), 
#             list(alpha = 2, tau.temp = 1.5))
coda1 = coda.samples(jm1,variable.names=c("deviance","b0","b1","sig","sig.eps", "mean.eps","eps"),
                   n.iter=n.iter)

### Second model, sig.eps = sig
mod <- textConnection(mod.string)
jm2 =jags.model(mod,
                data=data2,
                n.chains=n.chains,
                n.adapt=n.adapt)
# If using folded-Cauchy prior, may need to provide initials, as indicated above.
coda2 = coda.samples(jm2,variable.names=c("deviance","b0","b1","sig","sig.eps", "mean.eps","eps"),
                     n.iter=n.iter)

### Third model, sig.eps = 10sig
mod <- textConnection(mod.string)
jm3 =jags.model(mod,
                data=data3,
                n.chains=n.chains,
                n.adapt=n.adapt)
# If using folded-Cauchy prior, may need to provide initials, as indicated above.
coda3 = coda.samples(jm3,variable.names=c("deviance","b0","b1","sig","sig.eps", "mean.eps","eps"),
                     n.iter=n.iter)

# Compute Raftery diagnostic to determine number of MCMC iterations required:
raft1<-raftery.diag(coda1)
raft2<-raftery.diag(coda2)
raft3<-raftery.diag(coda3)


############################################################
# Update models according to number of samples required (based on Raftery):
n.iter1 = 190000/3
coda1up = coda.samples(jm1,variable.names=c("deviance","b0","b1","sig","sig.eps", "mean.eps", "eps"),
                     n.iter=n.iter1)
n.iter2 = 80000/3
coda2up = coda.samples(jm2,variable.names=c("deviance","b0","b1","sig","sig.eps", "mean.eps","eps"),
                       n.iter=n.iter2)
n.iter3 = 400000/3
coda3up = coda.samples(jm3,variable.names=c("deviance","b0","b1","sig","sig.eps", "mean.eps","eps"),
                       n.iter=n.iter3)

# MCMC plots
library(mcmcplots)
mcmcplot(window(coda1up,thin=100))
mcmcplot(window(coda2up,thin=40))
mcmcplot(window(coda3up,thin=100))

# Compute posterior summary statistics
stats1 <- summary(coda1up)
sum_tab1 <- cbind.data.frame(mean=stats1$statistics[,1],
                           lower2.5=stats1$quantiles[,1],
                           upper97.5=stats1$quantiles[,5])
sum_tab1<-round(sum_tab1,4)

stats2 <- summary(coda2up)
sum_tab2 <- cbind.data.frame(mean=stats2$statistics[,1],
                             lower2.5=stats2$quantiles[,1],
                             upper97.5=stats2$quantiles[,5])
sum_tab2 <- round(sum_tab2,4)

stats3 <- summary(coda3up)
sum_tab3 <- cbind.data.frame(mean=stats3$statistics[,1],
                             lower2.5=stats3$quantiles[,1],
                             upper97.5=stats3$quantiles[,5])
sum_tab3 <- round(sum_tab3,4)