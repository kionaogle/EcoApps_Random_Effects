## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
fig.path='sweavefigs/',
prompt=TRUE,
comment=NA,
tidy=FALSE,
size="normalsize",
cache=FALSE)


## ------------------------------------------------------------------------
load("../Data S1.RData")
ls()
s2.df<- cbind.data.frame(y1=as.vector(y1), 
                         y2=as.vector(y2),
                         y3=as.vector(y3),
                         x=as.vector(x),
                         unit=rep(1:10,rep(10,10)))
s2y1.gD<- nlme::groupedData(y1 ~ x | unit,
                            data=s2.df,
                            labels=list(x="x", y1="y1"))
s2y2.gD<- nlme::groupedData(y2 ~ x | unit,
                            data=s2.df,
                            labels=list(x="x", y2="y2"))
s2y3.gD<- nlme::groupedData(y3 ~ x | unit,
                            data=s2.df,
                            labels=list(x="x", y3="y3"))


## ------------------------------------------------------------------------
writeLines(readLines("S2.v2.stan"))


## ------------------------------------------------------------------------
library(rstan, quietly=TRUE)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
S2.v2.stanc<- stanc(file="S2.v2.stan")


## ------------------------------------------------------------------------
S2.v2.stanmod<- stan_model(stanc_ret=S2.v2.stanc)


## ------------------------------------------------------------------------
## y ## to come later, with particular scenarios
X<- model.matrix(y1 ~ x, random = ~1 | unit, data=s2y1.gD)
Z<- model.matrix(nlme::reStruct(list(unit=~1), data=s2y1.gD), 
                 data=s2y1.gD)
(N<- dim(X)[1])
(p<- dim(X)[2])
(ni<- table(as.numeric(as.character((s2y1.gD$unit)))))
(m<- length(ni))
sum(ni) == N ## check
(q<- dim(Z)[2])

mu_b<- rep(0,p);
Sig_b<- 10e6 * diag(p)

tau_eps_a <- tau_eps_b <- rep(0.1, q)
mu_eps <- rep(0, q)

tau_a<- tau_b<- 0.1

nchains<- 3


## ------------------------------------------------------------------------
S2.v2.1.stan.dat<- list(
    m = m,
    N = N,
    ni = ni,
    p = p,
    q = q,
    y = s2y1.gD$y1, ## sig.eps = 0.1 * sig
    X = X,
    Z = Z,
    mu_b = mu_b,
    Sig_b = Sig_b,
    tau_eps_a = as.array(tau_eps_a,1),
    tau_eps_b = as.array(tau_eps_b,1),
    mu_eps = as.array(mu_eps,1),
    tau_a = tau_a,
    tau_b = tau_b)


## ------------------------------------------------------------------------
s2y1.mod<- nlme::lme(y1 ~ x, 
                     random = ~1 | unit,
                     data=s2y1.gD)


## ------------------------------------------------------------------------
set.seed(20500 + 5150 + 24601 + 1)

## b[p]
b<- mvtnorm::rmvnorm(n=nchains,
                     nlme::fixef(s2y1.mod),
                     3*vcov(s2y1.mod))

## eps[q,m]
mu_eps<- t(nlme::ranef(s2y1.mod))
eps<- list(matrix(NA,q,m),matrix(NA,q,m),matrix(NA,q,m))
Sig_eps<- unclass(nlme::getVarCov(s2y1.mod, type="random.effects"))
for(c in 1:nchains) {
    for(i in 1:m) {
        eps[[c]][,i]<- mvtnorm::rmvnorm(n=q,mu_eps[,i],3*Sig_eps)
    }
}

## tau_eps[q] = 1/sig.eps^2
tau_eps<- 1000

## tau
tau<- 1.0

S2.v2.1.stan.inits<- list(
    list(b=b[1,],
         eps=eps[[1]],
         tau_eps=as.array(tau_eps,1),
         tau=tau),
    list(b=b[2,],
         eps=eps[[2]],
         tau_eps=as.array(tau_eps,1),
         tau=tau),
    list(b=b[3,],
         eps=eps[[3]],
         tau_eps=as.array(tau_eps,1),
         tau=tau))


## ------------------------------------------------------------------------
S2.v2.1.stanfit<- sampling(S2.v2.stanmod,
                           data=S2.v2.1.stan.dat,
                           init=S2.v2.1.stan.inits,
                           chains=nchains,
                           iter=10000,
                           warmup=5000,
                           seed=c(20500+5150+24601+1),
                           refresh=1000,
                           sample_file="S2.v2.1.chain")


## ------------------------------------------------------------------------
summary(S2.v2.1.stanfit)$summary


## ------------------------------------------------------------------------
traceplot(S2.v2.1.stanfit,
          pars=names(S2.v2.1.stanfit),
          inc_warmup=TRUE)


## ------------------------------------------------------------------------
S2.v2.2.stan.dat<- S2.v2.1.stan.dat
S2.v2.2.stan.dat$y<- s2y2.gD$y2 ## sig.eps = 1 * sig


## ------------------------------------------------------------------------
s2y2.mod<- nlme::lme(y2 ~ x, 
                     random = ~1 | unit,
                     data=s2y2.gD)


## ------------------------------------------------------------------------
set.seed(20500 + 5150 + 24601 + 2)

## b[p]
b<- mvtnorm::rmvnorm(n=nchains,
                     nlme::fixef(s2y2.mod),
                     3*vcov(s2y2.mod))

## eps[q,m]
mu_eps<- t(nlme::ranef(s2y2.mod))
eps<- list(matrix(NA,q,m),matrix(NA,q,m),matrix(NA,q,m))
Sig_eps<- unclass(nlme::getVarCov(s2y2.mod, type="random.effects"))
for(c in 1:nchains) {
    for(i in 1:m) {
        eps[[c]][,i]<- mvtnorm::rmvnorm(n=q,mu_eps[,i],3*Sig_eps)
    }
}

## tau_eps[q] = 1/sig.eps^2
tau_eps<- 1.0

## tau
tau<- 1.0

S2.v2.2.stan.inits<- list(
    list(b=b[1,],
         eps=eps[[1]],
         tau_eps=as.array(tau_eps,1),
         tau=tau),
    list(b=b[2,],
         eps=eps[[2]],
         tau_eps=as.array(tau_eps,1),
         tau=tau),
    list(b=b[3,],
         eps=eps[[3]],
         tau_eps=as.array(tau_eps,1),
         tau=tau))


## ------------------------------------------------------------------------
S2.v2.2.stanfit<- sampling(S2.v2.stanmod,
                           data=S2.v2.2.stan.dat,
                           init=S2.v2.2.stan.inits,
                           chains=nchains,
                           iter=10000,
                           warmup=5000,
                           seed=c(20500+5150+24601+2),
                           refresh=1000,
                           sample_file="S2.v2.2.chain")


## ------------------------------------------------------------------------
summary(S2.v2.2.stanfit)$summary


## ------------------------------------------------------------------------
traceplot(S2.v2.2.stanfit,
          pars=names(S2.v2.2.stanfit),
          inc_warmup=TRUE)


## ------------------------------------------------------------------------
S2.v2.3.stan.dat<- S2.v2.1.stan.dat
S2.v2.3.stan.dat$y<- s2y2.gD$y3 ## sig.eps = 10 * sig


## ------------------------------------------------------------------------
s2y3.mod<- nlme::lme(y3 ~ x, 
                     random = ~1 | unit,
                     data=s2y3.gD)


## ------------------------------------------------------------------------
set.seed(20500 + 5150 + 24601 + 3)

## b[p]
b<- mvtnorm::rmvnorm(n=nchains,
                     nlme::fixef(s2y3.mod),
                     3*vcov(s2y3.mod))

## eps[q,m]
mu_eps<- t(nlme::ranef(s2y3.mod))
eps<- list(matrix(NA,q,m),matrix(NA,q,m),matrix(NA,q,m))
Sig_eps<- unclass(nlme::getVarCov(s2y3.mod, type="random.effects"))
for(c in 1:nchains) {
    for(i in 1:m) {
        eps[[c]][,i]<- mvtnorm::rmvnorm(n=q,mu_eps[,i],3*Sig_eps)
    }
}

## tau_eps[q] = 1/sig.eps^2
tau_eps<- 1/100

## tau
tau<- 1.0

S2.v2.3.stan.inits<- list(
    list(b=b[1,],
         eps=eps[[1]],
         tau_eps=as.array(tau_eps,1),
         tau=tau),
    list(b=b[2,],
         eps=eps[[2]],
         tau_eps=as.array(tau_eps,1),
         tau=tau),
    list(b=b[3,],
         eps=eps[[3]],
         tau_eps=as.array(tau_eps,1),
         tau=tau))


## ------------------------------------------------------------------------
S2.v2.3.stanfit<- sampling(S2.v2.stanmod,
                           data=S2.v2.3.stan.dat,
                           init=S2.v2.3.stan.inits,
                           chains=nchains,
                           iter=10000,
                           warmup=5000,
                           seed=c(20500+5150+24601+3),
                           refresh=1000,
                           sample_file="S2.v2.3.chain")


## ------------------------------------------------------------------------
summary(S2.v2.3.stanfit)$summary


## ------------------------------------------------------------------------
traceplot(S2.v2.3.stanfit,
          pars=names(S2.v2.3.stanfit),
          inc_warmup=TRUE)


## ------------------------------------------------------------------------
S2.v2.1.stanfit.coda<- As.mcmc.list(S2.v2.1.stanfit)
S2.v2.2.stanfit.coda<- As.mcmc.list(S2.v2.2.stanfit)
S2.v2.3.stanfit.coda<- As.mcmc.list(S2.v2.3.stanfit)
save(list=c("S2.v2.1.stanfit.coda",
            "S2.v2.2.stanfit.coda",
            "S2.v2.3.stanfit.coda"),
     file="S2.v2.stanfit.coda.RData")


## ------------------------------------------------------------------------
## clean up
detach(package:rstan)
detach(package:StanHeaders)
detach(package:ggplot2)
rm(list=ls()) ## careful

