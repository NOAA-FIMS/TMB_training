## modified code from: https://github.com/kaskr/adcomp/blob/master/TMB/inst/examples/ar1xar1.R
library(TMB)
#compile model
compile("2022_TMB_Session_II/src/ar1xar1_modular.cpp")

set.seed(123)
n <- 20 ## Size of problem = n*n

## ======================= Simulate separable 2D GMRF 
## - With exponential correlation in both directions
## - phi1 = 1-lag correlation in 1st direction
## - phi2 = 1-lag correlation in 2nd direction
ar1corr <- function(n,phi){
  phi^abs(outer(1:n,1:n,"-"))
}
simgmrf <- function(n1,n2,phi1,phi2){
  u <- matrix(rnorm(n1*n2),n1,n2)
  L1 <- t(chol(ar1corr(n1,phi1)))
  L2 <- t(chol(ar1corr(n2,phi2)))
  x <- L1%*%u         ## phi1 in 1st direction (fastest)
  x <- t(L2%*%t(x))   ## phi2 in 2nd direction
  x
}

## ======================= Simulate data
phi1=exp(-1/(.1*n)) ## Correlation range=10% of grid size first dimension
phi2=exp(-1/(.2*n)) ## Correlation range=20% of grid size second dimension
eta <- simgmrf(n,n,phi1,phi2)
y <- rpois(length(eta),exp(eta))
d <- expand.grid(x=factor(1:n),y=factor(1:n))
d$y <- y

## ======================= Parameterization of phi
f <- function(x) 2/(1 + exp(-2 * x)) - 1
invf <- function(y) -0.5 * log(2/(y + 1) - 1)

## ======================= Fit model
dyn.load(dynlib("2022_TMB_Session_II/src/ar1xar1_modular"))
obj <- MakeADFun(data=list(y=y),
                 parameters=list(
                   eta=matrix(0,n,n),
                   transf_phi1=invf(0.5),
                   transf_phi2=invf(0.5)),
                 random=c("eta"),
                 DLL="ar1xar1_modular")
#runSymbolicAnalysis(obj)
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par
sdr <- sdreport(obj)
summary(sdr, 'report')
phi1;phi2

#test TMB validation
#oneStepGeneric slower than cdf
# osa.gen <- oneStepPredict(obj, 'y', data.term.indicator = 'keep',
#                          method = 'oneStepGeneric', discrete = TRUE,
#                          range = c(0,Inf))
# osa.cdf <- oneStepPredict(obj, 'y', data.term.indicator = 'keep',
#                          method = 'cdf', discrete = TRUE,
#                          range = c(0,Inf))
# save(osa.gen, osa.cdf, file = "2022_TMB_Session_II/R/osa_ar1xar1.RData")
load("2022_TMB_Session_II/R/osa_ar1xar1.RData")
qqnorm(osa.gen$residual);abline(0,1);ks.test(osa.gen$residual, 'pnorm')
qqnorm(osa.cdf$residual);abline(0,1);ks.test(osa.cdf$residual, 'pnorm')

