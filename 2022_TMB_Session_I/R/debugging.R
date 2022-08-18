library(TMB)

#Debug using RStudio:
TMB:::setupRStudio()
compile('2022_TMB_Session_I/src/debug1.cpp')

#TMB 1.8.0: disable Eigen warnings
compile('2022_TMB_Session_I/src/debug1.cpp', eigen.disable.warnings = TRUE)


#------------------------------------------------------------------------
# TMB-R linkage errors
#------------------------------------------------------------------------
library(TMB)
#Simulate data 
set.seed(123)
sig <- 2
beta <- c(0.5,2)
X <- cbind(rep(1,50), 1:50)
mu <- X%*%beta
y = rnorm(50, mu, sig)


# #Code will crash RStudio
# Data <- list(y = c(-2.1, 3.3, 4.2),
#              X = rbind(c(1,1),c(4,5))
# )
# 
# Pars <- list(beta = 0)
# obj <- MakeADFun(data = Data, 
#                  parameters = Pars,  
#                  DLL = "linReg")

#Debug using gdbsource
#hangs on Windows
#gdbsource("2022_TMB_Session_I/R/debug_gdbsource.R")
#works on Windows with interactive mode
gdbsource("2022_TMB_Session_I/R/debug_gdbsource.R", 
          interactive = TRUE)

compile('2022_TMB_Session_I/src/linReg.cpp')
dyn.load(dynlib('2022_TMB_Session_I/src/linReg'))
Data <- list(y = c(-2.1, 3.3, 4.2),
             X = cbind(c(1,1,1),c(4,5,6))
)
Pars <- list(beta = c(0,0), lnSigma = 0)
obj <- MakeADFun(data = Data, 
                 parameters = Pars,  
                 DLL = "linReg")

dyn.unload(dynlib('2022_TMB_Session_I/src/linReg'))
#---------------------------------------------------------------------
# Run-time errors
#--------------------------------------------------------------------

compile("2022_TMB_Session_I/src/debug2.cpp")
dyn.load(dynlib("2022_TMB_Session_I/src/debug2"))
Data <- list(y = c(2.1, 3.3, 4.2),
             X = cbind(c(1,1,1),c(4,5,2)))
Pars <- list(beta = c(0,0), lnSigma = 0)
obj <- MakeADFun(data = Data, 
                 parameters = Pars, 
                 DLL = "debug2")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$convergence

dyn.unload(dynlib("2022_TMB_Session_I/src/debug2"))

#----------------------------------------------------------------------
# Convergence errors
#----------------------------------------------------------------------

compile("2022_TMB_Session_I/src/debug3.cpp")
dyn.load(dynlib("2022_TMB_Session_I/src/debug3"))

set.seed(123)
sig <- 2
n <- 100
beta <- c(0.5,2)
X <- cbind(rep(1,n), 1:n)
mu <- X%*%beta
y = rnorm(n, mu, sig)
plot(X[,2], y)

Data <- list(y = y, X = X)
Pars <- list(beta = c(0,0), lnSigma = rep(0,n))
obj <- MakeADFun(data = Data, 
                 parameters = Pars, 
                 DLL = "debug3")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$message
opt$iterations
opt$evaluations
opt <- nlminb(obj$par, obj$fn, obj$gr, 
              control = list(iter.max = 1000, eval.max = 1000))
opt$message
sdr <- sdreport(obj)
sdr$pdHess
length(opt$par)
summary(sdr, "report")
dyn.unload(dynlib("2022_TMB_Session_I/src/debug3"))

#------------------------------------------------------------------
# Model Validation
#------------------------------------------------------------------

dyn.load(dynlib("2022_TMB_Session_I/src/linReg"))

set.seed(123)
sig <- 1
n <- 100
beta <- c(0.5,2)
X <- cbind(rep(1,n), 1:n)
mu <- X%*%beta
y = rnorm(n, mu, sig)

plot(X[,2], y)

#correctly specified model
Data0 <- list(y = y, 
              X = cbind(rep(1,n), 1:n))
Pars0 <- list(beta = c(0,0), lnSigma = 0)
obj0 <- MakeADFun(data = Data0, 
                  parameters = Pars0, 
                  DLL = "linReg")
opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr, 
               control = list(iter.max = 1000, eval.max = 1000))
opt0$convergence
obj0$gr(opt0$par)
sdr0 <- sdreport(obj0)
sdr0$pdHess
sdr0
osa0 <- oneStepPredict(obj0, observation.name = "y", method = "fullGaussian")
qqnorm(osa0$residual);abline(0,1)
#calculate AIC
2 * (opt0$objective + length(opt0$par))


#mis-specified model
Data1 <- list(y = y, 
              X = as.matrix(rep(1,n)))
Pars1 <- list(beta = 0, lnSigma = 0)
obj1 <- MakeADFun(data = Data1, 
                  parameters = Pars1, 
                  DLL = "linReg")
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr, 
               control = list(iter.max = 1000, eval.max = 1000))
opt1$convergence
obj1$gr(opt1$par)
sdr1 <- sdreport(obj1)
sdr1$pdHess
sdr1
osa1 <- oneStepPredict(obj1, observation.name = "y", method = "fullGaussian")
qqnorm(osa1$residual);abline(0,1)
#calculate AIC
2 * (opt1$objective + length(opt1$par))
plot(X[,2], osa1$residual)
