library(TMB)

#Debug using RStudio:
TMB:::setupRStudio()
compile('2022_TMB_Session_I/src/debug1.cpp')

#TMB 1.8.0: disable Eigen warnings
compile('2022_TMB_Session_I/src/debug1.cpp', eigen.disable.warnings = TRUE)

#Debug using gdbsource
#hangs on Windows
gdbsource("2022_TMB_Session_I/R/debug_gdbsource.R")
#works on Windows with interactive mode
gdbsource("2022_TMB_Session_I/R/debug_gdbsource.R", 
                 interactive = TRUE)

dyn.load(dynlib('2022_TMB_Session_I/src/linReg'))


#------------------------------------------------------------------------
# TMB-R linkage errors
#------------------------------------------------------------------------

#Simulate data 
set.seed(123)
sig <- 2
beta <- c(0.5,2)
X <- cbind(rep(1,50), 1:50)
mu <- X%*%beta
y = rnorm(50, mu, sig)


compile('2022_TMB_Session_I/src/linReg.cpp')
Data <- list(y = c(-2.1, 3.3, 4.2),
             X = data.frame(c(1,1),c(4,5))
)
Pars <- list(beta = 0)
obj <- MakeADFun(data = Data, 
                 parameters = Pars,  
                 DLL = "linReg")


#---------------------------------------------------------------------
# Run-time errors
#--------------------------------------------------------------------
Data <- list(y = c(-2.1, 3.3, 4.2),
             X = cbind(c(1,1,1),c(4,5,2))
)
Pars <- list(beta = c(0,0), sigma = 0)
obj <- MakeADFun(data = Data, 
                 parameters = Pars, 
                 DLL = "linReg")
opt <- nlminb(obj$par, obj$fn, obj$gr)

compile("2022_TMB_Session_I/src/debug2.cpp")
dyn.load(dynlib("2022_TMB_Session_I/src/debug2"))
Data <- list(y = c(-2.1, 3.3, 4.2),
             X = cbind(c(1,1,1),c(4,5,2)))
Pars <- list(beta = c(0,0), sigma = 0)
obj <- MakeADFun(data = Data, 
                 parameters = Pars, 
                 DLL = "debug2")
opt <- nlminb(obj$par, obj$fn, obj$gr)


