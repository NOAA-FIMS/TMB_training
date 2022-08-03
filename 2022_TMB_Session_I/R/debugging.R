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
Data <- list(y = c(-2.1, 3.3, 4.2),
             X = data.frame(c(1,1,1),c(4,5,6))
)
Pars <- list(beta = c(0,0))
obj <- MakeADFun(data = Data, 
                 parameters = Pars,  
                 DLL = "linReg")


#---------------------------------------------------------------------
# Run-time errors
#--------------------------------------------------------------------

compile("2022_TMB_Session_I/src/debug2.cpp")
dyn.load(dynlib("2022_TMB_Session_I/src/debug2"))
Data <- list(y = c(-2.1, 3.3, 4.2),
             X = cbind(c(1,1,1),c(4,5,2)))
Pars <- list(beta = c(0,0), sigma = 0)
obj <- MakeADFun(data = Data, 
                 parameters = Pars, 
                 DLL = "debug2")
opt <- nlminb(obj$par, obj$fn, obj$gr)


