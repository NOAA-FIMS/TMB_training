# Building a TMB model in R ================================

#-----------------------------------------------------------
# All C++ models, including TMB, need to be compiled
# In TMB, this creates a dynamic link library (.dll/.so)
# The .dll/.so file needs to be loaded for R to access
#-----------------------------------------------------------

library(TMB)
compile('2022_TMB_Session_I/src/linReg.cpp')
dyn.load(dynlib('2022_TMB_Session_I/src/linReg'))


#------------------------------------------------------------------------
# TMB additionally needs to build a computational graph
# - This sets up derivative chain used to determine the gradiant functions
# - TMB has a static graph, so this step happens before optimization
# - The R function, MakeADFun() completes this step
#------------------------------------------------------------------------

#Simulate data 
set.seed(123)
sig <- 2
beta <- c(0.5,2)
X <- cbind(rep(1,50), 1:50)
mu <- X%*%beta
y = rnorm(50, mu, sig)

Data <- list(y = y, 
             X = X)
Pars <- list(beta = c(0,0),lnSigma = 0)
obj <- MakeADFun(data = Data, 
                 parameters = Pars, 
                 DLL = 'linReg')

#--------------------------------------------------------------
# MakeADFun output:
# **par**: initial values
# **fn**: the likelihood function
# **gr**: the gradient function
# **report**: return values wrapped in REPORT()
# **env**: environment with access to all parts of structure
#
# The model is optimized outside of TMB using any R minimizer
# nlminb is a base R minimizer frequently used for TMB models
#--------------------------------------------------------------

# Pass initial values, function, and gradient functions
opt <- nlminb(obj$par, obj$fn, obj$gr)

#---------------------------------------------------------------------
# nlminb output:
# **par**: MLEs
# **objective**: the nll associated with the MLEs
# **convergence**: 0 indicates successful convergence
# **message**: convergence message
# **iterations**: number of iterations
# **evaluations**: number of objective and gradient function evaluations
#----------------------------------------------------------------------
# Report objects back to R
report <- obj$report()
report$sig2; 

## Check model convergence

#Check maximum gradient < 0.0001
obj$gr()
#Check convergence error = 0
opt$convergence
#Check convergence message
opt$message

## Uncertainty 
# calculates standard errors for all parameters and anything within ADREPORT()
sdr <- sdreport(obj)

#uses the Hessian to calculate parameter uncertainty
summary(sdr, "fixed")
#uses the Delta method to calculate undertainty for derived quantities
summary(sdr, "report")

#Check Hessian is positive definite
sdr$pdHess

