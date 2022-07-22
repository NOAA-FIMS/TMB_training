library(TMB)
# 1. Specify model in C++

# 2. Compile the model and link to R
compile("src/linreg.cpp")
dyn.load(dynlib("src/linreg"))

# Simulate data
set.seed(123)
data <- list(y = rnorm(10) + 1:10, X=cbind(rep(1,10),1:10))

# 3. Construct the C++ objective function with derivatives
parameters <- list(beta = c(0,0), lnSigma=0)
obj <- MakeADFun(data, parameters, DLL="linreg")

# Pass initial values, function, and gradient functions
# to R minimizer, nlminb
opt <- nlminb(obj$par, obj$fn, obj$gr)
# nlminb is what minimizes the function and finds the MLEs