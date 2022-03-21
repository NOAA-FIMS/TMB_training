library(TMB)
#precompile()
setwd("maturity_examples")
dyn.unload(dynlib("maturity_bb_re"))
compile("maturity_bb_re.cpp")
dyn.load(dynlib("maturity_bb_re"))

a50 = 5
slope = 0.5
beta = c(-slope*a50,slope)
phi = 1000 #as -> Inf betabin -> binomial
#phi = 0.1 #as -> Inf betabin -> binomial
ages = t(matrix(1:20, 20, 40))
years = matrix(1:40, 40,20)
mat = 1/(1+exp(-(beta[1] + beta[2]*ages)))
set.seed(123)
N = matrix(sample(50:100, length(ages),replace=TRUE), 40,20)
mat_bb = rbeta(length(N),mat*phi, (1-mat)*phi)
Y_bb = matrix(rbinom(length(N), N, mat_bb), nrow = 40,20)

#mat_bb = rbeta(matrix(100,40,20),mat*phi, (1-mat)*phi)
#matrix(mat_bb,40,20)

gen.logit <- function(x, low, upp) return(log((x-low)/(upp-x)))

input = list(data=list(),par=list())

input$par$beta = beta
input$par$log_phi = log(phi)
input$par$AR_pars = c(log(1), gen.logit(0.5,-1,1))
input$par$re = rep(0, NROW(Y_bb))

input$data$Y = c(Y_bb)
input$data$N = c(N)
input$data$re_ind = c(years)-1
input$data$age_obs = c(ages)
input$data$max_age = max(ages)


mod = MakeADFun(input$data, input$par, random = "re", DLL = "maturity_bb_re")
input$data = mod$simulate(complete=TRUE)
mod = MakeADFun(input$data, input$par, random = "re", DLL = "maturity_bb_re")


opt = nlminb(mod$par, mod$fn, mod$gr)
mod$rep = mod$report()
matrix(mod$rep$mat, 40)
mod$sdrep = sdreport(mod)
summary(mod$sdrep)
mod$env$parList()$re
input$data$re

#check Laplace Approximation
checkConsistency(mod)

#profile for AR1 sigma parameter (not good)
x = tmbprofile(mod,4)
#x = tmbprofile(mod,4, parm.range = c(-10,4))
plot(x)

#bias correction of standard errors for derived variables
as.list(sdreport(mod,bias.correct=TRUE), what = "Estimate", report=TRUE)$a50
as.list(mod$sdrep, what = "Estimate", report=TRUE)$a50
#need a nonlinear function of the parameters?

#instead use namespace density (multivariate normal NEGATIVE log-likelihoods)
dyn.unload(dynlib("maturity_bb_re_density"))
compile("maturity_bb_re_density.cpp")
dyn.load(dynlib("maturity_bb_re_density"))

mod = MakeADFun(input$data, input$par, random = "re", DLL = "maturity_bb_re_density")

opt2 = nlminb(mod$par, mod$fn, mod$gr)

#same results (more or less)
opt
opt2



dyn.unload(dynlib("maturity_bb_re_osa"))
compile("maturity_bb_re_osa.cpp")
dyn.load(dynlib("maturity_bb_re_osa"))

inputosa = input
#inputosa$par = mod$env$parList()
modosa = MakeADFun(inputosa$data, inputosa$par, random = "re", DLL = "maturity_bb_re_osa")
modosa$fn()
opt3 = nlminb(modosa$par, modosa$fn, modosa$gr)
c(opt$obj,opt2$obj,opt3$obj)

#takes a long time.
residuals = oneStepPredict(modosa, observation.name = "Y", data.term.indicator="keep", method = "oneStepGeneric", discrete=TRUE)
plot(input$data$age_obs, residuals$residual)
qqnorm(residuals$residual)
qqline(residuals$residual)

