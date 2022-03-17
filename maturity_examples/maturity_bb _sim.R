library(TMB)
#precompile()

dyn.unload(dynlib("maturity_bb_sim"))
compile("maturity_bb_sim.cpp")
dyn.load(dynlib("maturity_bb_sim"))

a50 = 5
slope = 2
beta = c(-slope*a50,slope)
phi = 0.1 #as -> Inf betabin -> binomial
ages = t(matrix(1:20, 20, 40))
mat = 1/(1+exp(-(beta[1] + beta[2]*ages)))
set.seed(123)
N = matrix(sample(50:100, length(ages),replace=TRUE), nrow = 20)
mat_bb = rbeta(length(N),mat*phi, (1-mat)*phi)
Y_bb = matrix(rbinom(length(N), N, mat_bb), nrow = 20)


input = list(data=list(),par=list())
input$par$beta = c(0,0)
input$par$log_phi = 0
input$data$Y = c(Y_bb)
input$data$N = c(N)
input$data$age_obs = c(ages)
input$data$max_age = max(ages)

mod = MakeADFun(input$data, input$par, DLL = "maturity_bb_sim")

opt = nlminb(mod$par, mod$fn, mod$gr)
mod$rep = mod$report()
mod$sdrep = sdreport(mod)
summary(mod$sdrep)

set.seed(456)
newdat = mod$simulate(complete = TRUE)
newinput = input
newinput$data = newdat

newmod = MakeADFun(newinput$data, newinput$par, DLL = "maturity_bb_sim")

newopt = nlminb(newmod$par, newmod$fn, newmod$gr)

