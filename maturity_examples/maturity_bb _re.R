library(TMB)
#precompile()

dyn.unload(dynlib("maturity_bb_re"))
compile("maturity_bb_re.cpp")
dyn.load(dynlib("maturity_bb_re"))

a50 = 5
slope = 0.5
beta = c(-slope*a50,slope)
phi = 0.1 #as -> Inf betabin -> binomial
ages = t(matrix(1:20, 20, 40))
years = matrix(1:40, 40,20)
mat = 1/(1+exp(-(beta[1] + beta[2]*ages)))
set.seed(123)
N = matrix(sample(50:100, length(ages),replace=TRUE), 40,20)
mat_bb = rbeta(length(N),mat*phi, (1-mat)*phi)
Y_bb = matrix(rbinom(length(N), N, mat_bb), nrow = 40,20)

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


opt = nlminb(mod$par, mod$fn, mod$gr)
mod$rep = mod$report()
matrix(mod$rep$mat, 40)
mod$sdrep = sdreport(mod)
summary(mod$sdrep)
mod$env$parList()$re
input$data$re

set.seed(456)
newdat = mod$simulate(complete = TRUE)
newinput = input
newinput$data = newdat

newmod = MakeADFun(newinput$data, newinput$par, DLL = "maturity_bb_sim")

newopt = nlminb(newmod$par, newmod$fn, newmod$gr)

