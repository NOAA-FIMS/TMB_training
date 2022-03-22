library(TMB)
#precompile()

dyn.unload(dynlib("maturity_bb"))
compile("maturity_bb.cpp")
dyn.load(dynlib("maturity_bb"))

a50 = 5
slope = 2
beta = c(-slope*a50,slope)
phi = 0.1 #as -> Inf betabin -> binomial
ages = t(matrix(1:20, 20, 40))
mat = 1/(1+exp(-(beta[1] + beta[2]*ages)))
set.seed(123)
N = matrix(sample(50:100, length(ages),replace=TRUE), 40,20)
mat_bb = rbeta(length(N),mat*phi, (1-mat)*phi)
Y_bb = matrix(rbinom(length(N), N, mat_bb), 40,20)


input = list(data=list(),par=list())
input$par$beta = c(0,0)
input$par$log_phi = 0
input$data$Y = c(Y_bb)
input$data$N = c(N)
input$data$age_obs = c(ages)
input$data$max_age = max(ages)

mod = MakeADFun(input$data, input$par, DLL = "maturity_bb")

opt = nlminb(mod$par, mod$fn, mod$gr)
mod$rep = mod$report()
mod$sdrep = sdreport(mod)
summary(mod$sdrep)
logit_mat = as.list(mod$sdrep, "Estimate", report = TRUE)$logit_mat_at_age
logit_mat.se =as.list(mod$sdrep, "Std. Error", report = TRUE)$logit_mat_at_age
logit_mat.ci =matrix(logit_mat, 20, 2) + cbind(-qnorm(0.975)*logit_mat.se, qnorm(0.975)*logit_mat.se)
plot(ages[1,], 1/(1+exp(-logit_mat)), type = 'l', ylim = c(0,1))
polygon(c(ages[1,],rev(ages[1,])), 1/(1 + exp(-c(logit_mat.ci[,1],rev(logit_mat.ci[,2])))), lty = 2)
