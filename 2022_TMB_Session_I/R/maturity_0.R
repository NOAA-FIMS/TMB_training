library(TMB)

compile("2022_TMB_Session_I/src/maturity_0.cpp")
dyn.load(dynlib("2022_TMB_Session_I/src/maturity_0"))

a50 = 5
slope = 2
beta = c(-slope*a50,slope)
ages = t(matrix(1:20, 20, 40))
mat = 1/(1+exp(-(beta[1] + beta[2]*ages)))
set.seed(123)
N = matrix(sample(5:10, length(ages),replace=TRUE), 40, 20)
Y = matrix(rbinom(length(N), N, mat), 40, 20)


input = list(data=list(),par=list())
input$par$beta = c(0,0)
input$data$Y = c(Y)
input$data$N = c(N)
input$data$age_obs = c(ages)
input$data$max_age = max(ages)

mod = MakeADFun(input$data, input$par, DLL = "maturity_0")

opt = nlminb(mod$par, mod$fn, mod$gr)
mod$rep = mod$report()
mod$sdrep = sdreport(mod)
summary(mod$sdrep, "fixed")
beta

logit_mat = as.list(mod$sdrep, "Estimate", report = TRUE)$logit_mat_at_age
logit_mat.se =as.list(mod$sdrep, "Std. Error", report = TRUE)$logit_mat_at_age
logit_mat.ci =matrix(logit_mat, 20, 2) + cbind(-qnorm(0.975)*logit_mat.se, qnorm(0.975)*logit_mat.se)
plot(ages[1,], 1/(1+exp(-logit_mat)), type = 'l', ylim = c(0,1))
polygon(c(ages[1,],rev(ages[1,])), 1/(1 + exp(-c(logit_mat.ci[,1],rev(logit_mat.ci[,2])))), lty = 2)
  
## Map parameters
input$map$beta = factor(c(1,NA)) 
mod = MakeADFun(input$data, input$par, map = input$map,  DLL = "maturity_0")
opt = nlminb(mod$par, mod$fn, mod$gr)
opt$par

#intercept = slope
input$map$beta = factor(c(1,1)) 
mod = MakeADFun(input$data, input$par, map = input$map,  DLL = "maturity_0")
opt = nlminb(mod$par, mod$fn, mod$gr)
opt$par
