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
years = matrix(1:40, 20,40)

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

input$map$beta = factor(c(1,NA)) #fix slope = 0 : constant proportion with age
mod = MakeADFun(input$data, input$par, map = input$map, DLL = "maturity_bb")

opt = nlminb(mod$par, mod$fn, mod$gr)
opt$par
mod$sdrep = sdreport(mod)
as.list(mod$sdrep, "Estimate", report = TRUE)$mat_at_age


dyn.unload(dynlib("maturity_bb_reg"))
compile("maturity_bb_reg.cpp")
dyn.load(dynlib("maturity_bb_reg"))

dd = data.frame(age = c(ages), year = c(years))
X = cbind(dd$year<=20, dd$year>20, dd$age)
head(X)
a50_reg = c(5,3)
beta = c(-slope*a50_reg,slope)
mat = 1/(1+exp(-X %*% beta))
set.seed(123)
N = matrix(sample(50:100, length(ages),replace=TRUE), nrow = 20)
mat_bb = rbeta(length(N),mat*phi, (1-mat)*phi)
Y_bb = matrix(rbinom(length(N), N, mat_bb), nrow = 20)

input = list(data=list(),par=list())
input$par$beta = rep(0, NCOL(X))
input$par$log_phi = 0
input$data$Y = c(Y_bb)
input$data$N = c(N)
input$data$X = X
input$data$age_col = NCOL(X)
input$data$max_age = max(ages)

mod = MakeADFun(input$data, input$par, map = input$map, DLL = "maturity_bb_reg")

opt = nlminb(mod$par, mod$fn, mod$gr)
opt$par
mod$sdrep = sdreport(mod)
as.list(mod$sdrep, "Estimate", report = TRUE)$a50

input$map$beta = factor(c(1,1,2))
mod = MakeADFun(input$data, input$par, map = input$map, DLL = "maturity_bb_reg")

opt = nlminb(mod$par, mod$fn, mod$gr)
opt$par
mod$sdrep = sdreport(mod)
as.list(mod$sdrep, "Estimate", report = TRUE)$a50
as.list(mod$sdrep, "Estimate")$beta

