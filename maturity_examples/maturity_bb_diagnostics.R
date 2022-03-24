library(TMB)
#precompile()

#dyn.unload(dynlib("maturity_examples/maturity_bb"))
compile("maturity_examples/maturity_bb.cpp")
dyn.load(dynlib("maturity_examples/maturity_bb"))

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

# Runtime errors
input$data$N[100] = -1
mod = MakeADFun(input$data, input$par, DLL = "maturity_bb")
opt = nlminb(mod$par, mod$fn, mod$gr)

dyn.unload(dynlib("maturity_examples/maturity_bb"))
compile("maturity_examples/maturity_bb.cpp")
dyn.load(dynlib("maturity_examples/maturity_bb"))

input = list(data=list(),par=list())
input$par$beta = c(0,0)
input$par$log_phi = 0
input$data$Y = c(Y_bb)
input$data$N = c(N)
input$data$age_obs = c(ages)
input$data$max_age = max(ages)

mod = MakeADFun(input$data, input$par, DLL = "maturity_bb")

opt = nlminb(mod$par, mod$fn, mod$gr)

# Convergence Errors
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

input$par$beta = c(0,0,0)
mod = MakeADFun(input$data, input$par, DLL = "maturity_bb")

opt = nlminb(mod$par, mod$fn, mod$gr)
opt$par
opt$message
mod$gr()
sdr <- sdreport(mod)
sdr
sdr$pdHess

#Validation
compile('src/simple_re.cpp')
dyn.load(dynlib('src/simple_re'))

## Simulate data
set.seed(123)
yr <- rep(1900:2010,each=2)
year <- factor(yr)
quarter <- factor(rep(1:4,length.out=length(year)))
period <- factor((yr > mean(yr))+1)
## Random year+quarter effect, fixed period effect:
B <- model.matrix(~year+quarter-1)
X <- model.matrix(~period-1)
B <- as(B,"dgTMatrix")
X <- as(X,"dgTMatrix")
u <- rnorm(ncol(B), sd=3) ## logsdu=.69
beta <- rnorm(ncol(X))*100
eps <- rnorm(nrow(B),sd=1) ## logsd0=0
y <- as.numeric( X %*% beta + B %*% u + eps )

## Fit model
obj <- MakeADFun(data=list(y=y, B=B, X=X),
                 parameters=list(u=u*0, beta=beta*0, lnSDu=1, lnSDy=1),
                 random="u",
                 DLL="simple_re",
                 silent=TRUE
)
opt <- nlminb(obj$par, obj$fn, obj$gr)
osa <- oneStepPredict(obj, 
                      observation = 'y',
                      data.vector.indicator = 'keep',
                      method = 'fullGaussian')
qqnorm(osa$residual);qqline(osa$residual)

#Mis-specify the model by fitting without the re term
obj <- MakeADFun(data=list(y=y, B=B, X=X),
                 parameters=list(u=u*0, beta=beta*0, lnSDu=1, lnSDy=1),
                 map = list(u = rep(factor(NA),length(u)), lnSDu = factor(NA)),
                 DLL="simple_re",
                 silent=TRUE
)
opt <- nlminb(obj$par, obj$fn, obj$gr)
osa <- oneStepPredict(obj, 
                      observation = 'y',
                      method = 'fullGaussian')
qqnorm(osa$residual);qqline(osa$residual)

#Other methods:
#oneStepGaussian
osa <- oneStepPredict(obj, 
                      observation = 'y',
                      data.term.indicator = 'keep',
                      method = 'oneStepGaussian')
qqnorm(osa$residual);qqline(osa$residual)
#oneStepGeneric
osa <- oneStepPredict(obj, 
                      observation = 'y',
                      data.term.indicator = 'keep',
                      method = 'oneStepGeneric')
qqnorm(osa$residual);qqline(osa$residual)
