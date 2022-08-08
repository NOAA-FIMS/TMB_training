require(TMB)
compile("2022_TMB_Session_I/src/simple_re.cpp")
dyn.load(dynlib("2022_TMB_Session_I/src/simple_re"))

## Test data
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
u <- rnorm(ncol(B)) ## logsdu=0
beta <- rnorm(ncol(X))*100
eps <- rnorm(nrow(B),sd=1) ## logsd0=0
y <- as.numeric( X %*% beta + B %*% u + eps )

#visualize data
plot(X%*%beta, y)
plot(B%*%u, y)

## Fit model
obj <- MakeADFun(data=list(y=y, B=B, X=X),
                 parameters=list(u=u*0, beta=beta*0, lnSDu=1, lnSDy=1),
                 random="u",
                 DLL="simple_re"
                 )
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par
beta
#marginal likelihood (via Laplace approximation)
obj$fn(opt$par)
report <- obj$report()
report$nll1 #random effect likelihood
report$nll2 #data likelihood 
report$nll #joint likelihood

sdr <- sdreport(obj)
summary(sdr, "fixed")
summary(sdr, "random")
sdr$pdHess
