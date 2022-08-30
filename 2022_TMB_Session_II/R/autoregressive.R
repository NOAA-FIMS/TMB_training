library(TMB)
#devtools::install_github("eddelbuettel/rbenchmark")
library(rbenchmark)

###################
# Visualize autocorrelation
###################

Rho <- 0.8
u <- 0:10
Corr <- Rho^u

plot( x=u, y=Corr, xlab="distance", ylab="Correlation", ylim=c(0,1), pch = 16 )


###################
# Visualize autocorrelation vs. random-walk
###################

Rho = c(-0.5,0,0.5,0.99)
Marginal_Sigma = 1
Conditional_Sigma = sqrt( 1-Rho^2 ) * Marginal_Sigma

Y = X = 1:1e2

par( mfrow=c(4,1), mgp=c(2,0.5,0), mar=c(2,2,1,0), xaxs="i")
for(i in 1:4){
  Y[1] = 0
  for(s in 2:length(X)) Y[s] = Y[s-1]*Rho[i] + rnorm(1, mean=0, sd=Conditional_Sigma[i])
  plot( x=X, y=Y, xlab="location", ylab="value", ylim=c(-3,3), type="l", main=paste0("Rho = ",Rho[i]) )
}
dev.off()
###################
# Math check on inverse-covariance matrix
###################
Rho <- 0.5
Sigma2 <- 2 ^ 2
x <- 1:10
Dist <- outer(x, x, FUN=function(a,b){abs(a-b)})

# Calculate Q directly
Cov <- Sigma2 / (1-Rho^2) * Rho^Dist
Q <- solve(Cov)
Cov
Q

# Calculate Q analytically
M <- diag(length(x)) * (1+Rho^2)
M[ cbind(1:9,2:10) ] <- -Rho
M[ cbind(2:10,1:9) ] <- -Rho
Q2 = 1/Sigma2 * M
Q2

###################
# Equal time-step autoregressive
###################
TMB::compile( "2022_TMB_Session_II/src/autoregressive.cpp" )
dyn.load( TMB::dynlib("2022_TMB_Session_II/src/autoregressive") )

n <- 100
x <- 1:n
Rho <- 0.8
Sigma2 <- (0.5) ^ 2
beta0 <- 3

# Simulate AR1 process ======================================
set.seed(123)
u <- rep(NA, n)
u[1] <- rnorm(1, mean=0, sd=sqrt(Sigma2))
for(i in 2:n){
  u[i] <- Rho*u[i-1] + rnorm(1, mean=0, sd=sqrt(Sigma2))
} 

# Simulate counts
y <- rpois(n, exp(beta0 + u) )

plot(x, y)

# Compile
Pars <- list( beta0 = 0, 
               ln_sigma2 = 0, 
               logit_rho = 0, 
               u = rnorm(n) )
Data <- list( Options_vec = 0, 
              y = y )

par.out <- rep.out <- he.out <- list()
for(i in 1:5){
  Data$Options_vec <- i-1
  Obj <- TMB::MakeADFun( data=Data, parameters=Pars, random="u", DLL="autoregressive" )
  Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
  par.out[[i]] <- Opt$par
  rep.out[[i]] <- Obj$report()
  he.out[[i]] <-  Obj$env$spHess( random=TRUE )
}

par.est <- data.frame(
  beta0 = c(beta0, sapply(par.out, function(x) x['beta0'])),
  rho = c(Rho, unlist(sapply(rep.out, function(x) x['rho']))),
  sigma2 = c(Sigma2, unlist(sapply(rep.out, function(x) x['sigma2'])))
)

cov.nms <- c("Conditional Independence",
             "Analytic Precision",
             "Built-in GMRF",
             "Covariance and MVNORM",
             "Built-in AR1")
row.names(par.est) <- c("True", cov.nms)
par.est

library(INLA)
for(i in 1:5){
  print(image(he.out[[i]], main = cov.nms[i]))
}

# Benchmark Analysis ================================================

n <- 10
set.seed(123)
Pars$u <- rnorm(n)
bnmk_n10 <- rbenchmark::benchmark(
  "Conditional Independence" = { 
    u <- rep(NA, n)
    u[1] <- rnorm(1, mean=0, sd=sqrt(Sigma2))
    for(i in 2:n){
      u[i] <- Rho*u[i-1] + rnorm(1, mean=0, sd=sqrt(Sigma2))
    } 
    y <- rpois(n, exp(beta0 + u) )
    Data <- list( Options_vec = 0, 
                  y = y )
    Obj <- TMB::MakeADFun( data=Data, parameters=Pars, random="u", DLL="autoregressive" )
    # Optimize
    Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
  },
  "Analytic Precision" = {
    u <- rep(NA, n)
    u[1] <- rnorm(1, mean=0, sd=sqrt(Sigma2))
    for(i in 2:n){
      u[i] <- Rho*u[i-1] + rnorm(1, mean=0, sd=sqrt(Sigma2))
    } 
    y <- rpois(n, exp(beta0 + u) )
    Data <- list( Options_vec = 1, 
                  y = y )
    Obj <- TMB::MakeADFun( data=Data, parameters=Pars, random="u", DLL="autoregressive" )
    # Optimize
    Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
  },
  "Built-in GMRF" = {
    u <- rep(NA, n)
    u[1] <- rnorm(1, mean=0, sd=sqrt(Sigma2))
    for(i in 2:n){
      u[i] <- Rho*u[i-1] + rnorm(1, mean=0, sd=sqrt(Sigma2))
    } 
    y <- rpois(n, exp(beta0 + u) )
    Data <- list( Options_vec = 2, 
                 y = y )
    Obj <- TMB::MakeADFun( data=Data, parameters=Pars, random="u", DLL="autoregressive" )
    # Optimize
    Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
  },
  "Covariance and MVNORM" = {
    u <- rep(NA, n)
    u[1] <- rnorm(1, mean=0, sd=sqrt(Sigma2))
    for(i in 2:n){
      u[i] <- Rho*u[i-1] + rnorm(1, mean=0, sd=sqrt(Sigma2))
    } 
    y <- rpois(n, exp(beta0 + u) )
    Data <- list( Options_vec = 3, 
                 y = y )
    Obj <- TMB::MakeADFun( data=Data, parameters=Pars, random="u", DLL="autoregressive" )
    # Optimize
    Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
  },
  "Built-in AR1" = {
    u <- rep(NA, n)
    u[1] <- rnorm(1, mean=0, sd=sqrt(Sigma2))
    for(i in 2:n){
      u[i] <- Rho*u[i-1] + rnorm(1, mean=0, sd=sqrt(Sigma2))
    } 
    y <- rpois(n, exp(beta0 + u) )
    Data <- list( Options_vec = 4, 
                 y = y )
    Obj <- TMB::MakeADFun( data=Data, parameters=Pars, random="u", DLL="autoregressive" )
    # Optimize
    Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
  },
  replications = 100,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)
bnmk_n10[order(bnmk_n10$elapsed),]

#save(bnmk_n100, file = "2022_TMB_Session_II/R/bnmk_n100.Rdata")

load("2022_TMB_Session_II/R/bnmk_n100.Rdata")
bnmk_n100[order(bnmk_n100$elapsed),]


#Simulation Analysis ==============================================

n.sim <- 1000
n <- 100
Pars <- list( beta0 = 0, 
              ln_sigma2 = 0, 
              logit_rho = 0, 
              u = rnorm(n) )
nm <- paste0("SimStudy_n", ii)
beta0.est <- rho.est <- sigma2.est <- matrix(0,n.sim,5)

for(s in 1:n.sim){
  u <- rep(NA, n)
  u[1] <- rnorm(1, mean=0, sd=sqrt(Sigma2))
  for(i in 2:n){
    u[i] <- Rho*u[i-1] + rnorm(1, mean=0, sd=sqrt(Sigma2))
  } 
  y <- rpois(n, exp(beta0 + u) )
  Data <- list( Options_vec = 0, 
                y = y )
  for(i in 1:5){
    Data$Options_vec <- i-1
    Obj <- TMB::MakeADFun( data=Data, parameters=Pars, random="u", DLL="autoregressive" )
    Opt <- nlminb(Obj$par, Obj$fn, Obj$gr)
    Rep <- Obj$report()
    beta0.est[s,i] <- Opt$par['beta0']
    rho.est[s,i] <- Rep$rho
    sigma2.est[s,i] <- Rep$sigma2
  }
}
save(beta0.est, rho.est, sigma2.est, file = paste0("SimStudy_n1.RData"))

load("SimStudy_n1.RData")

beta0.mean <- apply(beta0.est,2,mean)
names(beta0.mean) <- cov.nms
beta0
beta0.mean

rho.mean <- apply(rho.est,2,mean)
names(rho.mean) <- cov.nms
Rho
rho.mean

sig2.mean <- apply(sigma2.est, 2, mean)
names(sig2.mean) <- cov.nms
Sigma2
sig2.mean
