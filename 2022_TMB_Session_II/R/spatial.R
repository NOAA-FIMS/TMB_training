#package dependencies: mvtnorm
library(ggplot2)

# Matern covariance function 
cMatern <- function(H, Nu, R) {
  ifelse(H > 0, besselK(H/R, Nu) * (H/R)^Nu, 1) / gamma(Nu) * 2^(1-Nu)
}


## Simulate Data========================================

Dim = c("n_x"=20, "n_y"=20)
loc_xy = expand.grid("x"=1:Dim['n_x'], "y"=1:Dim['n_y'])
# distance matrix
Dmat <- as.matrix(dist(loc_xy), diag = TRUE, upper = TRUE)
sig2 <- 2
beta0 <- 1
Range <- 2 #related to distance at which spatial correlation ~10% 
Phi <- 1
Power <- 1.4

## Simulate spatial process
Sigma <- sig2 * cMatern(Dmat,1,Range)
plot(as.vector(Dmat), as.vector(Sigma))

set.seed(123)
omega <- t(mvtnorm::rmvnorm(1, rep(0,nrow(Sigma)), 
                        sigma = Sigma, method = 'chol'))
df <- data.frame(x = loc_xy$x, y = loc_xy$y, omega = omega)

ggplot(df, aes(x=x, y=y, fill = omega)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme_classic()
dev.off()
## Simulate Tweedie weights
mu  <- as.vector(exp(beta0 + omega))
y <- tweedie::rtweedie(nrow(loc_xy), mu = mu, 
                       phi = Phi, power = Power)
df$weight <- y
hist(y)

ggplot(df, aes(x=x, y=y, fill = weight)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme_classic()

total_weight <- sum(mu)

## Fit Model ==================================
library(TMB)
compile('2022_TMB_Session_II/src/spatial.cpp', framework = "TMBad")
dyn.load(dynlib('2022_TMB_Session_II/src/spatial'))

Data <- list(
  y = y, 
  D = Dmat
)

Pars <- list(
  beta0 = 0, 
  ln_phi = 0,
  logit_power = 0,
  ln_range = 0,
  ln_sigma2 = 0,
  omega = rep(0, length(y))
)
  
obj <- TMB::MakeADFun( 
  data=Data, 
  parameters=Pars, 
  random="omega", DLL="spatial" )
# Optimize
opt <- nlminb(obj$par, obj$fn, obj$gr)
report <- obj$report()

## explicitly calculate the laplace approximation
Hess <- obj$env$spHess(random = TRUE)
joint.nll <- report$nll
#log laplace approximation
-0.5 * log(2*pi)*length(y) + 0.5*log(det(as.matrix(Hess))) + joint.nll
-0.5 * log(2*pi)*length(y) + 0.5*sum(log(diag(as.matrix(Hess)))) + joint.nll
opt$objective


## check laplace approximation ===================
chk <- checkConsistency(obj)
summary(chk)


## Bias Correction ==============================
# Calculate SEs
sdr <- sdreport(obj)
summary(sdr,"fixed")
cbind(summary(sdr, "report"),  true = total_weight)

# Calculate with bias correction
sdr <- sdreport(obj, bias.correct = TRUE,
                bias.correct.control = list(sd = TRUE))
cbind(summary(sdr, "report"),  true = total_weight)

#Bias without bias correction
summary(sdr, "report")[1] - total_weight
#Bias eith bias correction
summary(sdr, "report")[3] - total_weight


conf.lower <- c(summary(sdr, "report")[c(1,3)]) - 
  1.96 * c(summary(sdr, "report")[c(2,4)])
conf.upper <- c(summary(sdr, "report")[c(1,3)]) + 
  1.96 * c(summary(sdr, "report")[c(2,4)])
df <- data.frame(est = c(summary(sdr, "report")[c(1,3)]),
                 lower = conf.lower,
                 upper = conf.upper,
                 type = c('Estimated', 'Bias Corrected'))
ggplot(df, aes(x = est, y = type )) + 
  geom_point() + 
  geom_errorbar(aes(xmin = lower, xmax = upper, y = type)) + 
  geom_vline(aes(xintercept = total_weight), col = 'red', show.legend = TRUE) +
 scale_color_manual(name = "", values = c(true = "red")) +
  theme_classic() + xlab('total weight')
