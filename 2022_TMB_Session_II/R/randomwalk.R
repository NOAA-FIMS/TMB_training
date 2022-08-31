#randomwalk example
#dependencies: DHARMa, gap
#install.packages("DHARMa")
#install.packages("gap")
library(TMB)
N <- 100
mu <- 2
sd.vec <- c(1,1)

# simulate random measurements and data ====================
set.seed(123)
u <- c(0,cumsum(rnorm(N-1,mean=mu,sd=sd.vec[2])))
y <- u + rnorm(N,sd=sd.vec[1])


# setup TMB model ==========================================
compile('2022_TMB_Session_II/src/randomwalk.cpp')
dyn.load(dynlib('2022_TMB_Session_II/src/randomwalk'))

Dat <- list(y = y, sim_re = 0)
Pars <- list( u = rep(1, N), mu = 0, 
              ln_sig = 0, ln_tau = 0
            )
# Correctly specified model
obj0 <- TMB::MakeADFun( data=Dat, parameters=Pars, 
                       random="u", DLL="randomwalk" )
opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)

# Mis-specified model: fix mu=0
obj1 <- TMB::MakeADFun( data=Dat, parameters=Pars,
                        map = list(mu = factor(NA)),
                        random="u", DLL="randomwalk" )
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

#OSA methods
osa.fg.h0 <- oneStepPredict(obj0, observation.name = "y",
                         method = "fullGaussian")
osa.osg <- oneStepPredict(obj0, observation.name = "y",
                          data.term.indicator = "keep",
                          method = "oneStepGaussian")
osa.cdf <- oneStepPredict(obj0, observation.name = "y",
                          data.term.indicator = "keep",
                          method = "cdf")
osa.gen <- oneStepPredict(obj0, observation.name = "y",
                          data.term.indicator = "keep",
                          method = "oneStepGeneric")
osa.fg.h1 <- oneStepPredict(obj1, observation.name = "y",
                            method = "fullGaussian")
qqnorm(osa.fg.h0$residual);abline(0,1)
qqnorm(osa.fg.h1$residual);abline(0,1)

#DHARMa methods
#Only simulate from data model
set.seed(123)
obj0$env$data$sim_re <- 0
sim.y <- replicate(100, {obj0$simulate()$y})
ecdf.cond.h0 <- DHARMa::createDHARMa(sim.y, y)

obj1$env$data$sim_re <- 0
sim.y <- replicate(100, {obj1$simulate()$y})
ecdf.cond.h1 <- DHARMa::createDHARMa(sim.y, y)

#Simulate from data and RE model
obj0$env$data$sim_re <- 1
sim.y <- replicate(100, {obj0$simulate()$y})
ecdf.uncond.h0 <- DHARMa::createDHARMa(sim.y, y, rotation = "estimated")

obj1$env$data$sim_re <- 1
sim.y <- replicate(100, {obj1$simulate()$y})
ecdf.uncond.h1 <- DHARMa::createDHARMa(sim.y, y, rotation = "estimated")

#Correctly specified model
gap::qqunif(ecdf.cond.h0$scaledResiduals, logscale = FALSE);abline(0,1)
gap::qqunif(ecdf.uncond.h0$scaledResiduals, logscale = FALSE);abline(0,1)
#Mis-specified model
gap::qqunif(ecdf.cond.h1$scaledResiduals, logscale = FALSE);abline(0,1)
gap::qqunif(ecdf.uncond.h1$scaledResiduals, logscale = FALSE);abline(0,1)
