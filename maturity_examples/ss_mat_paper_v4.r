library(TMB)

#don't uncomment this, just paste it into R if you have used the debugger flags above.
#gdbsource("ss_mat_paper_v4.r")

#this allows dbs debugger to be used
dyn.unload("ss_mat_paper_v4.so")
compile("ss_mat_paper_v4.cpp","-O0 -g")
dyn.load("ss_mat_paper_v4.so")

mat.dat <- read.csv('../../cod_gb_19702014s_0yr_mvavg.csv', header = TRUE)
mat.dat = mat.dat[which(mat.dat$yc %in% 1963:2014),]
mat.dat = mat.dat[which(mat.dat$AGE > 0),]
fall.temp.dat = read.csv("../../1963-2014-fall-temp-for-state-space.csv", header = TRUE)
temp = aggregate(as.integer(!(mat.dat$MATURITY == 'I')), by = list(yc = mat.dat$yc, age = mat.dat$AGE), FUN = sum)
temp = cbind.data.frame(temp, aggregate(as.integer(!(mat.dat$MATURITY == 'I')), by = list(yc = mat.dat$yc, age = mat.dat$AGE), FUN = length)[[3]])
names(temp)[3:4] = c("Y","N")
x = list(Y = temp$Y)
x$N = temp$N
x$age_obs = temp$age
x$year_obs = temp$age + temp$yc - 1962
x$cy_obs = temp$yc - 1962
x$Ecov_maxages_k = rep(0,length(x$Y))
x$Ecov_maxages_a50 = rep(0,length(x$Y))
x$Ecov_obs = cbind(fall.temp.dat$anom_bt)
x$use_Ecov_obs = cbind(match(as.integer(!is.na(x$Ecov_obs)),c(0,1)) - 1)
x$Ecov_obs[which(is.na(x$Ecov_obs))] = -999
x$Ecov_obs_sigma = cbind(fall.temp.dat$sd1_bt)
x$Ecov_obs_sigma[which(is.na(x$Ecov_obs_sigma))] = -999
x$fit_k = 0
x$fit_a50 = 0
x$binomial = 1

y = list(beta_k = 0)
y$beta_a50 = 0
y$beta_Ecov_k = rep(0,max(1,x$Ecov_maxages_k+1))
y$beta_Ecov_a50 = rep(0,max(1,x$Ecov_maxages_a50+1))
y$Ecov_mu = 0
y$Ecov_AR_pars = rep(0, 2)
y$Ecov_re = cbind(c(x$Ecov_obs))
y$Ecov_re[which(x$use_Ecov_obs == 0)] = 0
y$beta_phi = 0
y$k_AR_pars = rep(0,2)
y$k_re = rep(0, length(y$Ecov_re))
y$a50_AR_pars = rep(0,2)
y$a50_re = rep(0, length(y$Ecov_re))

gbcod.mat.rev = list(dat = x, par = y)

gbcod.mat.rev$map = list(
    beta_Ecov_k = factor(rep(NA, length(gbcod.mat.rev$par$beta_Ecov_k))),
    beta_Ecov_a50 = factor(rep(NA, length(gbcod.mat.rev$par$beta_Ecov_a50))),
    beta_phi = factor(rep(NA, length(gbcod.mat.rev$par$beta_phi))),
    k_AR_pars = factor(rep(NA, length(gbcod.mat.rev$par$k_AR_pars))),
    k_re = factor(rep(NA, length(gbcod.mat.rev$par$k_re))),
    a50_AR_pars = factor(rep(NA, length(gbcod.mat.rev$par$a50_AR_pars))),
    a50_re = factor(rep(NA, length(gbcod.mat.rev$par$a50_re)))
    )

temp = gbcod.mat.rev
#no temperature effects, constant a50, k
x = MakeADFun(temp$dat,temp$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
mod0 = x
mod0$rep = mod0$report()
mod0$sdrep = sdreport(x)

#fall temperature effects at age 0 on k
temp = gbcod.mat.rev
temp$map = temp$map[which(names(temp$map) != "beta_Ecov_k")]
mod1 = MakeADFun(temp$dat,temp$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = temp$map)
mod1$opt = nlminb(mod1$par,mod1$fn,mod1$gr)
mod1$rep = mod1$report()
mod1$sdrep = sdreport(mod1)

#fall temperature effects at age 0 on a50
temp = gbcod.mat.rev
temp$map = temp$map[which(names(temp$map) != "beta_Ecov_a50")]
mod2 = MakeADFun(temp$dat,temp$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = temp$map)
mod2$opt = nlminb(mod2$par,mod2$fn,mod2$gr)
mod2$rep = mod2$report()
mod2$sdrep = sdreport(mod2)

#fall temperature effects at age 0 on k and a50
temp = gbcod.mat.rev
temp$map = temp$map[which(!(names(temp$map) %in% c("beta_Ecov_k","beta_Ecov_a50")))]
mod3 = MakeADFun(temp$dat,temp$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = temp$map)
mod3$opt = nlminb(mod3$par,mod3$fn,mod3$gr)
mod3$rep = mod3$report()
mod3$sdrep = sdreport(mod3)

#fall temperature effects at age 0 on a50 and 0,1 on k
gbcod.mat.rev.4 = gbcod.mat.rev
gbcod.mat.rev.4$dat$Ecov_maxages_k = rep(1,length(gbcod.mat.rev.4$dat$Y))
gbcod.mat.rev.4$dat$Ecov_maxages_k[gbcod.mat.rev.4$dat$age_obs==1] = 0
gbcod.mat.rev.4$par$beta_Ecov_k = rep(0,max(1,gbcod.mat.rev.4$dat$Ecov_maxages_k+1))
gbcod.mat.rev.4$map = gbcod.mat.rev.4$map[which(!(names(gbcod.mat.rev.4$map) %in% c("beta_Ecov_k","beta_Ecov_a50")))]
mod4 = MakeADFun(gbcod.mat.rev.4$dat,gbcod.mat.rev.4$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = gbcod.mat.rev.4$map)
mod4$opt = nlminb(mod4$par,mod4$fn,mod4$gr)
mod4$rep = mod4$report()
mod4$sdrep = sdreport(mod4)

#fall temperature effects at age 0,1 on a50 and 0 on k
#fall temperature effects at age 0 and 1 on k
gbcod.mat.rev.5 = gbcod.mat.rev.4
gbcod.mat.rev.5$dat$Ecov_maxages_k = rep(0,length(gbcod.mat.rev.5$dat$Y))
gbcod.mat.rev.5$par$beta_Ecov_k = rep(0,max(1,gbcod.mat.rev.5$dat$Ecov_maxages_k+1))
gbcod.mat.rev.5$dat$Ecov_maxages_a50 = rep(1,length(gbcod.mat.rev.5$dat$Y))
gbcod.mat.rev.5$dat$Ecov_maxages_a50[gbcod.mat.rev.5$dat$age_obs==1] = 0
gbcod.mat.rev.5$par$beta_Ecov_a50 = rep(0,max(1,gbcod.mat.rev.5$dat$Ecov_maxages_a50+1))
#gbcod.mat.rev.5$map = gbcod.mat.rev.5$map[which(!(names(gbcod.mat.rev.5$map) %in% c("beta_Ecov_k","beta_Ecov_a50")))]
mod5 = MakeADFun(gbcod.mat.rev.5$dat,gbcod.mat.rev.5$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = gbcod.mat.rev.5$map)
mod5$opt = nlminb(mod5$par,mod5$fn,mod5$gr)
mod5$rep = mod5$report()
mod5$sdrep = sdreport(mod5)

#fall temperature effects at ages 0,1 on a50 and k
gbcod.mat.rev.6 = gbcod.mat.rev.5
gbcod.mat.rev.6$dat$Ecov_maxages_k = rep(1,length(gbcod.mat.rev.6$dat$Y))
gbcod.mat.rev.6$dat$Ecov_maxages_k[gbcod.mat.rev.6$dat$age_obs==1] = 0
gbcod.mat.rev.6$par$beta_Ecov_k = rep(0,max(1,gbcod.mat.rev.6$dat$Ecov_maxages_k+1))
#gbcod.mat.rev.6$map = gbcod.mat.rev.6$map[which(!(names(gbcod.mat.rev.6$map) %in% c("beta_Ecov_k","beta_Ecov_a50")))]
mod6 = MakeADFun(gbcod.mat.rev.6$dat,gbcod.mat.rev.6$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = gbcod.mat.rev.6$map)
mod6$opt = nlminb(mod6$par,mod6$fn,mod6$gr)
mod6$rep = mod6$report()
mod6$sdrep = sdreport(mod6)

#fall temperature effects at ages 0,1 on a50 and 0,1,2 on k
gbcod.mat.rev.7 = gbcod.mat.rev.6
gbcod.mat.rev.7$dat$Ecov_maxages_k = rep(2,length(gbcod.mat.rev.7$dat$Y))
gbcod.mat.rev.7$dat$Ecov_maxages_k[gbcod.mat.rev.7$dat$age_obs==1] = 0
gbcod.mat.rev.7$dat$Ecov_maxages_k[gbcod.mat.rev.7$dat$age_obs==2] = 1
gbcod.mat.rev.7$par$beta_Ecov_k = rep(0,max(1,gbcod.mat.rev.7$dat$Ecov_maxages_k+1))
#gbcod.mat.rev.7$map = gbcod.mat.rev.7$map[which(!(names(gbcod.mat.rev.7$map) %in% c("beta_Ecov_k","beta_Ecov_a50")))]
mod7 = MakeADFun(gbcod.mat.rev.7$dat,gbcod.mat.rev.7$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = gbcod.mat.rev.7$map)
mod7$opt = nlminb(mod7$par,mod7$fn,mod7$gr)
mod7$rep = mod7$report()
mod7$sdrep = sdreport(mod7)

#fall temperature effects at ages 0,1,2 on a50 and 0,1 on k
gbcod.mat.rev.8 = gbcod.mat.rev.7
gbcod.mat.rev.8$dat$Ecov_maxages_k = rep(1,length(gbcod.mat.rev.8$dat$Y))
gbcod.mat.rev.8$dat$Ecov_maxages_k[gbcod.mat.rev.8$dat$age_obs==1] = 0
gbcod.mat.rev.8$par$beta_Ecov_k = rep(0,max(1,gbcod.mat.rev.8$dat$Ecov_maxages_k+1))
gbcod.mat.rev.8$dat$Ecov_maxages_a50 = rep(2,length(gbcod.mat.rev.8$dat$Y))
gbcod.mat.rev.8$dat$Ecov_maxages_a50[gbcod.mat.rev.8$dat$age_obs==1] = 0
gbcod.mat.rev.8$dat$Ecov_maxages_a50[gbcod.mat.rev.8$dat$age_obs==2] = 1
gbcod.mat.rev.8$par$beta_Ecov_a50 = rep(0,max(1,gbcod.mat.rev.8$dat$Ecov_maxages_a50+1))
mod8 = MakeADFun(gbcod.mat.rev.8$dat,gbcod.mat.rev.8$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = gbcod.mat.rev.8$map)
mod8$opt = nlminb(mod8$par,mod8$fn,mod8$gr)
mod8$rep = mod8$report()
mod8$sdrep = sdreport(mod8)

#fall temperature effects at ages 0,1,2 on a50 and 0,1,2 on k
gbcod.mat.rev.9 = gbcod.mat.rev.8
gbcod.mat.rev.9$dat$Ecov_maxages_k = rep(2,length(gbcod.mat.rev.9$dat$Y))
gbcod.mat.rev.9$dat$Ecov_maxages_k[gbcod.mat.rev.9$dat$age_obs==1] = 0
gbcod.mat.rev.9$dat$Ecov_maxages_k[gbcod.mat.rev.9$dat$age_obs==2] = 1
gbcod.mat.rev.9$par$beta_Ecov_k = rep(0,max(1,gbcod.mat.rev.9$dat$Ecov_maxages_k+1))
mod9 = MakeADFun(gbcod.mat.rev.9$dat,gbcod.mat.rev.9$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = gbcod.mat.rev.9$map)
mod9$opt = nlminb(mod9$par,mod9$fn,mod9$gr)
mod9$rep = mod9$report()
mod9$sdrep = sdreport(mod9)

#betabinomial models
dyn.unload("ss_mat_paper_v4.so")
compile("ss_mat_paper_v4.cpp","-O0 -g")
dyn.load("ss_mat_paper_v4.so")

temp = gbcod.mat.rev
temp$dat$binomial = 0
temp$map = temp$map[which(names(temp$map) != "beta_phi")]
x = MakeADFun(temp$dat,temp$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
mod0.bb = x
mod0.bb$rep = mod0.bb$report()
mod0.bb$sdrep = sdreport(x)

x = mod0.bb$simulate(complete=TRUE)
temp = list(mod0.bb$rep$mat, exp(mod0.bb$opt$par['beta_phi']))
temp = rbeta(length(temp[[1]]), temp[[1]]*temp[[2]], (1-temp[[1]])*temp[[2]])
temp = rbinom(length(temp), size = mod0.bb$env$data$N, prob = temp)
x = mod0.bb$env$data
x$Y = as.numeric(temp) 
x = MakeADFun(x,mod0.bb$env$parList(),random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = mod0.bb$env$map, silent = TRUE)
nlminb(x$par,x$fn,x$gr)$par
mod0.bb$opt$par

temp = gbcod.mat.rev
temp$dat$binomial = 0
temp$map = temp$map[which(names(temp$map) != "beta_Ecov_k")]
temp$map = temp$map[which(names(temp$map) != "beta_phi")]
x = MakeADFun(temp$dat,temp$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
mod1.bb = x
mod1.bb$rep = mod1.bb$report()
mod1.bb$sdrep = sdreport(x)

temp = gbcod.mat.rev
temp$dat$binomial = 0
temp$map = temp$map[which(names(temp$map) != "beta_Ecov_a50")]
temp$map = temp$map[which(names(temp$map) != "beta_phi")]
x = MakeADFun(temp$dat,temp$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
mod2.bb = x
mod2.bb$rep = mod2.bb$report()
mod2.bb$sdrep = sdreport(x)

temp = gbcod.mat.rev
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("beta_Ecov_k","beta_Ecov_a50","beta_phi")))]
x = MakeADFun(temp$dat,temp$par,random=c("Ecov_re"),DLL="ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
mod3.bb = x
mod3.bb$rep = mod3.bb$report()
mod3.bb$sdrep = sdreport(x)

#binomial models with AR1 in k or a50

#AR1 in k at age 0
temp = gbcod.mat.rev
#temp$dat$Ecov_maxages_k = rep(temp$dat$age_obs-1,length(temp$dat$Y))
#temp$par$beta_Ecov_k = rep(0, max(temp$dat$age_obs))
temp$dat$fit_k = 1
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re")))]
#temp$map$beta_Ecov_k = factor(rep(NA,length(temp$par$beta_Ecov_k)))
mod0.k.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod0.k.AR$opt = nlminb(mod0.k.AR$par,mod0.k.AR$fn,mod0.k.AR$gr)
#801.72

#AR1 in a50 at age 0
temp = gbcod.mat.rev
#temp$dat$Ecov_maxages_a50 = rep(temp$dat$age_obs-1,length(temp$dat$Y))
#temp$par$beta_Ecov_a50 = rep(0, max(temp$dat$age_obs))
temp$dat$fit_a50 = 1
temp$map = temp$map[which(!(names(temp$map) %in% c("a50_AR_pars","a50_re")))]
#temp$map$beta_Ecov_a50 = factor(rep(NA,length(temp$par$beta_Ecov_a50)))
mod0.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod0.a50.AR$opt = nlminb(mod0.a50.AR$par,mod0.a50.AR$fn,mod0.a50.AR$gr)
#665.97

#AR1 in k at age 0:1 and a50 age 0
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_k = rep(1,length(temp$dat$Y))
temp$dat$Ecov_maxages_k[temp$dat$age_obs==1] = 0
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par$beta_Ecov_k = rep(0,max(1,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
mod1.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod1.k.a50.AR$opt = nlminb(mod1.k.a50.AR$par,mod1.k.a50.AR$fn,mod1.k.a50.AR$gr)
#669.85

#AR1 in a50 at age 0:1 and k age 0
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(1,length(temp$dat$Y))
temp$dat$Ecov_maxages_a50[temp$dat$age_obs==1] = 0
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par$beta_Ecov_a50 = rep(0,max(1,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
mod2.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod2.k.a50.AR$opt = nlminb(mod2.k.a50.AR$par,mod2.k.a50.AR$fn,mod2.k.a50.AR$gr)
#665.02

#AR1 in k at age 0:1
#AR1 in a50 at age 0:1
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(1,length(temp$dat$Y))
temp$dat$Ecov_maxages_a50[temp$dat$age_obs==1] = 0
temp$dat$Ecov_maxages_k = rep(1,length(temp$dat$Y))
temp$dat$Ecov_maxages_k[temp$dat$age_obs==1] = 0
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par$beta_Ecov_a50 = rep(0,max(1,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(1,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
mod3.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod3.k.a50.AR$opt = nlminb(mod3.k.a50.AR$par,mod3.k.a50.AR$fn,mod3.k.a50.AR$gr)
#662.09

#AR1 in k at age 0:2
#AR1 in a50 at age 0:2
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(2,length(temp$dat$Y))
temp$dat$Ecov_maxages_a50[temp$dat$age_obs==2] = 1
temp$dat$Ecov_maxages_a50[temp$dat$age_obs==1] = 0
temp$dat$Ecov_maxages_k = rep(2,length(temp$dat$Y))
temp$dat$Ecov_maxages_k[temp$dat$age_obs==2] = 1
temp$dat$Ecov_maxages_k[temp$dat$age_obs==1] = 0
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par$beta_Ecov_a50 = rep(0,max(2,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(2,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
mod4.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod4.k.a50.AR$opt = nlminb(mod4.k.a50.AR$par,mod4.k.a50.AR$fn,mod4.k.a50.AR$gr)
mod4.k.a50.AR$rep = mod4.k.a50.AR$report()
x = mod4.k.a50.AR
mod4.k.a50.AR$sdrep = sdreport(x)
#626.39

temp = list(dat = mod4.k.a50.AR$env$data, par = mod4.k.a50.AR$env$parList(), map = mod4.k.a50.AR$env$map)
temp$map$beta_Ecov_k = factor(c(1,NA,NA))
#temp$map$beta_Ecov_a50 = factor(c(1,NA,NA))
mod4.k.a50.AR.Ecov.1 = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod4.k.a50.AR.Ecov.1$opt = nlminb(mod4.k.a50.AR.Ecov.1$par,mod4.k.a50.AR.Ecov.1$fn,mod4.k.a50.AR.Ecov.1$gr)
#625.93

temp = list(dat = mod4.k.a50.AR$env$data, par = mod4.k.a50.AR$env$parList(), map = mod4.k.a50.AR$env$map)
#temp$map$beta_Ecov_k = factor(c(1,NA,NA))
temp$map$beta_Ecov_a50 = factor(c(1,NA,NA))
mod4.k.a50.AR.Ecov.2 = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod4.k.a50.AR.Ecov.2$opt = nlminb(mod4.k.a50.AR.Ecov.2$par,mod4.k.a50.AR.Ecov.2$fn,mod4.k.a50.AR.Ecov.2$gr)
#626.03

temp = list(dat = mod4.k.a50.AR$env$data, par = mod4.k.a50.AR$env$parList(), map = mod4.k.a50.AR$env$map)
temp$map$beta_Ecov_k = factor(c(1,NA,NA))
temp$map$beta_Ecov_a50 = factor(c(1,NA,NA))
mod4.k.a50.AR.Ecov.3 = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod4.k.a50.AR.Ecov.3$opt = nlminb(mod4.k.a50.AR.Ecov.3$par,mod4.k.a50.AR.Ecov.3$fn,mod4.k.a50.AR.Ecov.3$gr)
#625.75


#AR1 in k at age 0:3 
#AR1 in a50 at age 0:3
maxage = 3
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
x = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
#610.08

#AR1 in k at age 0:4 
#AR1 in a50 at age 0:4
maxage = 4
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par = mod4.k.a50.AR$env$parList()
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
x = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
#602.73

#AR1 in k at age 0:5 
#AR1 in a50 at age 0:5
maxage = 5
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par = mod4.k.a50.AR$env$parList()
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
x = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
#601.68

#AR1 in k at age 0:6 
#AR1 in a50 at age 0:6
maxage = 6
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par = mod4.k.a50.AR$env$parList()
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
x = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
#601.02

#AR1 in k at age 0:7 
#AR1 in a50 at age 0:7
maxage = 7
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par = mod4.k.a50.AR$env$parList()
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
mod5.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod5.k.a50.AR$opt = nlminb(mod5.k.a50.AR$par,mod5.k.a50.AR$fn,mod5.k.a50.AR$gr)
#600.67

#AR1 in k at age 0:8 
#AR1 in a50 at age 0:7
temp = gbcod.mat.rev
temp$par = mod4.k.a50.AR$env$parList()
maxage = 7
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
maxage = 8
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
mod6.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod6.k.a50.AR$opt = nlminb(mod6.k.a50.AR$par,mod6.k.a50.AR$fn,mod6.k.a50.AR$gr)
#600.69

#AR1 in k at age 0:7 
#AR1 in a50 at age 0:8
temp = gbcod.mat.rev
temp$par = mod4.k.a50.AR$env$parList()
maxage = 8
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
maxage = 7
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
mod7.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod7.k.a50.AR$opt = nlminb(mod7.k.a50.AR$par,mod7.k.a50.AR$fn,mod7.k.a50.AR$gr)
mod7.k.a50.AR$rep = mod7.k.a50.AR$report()
x = mod7.k.a50.AR
mod7.k.a50.AR$sdrep = sdreport(x)
#600.66

#AR1 in k at age 0:8 
#AR1 in a50 at age 0:8
maxage = 8
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$par = mod4.k.a50.AR$env$parList()
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
mod8.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod8.k.a50.AR$opt = nlminb(mod8.k.a50.AR$par,mod8.k.a50.AR$fn,mod8.k.a50.AR$gr)
#600.68

#betabinomial models with AR1 in k or a50

#AR1 in k at age 0
temp = gbcod.mat.rev
#temp$dat$Ecov_maxages_k = rep(temp$dat$age_obs-1,length(temp$dat$Y))
#temp$par = mod0.k.AR$env$parList()
#temp$par$beta_Ecov_k = rep(0, max(temp$dat$age_obs))
temp$dat$fit_k = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re")))]
temp$map = temp$map[which(!(names(temp$map) %in% c("beta_phi")))]
#temp$map$beta_Ecov_k = factor(rep(NA,length(temp$par$beta_Ecov_k)))
mod0.bb.k.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod0.bb.k.AR$opt = nlminb(mod0.bb.k.AR$par,mod0.bb.k.AR$fn,mod0.bb.k.AR$gr)
temp$map = temp$map[which(!(names(temp$map) %in% c("beta_phi")))]
temp$par = mod0.bb.k.AR$env$parList()
mod0.bb.k.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod0.bb.k.AR$opt = nlminb(mod0.bb.k.AR$par,mod0.bb.k.AR$fn,mod0.bb.k.AR$gr)

#AR1 in a50 at age 0
temp = gbcod.mat.rev
#temp$dat$Ecov_maxages_a50 = rep(temp$dat$age_obs-1,length(temp$dat$Y))
#temp$par = mod0.a50.AR$env$parList()
#temp$par$beta_Ecov_a50 = rep(0, max(temp$dat$age_obs))
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("a50_AR_pars","a50_re")))]
#temp$map$beta_Ecov_a50 = factor(rep(NA,length(temp$par$beta_Ecov_a50)))
mod0.bb.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod0.bb.a50.AR$opt = nlminb(mod0.bb.a50.AR$par,mod0.bb.a50.AR$fn,mod0.bb.a50.AR$gr)
temp$map = temp$map[which(!(names(temp$map) %in% c("beta_phi")))]
temp$par = mod0.bb.a50.AR$env$parList()
mod0.bb.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod0.bb.a50.AR$opt = nlminb(mod0.bb.a50.AR$par,mod0.bb.a50.AR$fn,mod0.bb.a50.AR$gr)

#AR1 in k and a50 at age 0
temp = gbcod.mat.rev
#temp$dat$Ecov_maxages_k = rep(temp$dat$age_obs-1,length(temp$dat$Y))
#temp$dat$Ecov_maxages_a50 = rep(temp$dat$age_obs-1,length(temp$dat$Y))
#temp$par = mod0.bb.k.AR$env$parList()
#temp$par$beta_Ecov_k = rep(0, max(temp$dat$age_obs))
#temp$par$beta_Ecov_a50 = rep(0, max(temp$dat$age_obs))
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
temp$map = temp$map[which(!(names(temp$map) %in% c("beta_phi")))]
#temp$map$beta_Ecov_k = factor(rep(NA,length(temp$par$beta_Ecov_k)))
#temp$map$beta_Ecov_a50 = factor(rep(NA,length(temp$par$beta_Ecov_a50)))
mod0.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod0.bb.k.a50.AR$opt = nlminb(mod0.bb.k.a50.AR$par,mod0.bb.k.a50.AR$fn,mod0.bb.k.a50.AR$gr)
#temp$map = temp$map[which(!(names(temp$map) %in% c("beta_phi")))]
#temp$par = mod0.bb.k.a50.AR$env$parList()
#mod0.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
#mod0.bb.k.a50.AR$opt = nlminb(mod0.bb.k.a50.AR$par,mod0.bb.k.a50.AR$fn,mod0.bb.k.a50.AR$gr)
mod0.bb.k.a50.AR$rep = mod0.bb.k.a50.AR$report()
x = mod0.bb.k.a50.AR
mod0.bb.k.a50.AR$sdrep = sdreport(x)
#571.70

#AR1 in k at age 0:1 and a50 age 0
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_k = rep(1,length(temp$dat$Y))
temp$dat$Ecov_maxages_k[temp$dat$age_obs==1] = 0
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$par$beta_Ecov_k = rep(0,max(1,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re")))]
temp$map = temp$map[which(!(names(temp$map) %in% c("beta_phi")))]
mod1.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod1.bb.k.a50.AR$opt = nlminb(mod1.bb.k.a50.AR$par,mod1.bb.k.a50.AR$fn,mod1.bb.k.a50.AR$gr)
mod1.bb.k.a50.AR$rep = mod1.bb.k.a50.AR$report()
x = mod1.bb.k.a50.AR
mod1.bb.k.a50.AR$sdrep = sdreport(x)
#556.56

#AR1 in a50 at age 0:1 and k age 0
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(1,length(temp$dat$Y))
temp$dat$Ecov_maxages_a50[temp$dat$age_obs==1] = 0
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$par$beta_Ecov_a50 = rep(0,max(1,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
mod2.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod2.bb.k.a50.AR$opt = nlminb(mod2.bb.k.a50.AR$par,mod2.bb.k.a50.AR$fn,mod2.bb.k.a50.AR$gr)
mod2.bb.k.a50.AR$rep = mod2.bb.k.a50.AR$report()
x = mod2.bb.k.a50.AR
mod2.bb.k.a50.AR$sdrep = sdreport(x)
#554.73

#AR1 in k at age 0:1
#AR1 in a50 at age 0:1
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(1,length(temp$dat$Y))
temp$dat$Ecov_maxages_a50[temp$dat$age_obs==1] = 0
temp$dat$Ecov_maxages_k = rep(1,length(temp$dat$Y))
temp$dat$Ecov_maxages_k[temp$dat$age_obs==1] = 0
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$par$beta_Ecov_a50 = rep(0,max(1,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(1,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
mod3.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod3.bb.k.a50.AR$opt = nlminb(mod3.bb.k.a50.AR$par,mod3.bb.k.a50.AR$fn,mod3.bb.k.a50.AR$gr)
mod3.bb.k.a50.AR$rep = mod3.bb.k.a50.AR$report()
x = mod3.bb.k.a50.AR
mod3.bb.k.a50.AR$sdrep = sdreport(x)
#554.69

#AR1 in k at age 0:1
#AR1 in a50 at age 0:2
maxage = 2
temp = gbcod.mat.rev
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
maxage = 1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
x = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
mod4.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod4.bb.k.a50.AR$opt = nlminb(mod4.bb.k.a50.AR$par,mod4.bb.k.a50.AR$fn,mod4.bb.k.a50.AR$gr)
mod4.bb.k.a50.AR$rep = mod4.bb.k.a50.AR$report()
x = mod4.bb.k.a50.AR
mod4.bb.k.a50.AR$sdrep = sdreport(x)
#545.03
mod4.bb.k.a50.AR$parList = mod4.bb.k.a50.AR$env$parList()
save(mod4.bb.k.a50.AR, file = "mod4.bb.k.a50.AR.RData")

temp = list(dat = mod4.bb.k.a50.AR$env$data, par = mod4.bb.k.a50.AR$env$parList(), map = mod4.bb.k.a50.AR$env$map)
temp$map$beta_Ecov_k = factor(c(1,NA))
#temp$map$beta_Ecov_a50 = factor(c(1,NA,NA))
mod4.bb.k.a50.AR.Ecov.1 = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod4.bb.k.a50.AR.Ecov.1$opt = nlminb(mod4.bb.k.a50.AR.Ecov.1$par,mod4.bb.k.a50.AR.Ecov.1$fn,mod4.bb.k.a50.AR.Ecov.1$gr)
#544.97

temp = list(dat = mod4.bb.k.a50.AR$env$data, par = mod4.bb.k.a50.AR$env$parList(), map = mod4.bb.k.a50.AR$env$map)
#temp$map$beta_Ecov_k = factor(c(1,NA))
temp$map$beta_Ecov_a50 = factor(c(1,NA,NA))
mod4.bb.k.a50.AR.Ecov.2 = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod4.bb.k.a50.AR.Ecov.2$opt = nlminb(mod4.bb.k.a50.AR.Ecov.2$par,mod4.bb.k.a50.AR.Ecov.2$fn,mod4.bb.k.a50.AR.Ecov.2$gr)
#544.97

temp = list(dat = mod4.bb.k.a50.AR$env$data, par = mod4.bb.k.a50.AR$env$parList(), map = mod4.bb.k.a50.AR$env$map)
temp$map$beta_Ecov_k = factor(c(1,NA))
temp$map$beta_Ecov_a50 = factor(c(1,NA,NA))
mod4.bb.k.a50.AR.Ecov.3 = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod4.bb.k.a50.AR.Ecov.3$opt = nlminb(mod4.bb.k.a50.AR.Ecov.3$par,mod4.bb.k.a50.AR.Ecov.3$fn,mod4.bb.k.a50.AR.Ecov.3$gr)
#544.93


#AR1 in k at age 0:1
#AR1 in a50 at age 0:2
maxage = 1
temp = gbcod.mat.rev
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
maxage = 2
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
mod5.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod5.bb.k.a50.AR$opt = nlminb(mod5.bb.k.a50.AR$par,mod5.bb.k.a50.AR$fn,mod5.bb.k.a50.AR$gr)
mod5.bb.k.a50.AR$rep = mod5.bb.k.a50.AR$report()
x = mod5.bb.k.a50.AR
mod5.bb.k.a50.AR$sdrep = sdreport(x)
#554.52

#AR1 in k at age 0:2
#AR1 in a50 at age 0:2
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(2,length(temp$dat$Y))
temp$dat$Ecov_maxages_a50[temp$dat$age_obs==2] = 1
temp$dat$Ecov_maxages_a50[temp$dat$age_obs==1] = 0
temp$dat$Ecov_maxages_k = rep(2,length(temp$dat$Y))
temp$dat$Ecov_maxages_k[temp$dat$age_obs==2] = 1
temp$dat$Ecov_maxages_k[temp$dat$age_obs==1] = 0
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$par$beta_Ecov_a50 = rep(0,max(2,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(2,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
mod6.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod6.bb.k.a50.AR$opt = nlminb(mod6.bb.k.a50.AR$par,mod6.bb.k.a50.AR$fn,mod6.bb.k.a50.AR$gr)
mod6.bb.k.a50.AR$rep = mod6.bb.k.a50.AR$report()
x = mod6.bb.k.a50.AR
mod6.bb.k.a50.AR$sdrep = sdreport(x)
#545.01

#AR1 in k at age 0:2
#AR1 in a50 at age 0:3
maxage = 3
temp = gbcod.mat.rev
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
maxage = 2
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
x = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
x$opt = nlminb(x$par,x$fn,x$gr)
x$opt$obj
#548.49

#AR1 in k at age 0:3 
#AR1 in a50 at age 0:2
maxage = 2
temp = gbcod.mat.rev
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
maxage = 3
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
mod7.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod7.bb.k.a50.AR$opt = nlminb(mod7.bb.k.a50.AR$par,mod7.bb.k.a50.AR$fn,mod7.bb.k.a50.AR$gr)
mod7.bb.k.a50.AR$rep = mod7.bb.k.a50.AR$report()
x = mod7.bb.k.a50.AR
mod7.bb.k.a50.AR$sdrep = sdreport(x)
#544.97
mod7.bb.k.a50.AR$parList = mod7.bb.k.a50.AR$env$parList()
save(mod7.bb.k.a50.AR, file = "mod7.bb.k.a50.AR.RData")

x = summary(mod7.bb.k.a50.AR$sdrep)
x[rownames(x) %in% "a50_re",1]
x = x[rownames(x) %in% c("Ecov_y","k_re","a50_re","beta_a50","beta_k","beta_Ecov_k", "beta_Ecov_a50"),1]
ind = sum(names(x) == "Ecov_y")
ind = c(1:(sum(names(x) == "Ecov_y")-2),ind)
x = x[-which(names(x) == "Ecov_y")[ind]]
x = x[-which(names(x) == "k_re")[ind]]
x = x[-which(names(x) == "a50_re")[ind]]
x
k = exp(x["beta_k"] + x["k_re"])
a50 = exp(x["beta_a50"] + x["a50_re"])
1/(1 + exp(-k*(1-a50)))

x = summary(mod7.bb.k.a50.AR$sdrep)
x = x[rownames(x) %in% "a50_re",1]

exp(temp["beta_Linf"])*(1-exp(-k1*(1-temp["t0"])))


#AR1 in k at age 0:3 
#AR1 in a50 at age 0:3
maxage = 3
temp = gbcod.mat.rev
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
mod8.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod8.bb.k.a50.AR$opt = nlminb(mod8.bb.k.a50.AR$par,mod8.bb.k.a50.AR$fn,mod8.bb.k.a50.AR$gr)
mod8.bb.k.a50.AR$opt$obj
#546.14
#not as good

#AR1 in k at age 0:4 
#AR1 in a50 at age 0:2
maxage = 2
temp = gbcod.mat.rev
temp$par = mod0.bb.k.a50.AR$env$parList()
temp$dat$Ecov_maxages_a50 = rep(maxage,length(temp$dat$Y))
if(maxage>0) for(i in 1:maxage) temp$dat$Ecov_maxages_a50[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_a50 = rep(0,max(maxage,temp$dat$Ecov_maxages_a50+1))
temp$map$beta_Ecov_a50 = factor(rep(NA, length(temp$par$beta_Ecov_a50)))
maxage = 4
temp$dat$Ecov_maxages_k = rep(maxage,length(temp$dat$Y))
if(maxage>0) for(i in 1:maxage) temp$dat$Ecov_maxages_k[temp$dat$age_obs==i] = i-1
temp$par$beta_Ecov_k = rep(0,max(maxage,temp$dat$Ecov_maxages_k+1))
temp$map$beta_Ecov_k = factor(rep(NA, length(temp$par$beta_Ecov_k)))
temp$dat$fit_k = 1
temp$dat$fit_a50 = 1
temp$dat$binomial = 0
temp$map = temp$map[which(!(names(temp$map) %in% c("k_AR_pars","k_re","a50_AR_pars","a50_re","beta_phi")))]
mod9.bb.k.a50.AR = MakeADFun(temp$dat, temp$par, random = c("Ecov_re","k_re","a50_re"), DLL = "ss_mat_paper_v4", map = temp$map)
mod9.bb.k.a50.AR$opt = nlminb(mod9.bb.k.a50.AR$par,mod9.bb.k.a50.AR$fn,mod9.bb.k.a50.AR$gr)
mod9.bb.k.a50.AR$opt$obj
#545.06
#not as good


temp = list(mod0$opt,mod1$opt,mod2$opt,mod3$opt,mod4$opt,mod5$opt,mod6$opt,mod7$opt,mod8$opt,mod9$opt)
temp = sapply(temp, function(x) 2 * x$obj + 2 * length(x$par))
temp - min(temp)
#temp = list(mod0.bb$opt,mod1.bb$opt,mod2.bb$opt,mod3.bb$opt,mod4.bb$opt,mod5.bb$opt,mod6.bb$opt,mod7.bb$opt,mod8.bb$opt,mod9.bb$opt)
temp = list(mod0.bb$opt,mod1.bb$opt,mod2.bb$opt,mod3.bb$opt)
temp = sapply(temp, function(x) 2 * x$obj + 2 * length(x$par))
temp - min(temp)
#temp = list(mod0.bb$opt,mod1.bb$opt,mod2.bb$opt,mod3.bb$opt,mod4.bb$opt,mod5.bb$opt,mod6.bb$opt,mod7.bb$opt,mod8.bb$opt,mod9.bb$opt,
temp = list(mod0$opt,mod1$opt,mod2$opt,mod3$opt,mod4$opt,mod5$opt,mod6$opt,mod7$opt,mod8$opt,mod9$opt,
  mod0.bb$opt,mod1.bb$opt,mod2.bb$opt,mod3.bb$opt)
aic = sapply(temp, function(x) 2 * x$obj + 2 * length(x$par))
delta.aic = aic - min(aic)
aic.wts = exp(-0.5*delta.aic)/sum(exp(-0.5*delta.aic))

temp = list(
  mod0$opt,mod7$opt,mod4.k.a50.AR$opt,mod4.k.a50.AR.Ecov.1$opt,mod4.k.a50.AR.Ecov.2$opt,mod4.k.a50.AR.Ecov.3$opt,
  mod0.bb$opt,mod1.bb$opt,mod2.bb$opt,mod3.bb$opt,mod4.bb.k.a50.AR$opt,
  mod4.bb.k.a50.AR.Ecov.1$opt,mod4.bb.k.a50.AR.Ecov.2$opt,mod4.bb.k.a50.AR.Ecov.3$opt)
sapply(temp, function(x) length(x$par))
aic = sapply(temp, function(x) 2 * x$obj + 2 * length(x$par))
aic - min(aic)


library(plotrix)
tcol <- col2rgb('black')
tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')


mat.dat.table = table(mat.dat$YEAR, mat.dat$AGE)
temp = apply(mat.dat.table[,-(1:8)], 1, sum)
mat.dat.table = cbind(mat.dat.table[,1:8], "$>$ 8" = temp)
mat.dat.table = cbind(mat.dat.table, Weight = sapply(unique(mat.dat$YEAR), function(x) sum(mat.dat$YEAR == x & mat.dat$INDWT>0)))
#mat.dat.table = cbind(mat.dat.table, Length = sapply(unique(mat.dat$YEAR), function(x) sum(mat.dat$YEAR == x & mat.dat$LENGTH>0)))
mat.dat.table = cbind(mat.dat.table, Total = sapply(unique(mat.dat$YEAR), function(x) sum(mat.dat$YEAR == x)))
mat.dat.table = rbind(mat.dat.table, Total = apply(mat.dat.table, 2, sum))

library(Hmisc)
x = latex(mat.dat.table, file = '~/work/cod/ss_maturity/tex/mat_data_table_v2.tex', 
  rowlabel = 'Year', cgroup = c("Age", "", ""), n.cgroup = c(9,1,1), table.env = FALSE, rowlabel.just = "c", cgroupTexCmd = NULL)#, collabel.just)

x = cbind(gbcod.mat$dat$Ecov_obs, gbcod.mat$dat$Ecov_obs_sigma)
x[which(x < -99)] = NA
colnames(x) = c("Anomaly", "Standard Error")
rownames(x) = 1963:2014
temp = latex(x, file = '~/work/cod/ss_maturity/tex/bottom_temperature_anomalies.tex', 
  rowlabel = 'Year', table.env = FALSE, rowlabel.just = "c")#, collabel.just)
temp = latex(x[1963:2000-1962,], file = '~/work/cod/ss_maturity/tex/bottom_temperature_anomalies_1.tex', 
  rowlabel = 'Year', table.env = FALSE, rowlabel.just = "c")#, collabel.just)
temp = latex(x[2001:2014-1962,], file = '~/work/cod/ss_maturity/tex/bottom_temperature_anomalies_2.tex', 
  rowlabel = 'Year', table.env = FALSE, rowlabel.just = "c")#, collabel.just)

par(mfrow = c(1,1), mar = c(1,4,1,1), oma = c(4,1,0,0))
temp = summary(mod0.sdrep)
temp = temp[which(rownames(temp) == "Ecov_y"),]
temp = cbind(temp[,1], temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2]))
plot(1963:2014, temp[,1], type = 'n', axes = FALSE, ylim = range(temp), xlab = "", ylab = "")
grid(col = gray(0.7), lwd = 1)
lines(1963:2014, temp[,1], lwd = 2)
polygon(c(1963:2014,2014:1963), c(temp[,2],rev(temp[,3])), col = tcol, border = "transparent", lty = 2)
x = gbcod.mat$dat
temp <- which(x$use_Ecov_obs == 1)
plotCI((1963:2014)[temp], x$Ecov_obs[temp], li = (x$Ecov_obs - qnorm(0.975)*x$Ecov_obs_sigma)[temp], ui = (x$Ecov_obs + qnorm(0.975)*x$Ecov_obs_sigma)[temp], 
  add = TRUE, lwd = 2)
axis(1, labels = FALSE, lwd = 2)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
mtext(side = 2, outer = FALSE, line = 3, "Bottom temperature anomaly", cex = 1.5)

###########################################################
#plots comparing maturity predictions and emperical estimates

y = summary(mod7.bb.k.a50.AR$sdrep)
y = y[rownames(y) == "logit_pmat",]
z = matrix(1/(1+exp(-y[,1])), ncol = 3)
plot((1963:2014)[-(1:3)],z[-(1:3),3], ylim = c(0,1), xlab = '', ylab = '', type = 'l')
lines((1963:2014)[-1],z[-1,1], col = 'red')
plot((1963:2014)[-(1:2)],z[-(1:2),2], ylim = c(0,1), xlab = '', ylab = '', type = 'n')

temp.fn = function(whichage = 2, model = mod7.bb.k.a50.AR)
{
  library(plotrix)
  y = summary(model$sdrep)
  #y = summary(sdreport(x))
  y = y[rownames(y) == "logit_pmat",]
  y = cbind(y[,1], y[,1] +qnorm(0.975)*cbind(-y[,2],y[,2]))
  z = cbind(matrix(1/(1+exp(-y[,1])), ncol = 3)[,whichage], matrix(1/(1+exp(-y[,2])), ncol = 3)[,whichage],matrix(1/(1+exp(-y[,3])), ncol = 3)[,whichage])
  years = 1963:2014
  years.dat = 1970:2014
  ind = which(years %in% years.dat)#[-(1:whichage)]
  plot(years.dat,z[ind,2], ylim = c(0,1), xlab = '', ylab = '', type = 'n', axes = FALSE, xlim = range(years.dat))
  grid(col = gray(0.7), lwd = 1)
  lines(years.dat ,z[ind,1], col = 'black', lwd = 2)
  polygon(c(years.dat,rev(years.dat)), c(z[ind,2],rev(z[ind,3])), col = tcol, border = "transparent", lty = 2)

  temp = cbind.data.frame(mat = gbcod.mat.rev$dat$Y, n = gbcod.mat.rev$dat$N, age = gbcod.mat.rev$dat$age_obs, 
    cohort = gbcod.mat.rev$dat$cy_obs + 1962, year = gbcod.mat.rev$dat$year_obs + 1962)
  x = t(sapply(1970:2014, function(y)
  {
    #x = summary(glm(mat ~ 1, family = binomial, data = temp, subset = age == 3 & year == 1975))$coef[1:2]
    x = try(summary(glm(cbind(mat,n-mat) ~ 1, family = binomial, data = temp, subset = age == whichage & year == y))$coef[1:2])
    if(!is.character(x))
    {
      x = c(1/(1 + exp(-x[1])), 1/(1 + exp(-(x[1] + c(-1,1)*qnorm(0.975)*x[2]))))
#      print(x)
      plotCI(x = y, y = x[1], li = x[2], ui = x[3], add = TRUE, slty = 2, sfrac = 0, lwd = 2)
      return(x)
    }
    else return(rep(NA,3))
  }))
}

par(mfcol = c(2,3), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp.fn(3, mod4.k.a50.AR)
temp.fn(2, mod4.k.a50.AR)
temp.fn(3, mod7.k.a50.AR)
temp.fn(2, mod7.k.a50.AR)
temp.fn(3, model = mod4.bb.k.a50.AR)
temp.fn(2, model = mod4.bb.k.a50.AR)

cairo_pdf('~/work/cod/ss_maturity/tex/maturity_at_age.pdf', family = "Times", height = 10, width = 7)
#png(filename = '~/work/cod/ss_maturity/tex/maturity_at_age.png', width = 7*144, height = 10*144, res = 144, pointsize = 12, family = "Times")#,
par(mfrow = c(2,1), mar = c(1,1,1,1), oma = c(4,4,0,0))
temp.fn(3, model = mod4.bb.k.a50.AR)
axis(2, lwd = 2, cex.axis = 1.5)
axis(1, labels = FALSE, cex.axis = 1.5)
box(lwd =2)
temp.fn(2, model = mod4.bb.k.a50.AR)
axis(2, lwd = 2, cex.axis = 1.5)
axis(1, lwd = 2, cex.axis = 1.5)
box(lwd =2)
#temp.fn(1)
#axis(2, lwd = 2, cex.axis = 1.5)
#axis(1, lwd = 2, cex.axis = 1.5)
#box(lwd =2)
mtext(side = 1, line = 2, "Year", cex = 2, outer = TRUE)
mtext(side = 2, line = 2, "Proportion mature", cex = 2, outer = TRUE)
dev.off()


