library(TMB)

compile('2022_TMB_Session_I/src/linReg.cpp',"-O1 -g",DLLFLAGS="")
dyn.load(dynlib('2022_TMB_Session_I/src/linReg'))
Data <- list(y = c(-2.1, 3.3, 4.2),
             X = cbind(c(1,1,1),c(4,5,6))
)

Pars <- list(beta = c(0,0), lnSigma = 0)
obj <- MakeADFun(data = Data, 
                 parameters = Pars,  
                 DLL = "linReg")