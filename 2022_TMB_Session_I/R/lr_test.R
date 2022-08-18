## Illustrate map feature of TMB to perform likelihood ratio tests on a ragged array dataset.
library(TMB)
compile("2022_TMB_Session_I/src/lr_test.cpp")
dyn.load(dynlib("2022_TMB_Session_I/src/lr_test"))

ngroup <- 5
nrep <- c(5,8,11,13,2)  ## Number of samples per group
mu <- rep(0,ngroup)     ## Mean value per group
sd <- c(1,1,1,2,2)      ## Standard deviation per group

## Simulate data
set.seed(123)
raggedArray <- lapply(1:ngroup,function(i)rnorm(nrep[i],mu[i],sd[i]))
str(raggedArray)

## Prepare data for TMB (ragged array not available):
obs <- unlist(raggedArray)
group <- factor( rep(1:length(raggedArray),sapply(raggedArray,length)) )

plot(group, obs)

## Both mu's and sd's un-restricted.
full.model <- MakeADFun(data=list(obs=obs,group=group),
                        parameters=list(mu=rep(0,ngroup),sd=rep(1,ngroup)),
                        DLL="lr_test"
                        )
full.model$par

## mu's restricted to be equal
map <- list(mu=factor(c(1,1,1,1,1)))
restricted.model1 <- MakeADFun(data=list(obs=obs,group=group),
                               parameters=list(mu=rep(0,ngroup),sd=rep(1,ngroup)),
                               DLL="lr_test",
                               map=map
                               )
restricted.model1$par

##Both mu's and sd's restricted to be equal
map <- list(mu=factor(c(1,1,1,1,1)),sd=factor(c(1,1,1,1,1)))
map <- list(mu = factor(c(NA, NA, NA, NA, NA)), sd = factor(c(NA, NA, NA, NA, NA)))
map <- list(mu=factor(c(1 ,2,3,4,5)),sd=factor(c(1,2,3,4,5)))
restricted.model2 <- MakeADFun(data=list(obs=obs,group=group),
                               parameters=list(mu=rep(0,ngroup),sd=rep(1,ngroup)),
                               DLL="lr_test",
                               map=map
                               )
restricted.model2$par

## Run models:
opt.full <- do.call("optim",full.model)
opt.r1 <- do.call("optim",restricted.model1)
opt.r2 <- do.call("optim",restricted.model2)

LRtest <- function(full,restricted){
    statistic <- 2*(restricted$value - full$value)
    df <- length(full$par) - length(restricted$par)
    p.value <- 1-pchisq(statistic,df=df)
    data.frame(statistic,df,p.value)
}

## H1, full: Both mu's and sd's un-restricted.
## H0, r1: mu's restricted to be equal, sd's un-restricted
LRtest(opt.full, opt.r1)

## H1, r1: mu's restricted to be equal, sd's un-restricted
## H0, r2: Both mu's and sd's restricted to be equal
LRtest(opt.r1,opt.r2)

opt.r1$par

