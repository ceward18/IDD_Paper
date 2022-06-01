################################################################################
# time per iteration for IDD, Exp, PS for varying sizes of epidemics
# takes ~ 
################################################################################

library(BayesSEIR)
library(microbenchmark)

# initial conditions for all epidemics
I0 <- 1
E0 <- 0
tau <- 200
X <- matrix(rep(1, tau), ncol = 1)

# parameters used to generate highly transmissible epidemics
beta <- log(0.7)
rateE <- 1/7
rateI <- 1/7

# various population sizes to be considered
NVec <- c(200, 2000, 5000, 10000, 15000, 25000, 50000)

# MCMC settings (low iterations because just running for time per iter)
niter <- 100
nburn <- 80

# priors that will be consistent for all models
betaPrior <- function(x) {
    dnorm(x, mean = 0, sd = 4, log = T) 
}

rateEPrior <- function(x) {
    dgamma(x, 14, 100, log = T)
}

# exponential model specific priors
rateIPrior <- function(x) {
    dgamma(x, 14, 100, log = T)
}

priorListExp <- list(betaPrior=betaPrior,
                     rateEPrior=rateEPrior,
                     rateIPrior=rateIPrior)

# IDD specific set-up
iddParamsPrior <- function(x) {
    dgamma(x['shape'], 1, 1, log = T) +
        dgamma(x['rate'], 1, 1, log = T)
}

priorListIDD <- list(betaPrior=betaPrior,
                     rateEPrior=rateEPrior,
                     iddParamsPrior=iddParamsPrior)

iddFun <- 'dgammaIDD'
maxInf <- 16

# PS specific set-up
psParamsPrior <- function(x) {
    dgamma(x['shape'], 700, 10, log = T) +
        dgamma(x['rate'], 100, 10, log = T)
}

priorListPS <- list(betaPrior=betaPrior,
                    rateEPrior=rateEPrior,
                    psParamsPrior=psParamsPrior)

dist <- 'gamma'


################################################################################
# loop over each epidemic size

# record total number infected from simulation and compute time for each model
nInf <- rep(NA, length(NVec))
expTime <- iddTime <- psTime <- rep(NA, length(NVec))

for (i in 1:length(NVec)) {
    
    print(i)
    
    N <- NVec[i]
    S0 <- N - E0 - I0
    
    set.seed(i)
    repeat {
        # simulate data using exponential distribution
        simDat <- simSEIR(S0, E0, I0, N, tau,
                          beta, X, rateE,
                          infPeriodSpec = 'exp',
                          infExpParams = list(rateI = rateI))
        
        epiOver <- max(which(simDat$Istar > 0)) < (tau - maxInf - 1)
        
        if ((sum(simDat$Istar) > 0.8 * N) & epiOver) break
    }
    
    nInf[i] <- sum(simDat$Istar)
    
    
    ### fit exponential model
    print('Starting exponential model')
    initsList <- list(beta = rnorm(1, 0, 2),  
                      rateE = rgamma(1, 14, 100),
                      rateI = rgamma(1, 14, 100))
    
    expBM <- microbenchmark(mcmcSEIR(dat = simDat, X = X, 
                                     inits = initsList, 
                                     niter = niter, nburn = nburn,
                                     infPeriodSpec = 'exp',
                                     priors = priorListExp,
                                     EKnown = T), times = 3)
    
    # time is expressed in nanoseconds - convert to seconds
    expTime[i] <-  median(expBM$time / 1e9) / niter
    
    
    ### fit IDD model (gamma pdf)
    print('Starting IDD model')
    initsList <- list(beta = rnorm(1, 0, 2), 
                      rateE = rgamma(1, 14, 100),
                      iddParams = list(shape = rgamma(1, 1, 1),
                                       rate = rgamma(1, 1, 1)))
    
    iddBM <- microbenchmark(mcmcSEIR(dat = simDat, X = X, 
                                     inits = initsList, 
                                     niter = niter, nburn = nburn,
                                     infPeriodSpec = 'IDD',
                                     priors = priorListIDD,
                                     iddFun = iddFun, maxInf = maxInf,
                                     EKnown = T), times = 3)
    
    # time is expressed in nanoseconds - convert to seconds
    iddTime[i] <- median(iddBM$time / 1e9) / niter
    
    
    ### fit PS model
    print('Starting PS model')
    initsList <- list(beta = rnorm(1, 0, 2),
                      rateE = rgamma(1, 14, 100),
                      psParams = list(shape = rgamma(1, 700, 10),
                                      rate = rgamma(1, 100, 10)))
    
    psBM <- microbenchmark(mcmcSEIR(dat = simDat, X = X,
                                    inits = initsList,
                                    niter = niter, nburn = nburn,
                                    infPeriodSpec = 'PS',
                                    priors = priorListPS,
                                    dist = dist, maxInf = maxInf,
                                    EKnown = T), times = 3)
    
    # time is expressed in nanoseconds - convert to seconds
    psTime[i] <-  median(psBM$time / 1e9) / niter
  
    
} # end loop

fullRes <- cbind.data.frame(nInf=rep(nInf, 3),
                            sec = c(expTime, iddTime, psTime),
                            model = rep(c('exp', 'IDD', 'PS'), each = length(nInf)))

saveRDS(fullRes, 'timePerIter.rds')















