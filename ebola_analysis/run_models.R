################################################################################
# kikwit analysis
# fit 6 models to ebola data
#   exponential infectious period
#   path-specific infectious period
#   IDD: gamma pdf
#   IDD: log normal pdf
#   IDD: logistic decay
#   IDD: basis splines
################################################################################

args <- commandArgs(trailingOnly=TRUE)
idx <- gsub('\r', '', args)
idx <- as.numeric(idx)

### load libraries
library(BayesSEIR)
library(parallel)

# source helper functions
source('../helper_functions.R')

# matrix of 6 models to be run
modelsExpPS <- data.frame(infPeriodSpec = c('exp', 'PS'),
                         iddFun = NA,
                         stringsAsFactors = FALSE)
modelsIDD <- data.frame(infPeriodSpec = 'IDD',
                         iddFun = c('dgammaIDD', 'dlnormIDD', 'logitIDD', 'splineIDD'),
                         stringsAsFactors = FALSE)
allModels <- rbind.data.frame(modelsExpPS, modelsIDD)


# 4800 models to run (6 models * 4 data generation * 2 maxInfs * 100 sims)
# IDD 1 - 3200
# Exp 3201 - 4000
# PS 4001 - 4800
allModels <- rbind.data.frame(modelsIDD, modelsExp, modelsPS)


# load data from ABSEIR package
# devtools::install_git("https://github.com/grantbrown/ABSEIR.git")
ebola <- ABSEIR::Kikwit1995

# maximum infectious period length
maxInf <- 21

# trim ebola data to only include maxInf days after last infection time
lastInfTime <- max(which(ebola$Count > 0))
newTime <- lastInfTime + maxInf
    
ebola <- ebola[1:newTime,]

# format data and initial conditions for modeling
N <- 5363500
E0 <- 3
I0 <- 0
S0 <- N - E0 - I0

ebolaDat <- list(Istar = ebola$Count[-c(1:3)],
                 S0 = S0,
                 E0 = E0,
                 I0 = I0,
                 N = N)

# format matrix for the intensity process
tau <- length(ebolaDat$Istar)
interventionTime <- which(ebola$Date[-c(1:3)] == as.Date('1995-05-10'))

X <- cbind(1, cumsum(1:tau > interventionTime) / 100)


# obtain model specifications by array parameter
infPeriodSpec <- allModels$infPeriodSpec[idx]
iddFun <- allModels$iddFun[idx]


# get priors and initial values based on model/data generating scenario
priorsInits <- get_priors_inits(infPeriodSpec = infPeriodSpec, 
                                iddFun = iddFun, 
                                maxInf = maxInf) 

initsList<- priorsInits$initsList 
priorList<- priorsInits$priorList 

# run three chains in parallel
cl <- makeCluster(3)
clusterExport(cl, list('datList',  'X', 'initsList',
                       'priorList', 'infPeriodSpec_i', 'iddFun_i', 'maxInf_i'))

resThree <- parLapplyLB(cl, 1:3, function(x) {
    
    library(BayesSEIR)
    
    # MCMC specifications
    # total number of iterations to be run
    niter <- 500000
    
    # number of burn-in iterations to be discarded     
    nburn <- 50000
    
    set.seed(x)
    
    # start timing for MCMC efficiency
    startTime <- Sys.time()
    
    if (infPeriodSpec_i == 'exp') {
        
        res <-  mcmcSEIR(dat = datList, X = X, 
                         inits = initsList, 
                         niter = niter, nburn = nburn,
                         infPeriodSpec = infPeriodSpec_i,
                         priors = priorList,
                         WAIC = TRUE)
        
    } else if (infPeriodSpec_i == 'PS') {
        
        res <- mcmcSEIR(dat = datList, X = X, 
                        inits = initsList, 
                        niter = niter, nburn = nburn,
                        infPeriodSpec = infPeriodSpec_i,
                        priors = priorList,
                        dist = 'gamma', maxInf = maxInf_i,
                        WAIC = TRUE)
        
    } else if (infPeriodSpec_i == 'IDD') {
        
        res <-  mcmcSEIR(dat = datList, X = X, 
                         inits = initsList, 
                         niter = niter, nburn = nburn,
                         infPeriodSpec = infPeriodSpec_i,
                         priors = priorList,
                         iddFun = iddFun_i, maxInf = maxInf_i,
                         WAIC = TRUE)
    }
    
    endTime <- Sys.time()
    
    res
    
})
stopCluster(cl)
























