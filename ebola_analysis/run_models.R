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
library(ABSEIR)
library(coda)
library(splines)

# source helper functions
source('../helper_functions.R')
source('post_processing.R')

# matrix of 6 models to be run
modelsExpPS <- data.frame(infPeriodSpec = c('exp', 'PS'),
                         iddFun = NA,
                         stringsAsFactors = FALSE)
modelsIDD <- data.frame(infPeriodSpec = 'IDD',
                         iddFun = c('dgammaIDD', 'dlnormIDD', 'logitIDD', 'splineIDD'),
                         stringsAsFactors = FALSE)
allModels <- rbind.data.frame(modelsExpPS, modelsIDD)


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

# run three chains in parallel
cl <- makeCluster(3)
clusterExport(cl, list('ebolaDat',  'X', 'infPeriodSpec', 'iddFun', 'maxInf', 'idx'))

resThree <- parLapplyLB(cl, 1:3, function(x) {
    
    library(BayesSEIR)
    library(splines)
    
    # MCMC specifications
    niter <- 500000    # total number of iterations to be run
    nburn <- 100000     # number of burn-in iterations to be discarded     

    # get priors and initial values based on model/data generating scenario
    
    # set seed for reproducibility
    set.seed(x + idx)
    
    source('get_priors_inits.R')
    priorsInits <- get_priors_inits(infPeriodSpec = infPeriodSpec, 
                                    iddFun = iddFun, 
                                    maxInf = maxInf) 
    
    initsList<- priorsInits$initsList 
    priorList<- priorsInits$priorList 
    
    
    set.seed(x)
    
    if (infPeriodSpec == 'exp') {
        
        res <-  mcmcSEIR(dat = ebolaDat, X = X, 
                         inits = initsList, 
                         niter = niter, nburn = nburn,
                         infPeriodSpec = infPeriodSpec,
                         priors = priorList,
                         WAIC = TRUE)
        
    } else if (infPeriodSpec == 'PS') {
        
        res <- mcmcSEIR(dat = ebolaDat, X = X, 
                        inits = initsList, 
                        niter = niter, nburn = nburn,
                        infPeriodSpec = infPeriodSpec,
                        priors = priorList,
                        dist = 'gamma', maxInf = maxInf,
                        WAIC = TRUE)
        
    } else if (infPeriodSpec == 'IDD') {
        
        res <-  mcmcSEIR(dat = ebolaDat, X = X, 
                         inits = initsList, 
                         niter = niter, nburn = nburn,
                         infPeriodSpec = infPeriodSpec,
                         priors = priorList,
                         iddFun = iddFun, maxInf = maxInf,
                         WAIC = TRUE)
    }
    
    res
    
})
stopCluster(cl)

# get summaries from chains
postSummaries <- post_processing(modelOutput = resThree, 
                                 infPeriodSpec = infPeriodSpec, iddFun = iddFun, 
                                 maxInf = maxInf, X = X, N = N)


# save output in RDS form
saveRDS(postSummaries, paste0('./output/ebola_batch', idx, '.rds'))

