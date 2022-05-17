################################################################################
# Model fitting IDD transmissibility models with known exposure times
# Using four IDD model fits:
#   gamma pdf       - fit to IDD Peak and IDD Exp
#   log-normal pdf  - fit to IDD Peak and IDD Exp
#   logistic decay  - fit to IDD Logit
#   basis splines   - fit to all IDD data generation
# Using two maximum duration of infectious period:
#   15 days (truth)
#   20 days (misspecification)
# Run three chains in parallel
# For each model/simulation, want:
#   Gelman Rubin to ensure convergence
#   Posterior median estimated IDD curve
################################################################################

### take in array parameter from HPC batch job
args <- commandArgs(trailingOnly=TRUE)
idx <- gsub('\r', '', args)
idx <- as.numeric(idx)

### load libraries
library(parallel)
library(coda)
library(splines)
library(BayesSEIR)

source('../helper_functions.R')
source('post_processing.R')

# create data frame of all possible models to be fit

maxInfs <- c(15, 20)
nSim <- 100

modelsGamma <- expand.grid(iddFun = 'dgammaIDD',
                           simNumber = 1:nSim,
                           maxInf = maxInfs,
                           datGen = c('IDD_peak', 'IDD_exp'),
                           stringsAsFactors = FALSE)

modelsLognormal <- expand.grid(iddFun = 'dlnormIDD',
                               simNumber = 1:nSim,
                               maxInf = maxInfs,
                               datGen = c('IDD_peak', 'IDD_exp'),
                               stringsAsFactors = FALSE)

modelsLogit <- expand.grid(iddFun = 'logitIDD',
                           simNumber = 1:nSim,
                           maxInf = maxInfs,
                           datGen = 'IDD_logit',
                           stringsAsFactors = FALSE)

modelsSpline <- expand.grid(iddFun = 'splineIDD',
                            simNumber = 1:nSim,
                            maxInf = maxInfs,
                            datGen = c('IDD_peak', 'IDD_exp', 'IDD_logit'),
                            stringsAsFactors = FALSE)

# 1600 models to be fit
allModels <- rbind.data.frame(modelsGamma,
                              modelsLognormal,
                              modelsLogit,
                              modelsSpline)

# model specifications that are the same for all models
N <- 5363500
E0 <- 1
I0 <- 0
S0 <- N - E0 - I0

# intervention time
tstar <- 120

# fit models in batches of 100 (16 batches total)
# each batch is one model/data generation scenario
batchSize <- 100
batchIdx <- batchSize * (idx - 1) + 1:batchSize


for (i in batchIdx) {
    
    iddFun_i <- allModels$iddFun[i]
    simNumber_i <- allModels$simNumber[i]
    maxInf_i <- allModels$maxInf[i]
    datGen_i <- allModels$datGen[i]
    
    print(paste0('IDD Fun: ', iddFun_i,
                 ', data gen: ', datGen_i,
                 ', max inf: ', maxInf_i,
                 ', sim number: ', simNumber_i))
    
    ############################################################################
    ### set up data
    
    # load data
    dat <- readRDS(paste0('data/', datGen_i, '_data.rds'))
    
    # extract data for simulation of interest
    Istar <- dat$Istar[,simNumber_i]
    Estar <- dat$Estar[,simNumber_i]
    
    # trim/add to Istar and Estar in the case of excess 0's
    fullTime <- length(Istar)
    lastInfTime <- max(which(Istar > 0))
    if (lastInfTime + maxInf_i <= fullTime) {
        
        newTime <- lastInfTime + maxInf_i
        
        Istar <- Istar[1:newTime]
        Estar <- Estar[1:newTime]
        
    } else {
        
        zerosAdd <- lastInfTime + maxInf_i - fullTime
        
        Istar <- c(Istar, rep(0, zerosAdd))
        Estar <- c(Estar, rep(0, zerosAdd))
        newTime <- length(Istar)
        
    }
    
    # design matrix for intervention
    X <- getX(newTime, tstar)
    
    datList <- list(Istar = Istar,
                    Estar = Estar,
                    S0 = S0,
                    E0 = E0,
                    I0 = I0,
                    N = N)
    
    ############################################################################
    ### run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('datList',  'X', 'iddFun_i', 'datGen_i', 
                           'maxInf_i', 'i'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(BayesSEIR)
        
        # MCMC specifications
        niter <- 1e6
        nburn <- 2e5
        
        # set seed for reproducibility of initial values
        set.seed(x + i)
        
        # get priors and initial values based on model/data generating scenario
        source('get_priors_inits.R')
        priorsInits <- get_priors_inits(infPeriodSpec = 'IDD', 
                                        iddFun = iddFun_i, 
                                        datGen = datGen_i, 
                                        maxInf = maxInf_i) 
        
        initsList<- priorsInits$initsList 
        priorList<- priorsInits$priorList 
        
        set.seed(x)
        fullPost <- mcmcSEIR(dat = datList, X = X, 
                             inits = initsList, 
                             niter = niter, nburn = nburn,
                             infPeriodSpec = 'IDD',
                             priors = priorList,
                             iddFun = iddFun_i, maxInf = maxInf_i,
                             EKnown = TRUE)
        
        list(fullPost = fullPost)
    })
    stopCluster(cl)
    
    ### post processing to get Gelman-Rubin and IDD curve estimates
    postSummaries <- post_processing(modelOutput = resThree, EType = 'known',
                                     infPeriodSpec = "IDD", 
                                     datGen = datGen_i, iddFun = iddFun_i, 
                                     simNumber = simNumber_i, maxInf = maxInf_i,
                                     X = X, N = N)
    
    iddSummary <- postSummaries$iddSummary
    iddSummary$allConverge <- all(postSummaries$gdiag$gr < 1.1)
    
    # concatenate results across batches for output
    if (i == batchIdx[1]) {
        batchOutput <- iddSummary
    } else {
        batchOutput <- rbind.data.frame(batchOutput, iddSummary)
    }
    
} # end loop

# save output in RDS form
saveRDS(batchOutput, paste0('./batch_output/knownE_batch', idx, '.rds'))

