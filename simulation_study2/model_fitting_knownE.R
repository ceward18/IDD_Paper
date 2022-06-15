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
    ### run chains until convergence
    postSummaries <- runModels(datList = datList, X = X, infPeriodSpec_i = 'IDD', 
                               iddFun_i = iddFun_i, datGen_i = datGen_i, 
                               maxInf_i = maxInf_i, EType = 'known', i = i, 
                               niter = 5e5)

    
    # concatenate results across batches for output
    if (i == batchIdx[1]) {
        gdiag <- postSummaries$gdiag
        paramsSummary <- postSummaries$paramsSummary
        iddSummary <- postSummaries$iddSummary
        r0Summary <- postSummaries$r0Summary
        mcmcEffSummary <- postSummaries$mcmcEffSummary
    } else {
        gdiag <- rbind.data.frame(gdiag, postSummaries$gdiag)
        paramsSummary <- rbind.data.frame(paramsSummary, postSummaries$paramsSummary)
        iddSummary <- rbind.data.frame(iddSummary, postSummaries$iddSummary)
        r0Summary <- rbind.data.frame(r0Summary, postSummaries$r0Summary)
        mcmcEffSummary <- rbind.data.frame(mcmcEffSummary, postSummaries$mcmcEffSummary)
    }
    
} # end loop

# combine results back together
batchOutput <- list(gdiag = gdiag,
                    paramsSummary = paramsSummary,
                    iddSummary = iddSummary,
                    r0Summary = r0Summary,
                    mcmcEffSummary = mcmcEffSummary)

# save output in RDS form
saveRDS(batchOutput, paste0('./batch_output/knownE_batch', idx, '.rds'))

