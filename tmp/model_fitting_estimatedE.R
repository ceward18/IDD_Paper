################################################################################
# Model fitting IDD transmissibility models with unknown exposure times
#   (MAIN RESULTS)
# Using six model fits:
#   exponential infectious period
#   path-specific infections period
#   gamma pdf       
#   log-normal pdf 
#   logistic decay  
#   basis splines   
# All six model fits are fit to all four data generating scenarios
# Using two maximum duration of infectious period:
#   15 days (truth)
#   20 days (miss-specification)
# Run three chains in parallel
# For each model/simulation, want:
#   Gelman Rubin to ensure convergence
#   Posterior median estimated IDD curve
#   Posterior mean estimated R0(t)
#   MSE of R0
#   MCMC Efficiency
################################################################################

### take in array parameter from HPC batch job
args <- commandArgs(trailingOnly=TRUE)
idx <- gsub('\r', '', args)
idx <- as.numeric(idx)

if (file.exists(paste0('./batch_output/estimatedE_batch', idx, '.rds'))) {
  stop('this already ran!')
}

### load libraries
library(parallel)
library(coda)
library(splines)
library(BayesSEIR)

source('../helper_functions.R')
source('post_processing.R')

datGens <- c('PS', 'IDD_peak', 'IDD_exp', 'IDD_logit')
maxInfs <- c(15, 20)
nSim <- 100
chains <- 1:3

modelsPS <- expand.grid(datGen = datGens,
                        infPeriodSpec = 'PS',
                        iddFun = NA,
                        maxInf = maxInfs,
                        simNumber = 1:nSim,
                        chain = chains,
                        stringsAsFactors = FALSE)
modelsPS <- modelsPS[order(modelsPS$datGen, modelsPS$infPeriodSpec,
                           modelsPS$iddFun, modelsPS$maxInf, modelsPS$simNumber),]


# 4800 models to run (6 models * 4 data generation * 2 maxInfs * 100 sims)
# IDD 1 - 3200
# Exp 3201 - 4000
# PS 4001 - 4800
allModels <- modelsPS
rownames(allModels) <- NULL

# model specifications that are the same for all models
N <- 5363500
E0 <- 1
I0 <- 0
S0 <- N - E0 - I0

# intervention time
tstar <- 120

infPeriodSpec_i <- allModels$infPeriodSpec[idx]
iddFun_i <- allModels$iddFun[idx]
simNumber_i <- allModels$simNumber[idx]
maxInf_i <- allModels$maxInf[idx]
datGen_i <- allModels$datGen[idx]
x <- allModels$chain[idx]
i <- ceiling(idx/3)

print(paste0('Model: ', infPeriodSpec_i,
             ', IDD Fun: ', iddFun_i,
             ', data gen: ', datGen_i,
             ', max inf: ', maxInf_i,
             ', sim number: ', simNumber_i))

############################################################################
### set up data

# load data
dat <- readRDS(paste0('data/', datGen_i, '_data.rds'))

# extract data for simulation of interest
Istar <- dat$Istar[,simNumber_i]

# trim/add to Istar and Estar in the case of excess 0's
# need maxInf days after last infection start
fullTime <- length(Istar)
lastInfTime <- max(which(Istar > 0))
if (lastInfTime + maxInf_i <= fullTime) {
  
  newTime <- lastInfTime + maxInf_i
  
  Istar <- Istar[1:newTime]
  
} else {
  
  zerosAdd <- lastInfTime + maxInf_i - fullTime
  
  Istar <- c(Istar, rep(0, zerosAdd))
  newTime <- length(Istar)
  
}

# design matrix for intervention
X <- getX(newTime, tstar)

datList <- list(Istar = Istar,
                S0 = S0,
                E0 = E0,
                I0 = I0,
                N = N)

############################################################################

# MCMC specifications
# total number of iterations to be run
if (infPeriodSpec_i == 'exp') {
  niter <- 1e6 
} else if (infPeriodSpec_i == 'PS') {
  niter <- 5e5 
} else if (infPeriodSpec_i == 'IDD') {
  niter <- 1.3e6 
} 

# number of burn-in iterations to be discarded     
nburn <- 1e5

# set seed for reproducibility of initial values
set.seed(x + i)

# get priors and initial values based on model/data generating scenario
source('get_priors_inits.R')
priorsInits <- get_priors_inits(infPeriodSpec = infPeriodSpec_i, 
                                iddFun = iddFun_i, 
                                datGen = datGen_i, 
                                maxInf = maxInf_i) 

initsList<- priorsInits$initsList 
priorList<- priorsInits$priorList 

# set seed for reproducibility of chains
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

res$chainTime <- as.numeric(endTime - startTime, units = 'mins')

############################################################################
# post processing of the model output

# postSummaries <- post_processing(modelOutput = resThree, EType = 'estimated',
#                                  infPeriodSpec = infPeriodSpec_i, 
#                                  datGen = datGen_i, iddFun = iddFun_i, 
#                                  simNumber = simNumber_i, maxInf = maxInf_i,
#                                  X = X, N = N)

# save output in RDS form
saveRDS(res, paste0('./batch_output/estimatedE_batch', idx, '.rds'))










