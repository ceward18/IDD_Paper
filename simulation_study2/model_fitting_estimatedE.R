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


modelsIDD <- expand.grid(datGen = datGens,
                         infPeriodSpec = 'IDD',
                         iddFun = c('dgammaIDD', 'dlnormIDD', 'logitIDD', 'splineIDD'),
                         maxInf = maxInfs,
                         simNumber = 1:nSim,
                         stringsAsFactors = FALSE)
modelsIDD <- modelsIDD[order(modelsIDD$datGen, modelsIDD$infPeriodSpec,
                             modelsIDD$iddFun, modelsIDD$maxInf, modelsIDD$simNumber),]


modelsExp <- expand.grid(datGen = datGens,
                         infPeriodSpec = 'exp',
                         iddFun = NA,
                         maxInf = 20,
                         simNumber = 1:nSim,
                         stringsAsFactors = FALSE)
modelsExp <- modelsExp[order(modelsExp$datGen, modelsExp$infPeriodSpec,
                             modelsExp$iddFun, modelsExp$maxInf, modelsExp$simNumber),]


modelsPS <- expand.grid(datGen = datGens,
                        infPeriodSpec = 'PS',
                        iddFun = NA,
                        maxInf = maxInfs,
                        simNumber = 1:nSim,
                        stringsAsFactors = FALSE)
modelsPS <- modelsPS[order(modelsPS$datGen, modelsPS$infPeriodSpec,
                           modelsPS$iddFun, modelsPS$maxInf, modelsPS$simNumber),]


# 4400 models to run 
# IDD 1 - 3200
# Exp 3201 - 3600
# PS 3601 - 4400
allModels <- rbind.data.frame(modelsIDD, modelsExp, modelsPS)

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

################################################################################
### run chains until convergence

postSummaries <- runModels(datList = datList, X = X, 
                           infPeriodSpec_i = infPeriodSpec_i, 
                           iddFun_i = iddFun_i, datGen_i = datGen_i, 
                           maxInf_i = maxInf_i, EType = 'estimated', i = idx, 
                           niter = 4e5)

# save output in RDS form
saveRDS(postSummaries, paste0('./batch_output/estimatedE_batch', idx, '.rds'))










