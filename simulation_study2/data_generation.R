################################################################################
# Data generation
#
################################################################################

# load packages
library(BayesSEIR)
library(truncdist)

# source helper functions
source('../helper_functions.R')

# simulation characteristics

nSim <- 100

# parameters for all
N <- 5363500
E0 <- 1
I0 <- 0
S0 <- N - E0 - I0

tau <- 250
tstar <- 120
X <- getX(tau, tstar)
rateE <- 0.2
maxInf <- 15


################################################################################
# Path specific model

dist <- 'gamma'
psParams <- list(shape = 56, rate = 7)
betaPS <- c(-1.77, -7)

# want to save Estar and Istar from the simulated epidemics
store_PS_epidemics <- list(Estar = matrix(NA, nrow = tau, ncol = nSim),
                           Istar = matrix(NA, nrow = tau, ncol = nSim))

# set seed for reproducibility
set.seed(1)

# simulate until nSim epidemics meet the conditions
i <- 1
repeat {
    
    sim_epi <- simSEIR(S0 = S0, E0 = E0, I0 = I0, N = N, tau = tau,
                      beta = betaPS, X = X, rateE = rateE,
                      infPeriodSpec = 'PS',
                      infPSParams = list(dist = dist,
                                         psParams = psParams,
                                         maxInf = maxInf))
    
    # only keep realistic epidemics (> 250 cases lasting at least 140 days)
    numInf <- sum(sim_epi$Istar)
    epiLength <- max(which(sim_epi$Istar > 0))
    epiOver <- all(sim_epi$Istar[(tau-20):tau] == 0) 
    
    
    if (numInf > 250 & epiLength > 140 & epiOver) {
        store_PS_epidemics$Estar[,i] <- sim_epi$Estar
        store_PS_epidemics$Istar[,i] <- sim_epi$Istar
        i <- i + 1
    }
    
    if (i == nSim + 1) break
}

# save
saveRDS(store_PS_epidemics, 'data/PS_data.rds')

################################################################################
# IDD Peak
################################################################################

iddFun <- dgammaIDD
iddParams <- list(shape = 4, rate = 1)
betaIDDPeak <- c(0.25, -7)

# want to save Estar and Istar from the simulated epidemics
store_IDD_peak_epidemics <- list(Estar = matrix(NA, nrow = tau, ncol = nSim),
                           Istar = matrix(NA, nrow = tau, ncol = nSim))

# set seed for reproducibility
set.seed(1)

# simulate until nSim epidemics meet the conditions
i <- 1
repeat {
    
    sim_epi <- simSEIR(S0 = S0, E0 = E0, I0 = I0, N = N, tau = tau,
                       beta = betaIDDPeak, X = X, rateE = rateE,
                       infPeriodSpec = 'IDD',
                       infIDDParams = list(iddFun = iddFun,
                                           iddParams = iddParams,
                                           maxInf = maxInf))
    
    # only keep realistic epidemics (> 250 cases lasting between 140 - 230 days)
    numInf <- sum(sim_epi$Istar)
    epiLength <- max(which(sim_epi$Istar > 0))
    epiOver <- all(sim_epi$Istar[(tau-20):tau] == 0) 
    
    
    if (numInf > 250 & epiLength > 140 & epiOver) {
        store_IDD_peak_epidemics$Estar[,i] <- sim_epi$Estar
        store_IDD_peak_epidemics$Istar[,i] <- sim_epi$Istar
        i <- i + 1
    }
    
    if (i == nSim + 1) break
}

# save
saveRDS(store_IDD_peak_epidemics, 'data/IDD_peak_data.rds')

################################################################################
# IDD Exp
################################################################################

iddFun <- dgammaIDD
iddParams <- list(shape = 0.9, rate = 0.3)
betaIDDExp <- c(0.4, -7)

# want to save Estar and Istar from the simulated epidemics
store_IDD_exp_epidemics <- list(Estar = matrix(NA, nrow = tau, ncol = nSim),
                                 Istar = matrix(NA, nrow = tau, ncol = nSim))

# set seed for reproducibility
set.seed(1)

# simulate until nSim epidemics meet the conditions
i <- 1
repeat {
    
    sim_epi <- simSEIR(S0 = S0, E0 = E0, I0 = I0, N = N, tau = tau,
                       beta = betaIDDExp, X = X, rateE = rateE,
                       infPeriodSpec = 'IDD',
                       infIDDParams = list(iddFun = iddFun,
                                           iddParams = iddParams,
                                           maxInf = maxInf))
    
    # only keep realistic epidemics (> 250 cases lasting between 140 - 230 days)
    numInf <- sum(sim_epi$Istar)
    epiLength <- max(which(sim_epi$Istar > 0))
    epiOver <- all(sim_epi$Istar[(tau-20):tau] == 0) 
    
    
    if (numInf > 250 & epiLength > 140 & epiOver) {
        store_IDD_exp_epidemics$Estar[,i] <- sim_epi$Estar
        store_IDD_exp_epidemics$Istar[,i] <- sim_epi$Istar
        i <- i + 1
    }
    
    if (i == nSim + 1) break
}

# save
saveRDS(store_IDD_exp_epidemics, 'data/IDD_exp_data.rds')

################################################################################
# IDD Logit
################################################################################

iddFun <- logitIDD
iddParams <- list(rate = 1.5, mid = 8)
betaIDDLogit <- c(-1.77, -7)

# want to save Estar and Istar from the simulated epidemics
store_IDD_logit_epidemics <- list(Estar = matrix(NA, nrow = tau, ncol = nSim),
                                 Istar = matrix(NA, nrow = tau, ncol = nSim))

# set seed for reproducibility
set.seed(1)

# simulate until nSim epidemics meet the conditions
i <- 1
repeat {
    
    sim_epi <- simSEIR(S0 = S0, E0 = E0, I0 = I0, N = N, tau = tau,
                       beta = betaIDDLogit, X = X, rateE = rateE,
                       infPeriodSpec = 'IDD',
                       infIDDParams = list(iddFun = iddFun,
                                           iddParams = iddParams,
                                           maxInf = maxInf))
    
    # only keep realistic epidemics (> 250 cases lasting between 140 - 230 days)
    numInf <- sum(sim_epi$Istar)
    epiLength <- max(which(sim_epi$Istar > 0))
    epiOver <- all(sim_epi$Istar[(tau-20):tau] == 0) 
    
    
    if (numInf > 250 & epiLength > 140 & epiOver) {
        store_IDD_logit_epidemics$Estar[,i] <- sim_epi$Estar
        store_IDD_logit_epidemics$Istar[,i] <- sim_epi$Istar
        i <- i + 1
    }
    
    if (i == nSim + 1) break
}

# save
saveRDS(store_IDD_logit_epidemics, 'data/IDD_logit_data.rds')
