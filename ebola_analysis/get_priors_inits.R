################################################################################
# Function to set up priors and initial values depending on
#   model being fit
################################################################################

get_priors_inits <- function(infPeriodSpec, iddFun, maxInf) {
    
    
    # beta and rateE priors are the same for all models
    betaPrior <- function(x) {
        dnorm(x[1], mean = 0, sd = 4, log = T) + 
            dnorm(x[2], mean = 0, sd = 10, log = T)
    }
    
    rateEPrior <- function(x) {
        dgamma(x, 11, 100, log = T)
    }
    
    # beta and rateE initial values:
    # initial value on beta2 negative to reflect expected intervention effect
    betaInit <- c(rnorm(1, 0, 1),
                  runif(1, -20, -10))
    
    rateEInit <- rgamma(1, 11, 100)
    
    
    if (infPeriodSpec == 'exp') {
        # exponential infectious period
        
        rateIPrior <- function(x) {
            dgamma(x, 16, 100, log = T)
        }
        
        rateIInit <- rgamma(1, 16, 100)
        
        priorList <- list(betaPrior = betaPrior,
                          rateEPrior = rateEPrior,
                          rateIPrior = rateIPrior)
        
        initsList <- list(beta = round(betaInit, 8), 
                          rateE = round(rateEInit, 8), 
                          rateI = round(rateIInit, 8))
        
    } else if (infPeriodSpec == 'PS') {
        # path-specific infectious period
        
        psParamsPrior <- function(x) {
            dgamma(x['shape'], 180, 10, log = T) +
                dgamma(x['rate'], 30, 10, log = T)
        }
        psParamsInit <- list(shape = rgamma(1, 180, 10), 
                             rate = rgamma(1, 30, 10))
        
        
        priorList <- list(betaPrior = betaPrior,
                          rateEPrior = rateEPrior,
                          psParamsPrior = psParamsPrior)
        
        initsList <- list(beta = round(betaInit, 8), 
                          rateE = round(rateEInit, 8), 
                          psParams = lapply(psParamsInit, round, 8))
        
    } else if (infPeriodSpec == 'IDD') {
        # IDD transmissibility 
        
        if (iddFun == 'dgammaIDD') {
            
            iddParamsPrior <- function(x) {
                dgamma(x['shape'], 1, 1, log = T) +
                    dgamma(x['rate'], 1, 1, log = T)
            }
            
            iddParamsInit = list(shape = runif(1, 0.2, 5),
                                 rate = runif(1, 0.1, 2))
            
        } else if (iddFun == 'dlnormIDD') {
            
            iddParamsPrior <- function(x) {
                dnorm(x['meanlog'], log(2), 0.2, log = T) +
                    dgamma(x['sdlog'], 1, 1, log = T)
            }
            
            iddParamsInit = list(meanlog = rnorm(1, log(2), 0.2),
                                 sdlog = rgamma(1, 1, 1))
            
        } else if (iddFun == 'logitIDD') {
            
            iddParamsPrior <- function(x) {
                dnorm(x['mid'], 6, 1, log = T) +
                    dgamma(x['rate'], 1, 1, log = T)
            }
            
            iddParamsInit = list(mid = rnorm(1, 6, 1),
                                 rate = runif(1, 0.5, 2))
            
        } else if (iddFun == 'splineIDD') {
            
            iddParamsPrior <- function(x) {
                sum(dnorm(x[1:2], 0, 2, log = T)) +
                    dnorm(x[3], -2, 2, log = T) +
                    dnorm(x[4], -3, 2, log = T) +
                    dnorm(x[5], -4, 2, log = T)
            }
            
            iddParamsInit = list(b1 = rnorm(1, 0, 2),
                                 b2 = rnorm(1, 0, 2),
                                 b3 = rnorm(1, -2, 2),
                                 b4 = rnorm(1, -3, 2),
                                 b5 = rnorm(1, -4, 2),
                                 XBasis = bs(1:maxInf, df = 5))
        } 
        
        priorList <- list(betaPrior = betaPrior,
                          rateEPrior = rateEPrior,
                          iddParamsPrior = iddParamsPrior)
        
        initsList <- list(beta = round(betaInit, 8), 
                          rateE = round(rateEInit, 8), 
                          iddParams = lapply(iddParamsInit, round, 8))
        
    }
    
    list(priorList = priorList,
         initsList= initsList)
    
    
}
