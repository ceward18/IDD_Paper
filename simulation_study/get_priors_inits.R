################################################################################
# Function to set up priors and initial values depending on
#   model being fit
#   data generating scenario (only impacts log-normal model)
################################################################################

get_priors_inits <- function(infPeriodSpec, iddFun, datGen, maxInf) {
    
    
    # beta and rateE priors are the same for all models
    betaPrior <- function(x) {
        sum(dnorm(x, mean = 0, sd = 4, log = T)) 
    }
    
    rateEPrior <- function(x) {
        dgamma(x, 20, 100, log = T)
    }
    
    # beta and rateE initial values:
    # initial value on beta2 negative to reflect expected intervention effect
    betaInit <- c(rnorm(1, 0, 1),
                  runif(1, -10, -4))
    
    rateEInit <- rgamma(1, 20, 100)
    
    
    if (infPeriodSpec == 'exp') {
        # exponential infectious period
        
        rateIPrior <- function(x) {
            dgamma(x, 2, 16, log = T)
        }
        
        rateIInit <- rgamma(1, 2, 16)
        
        priorList <- list(betaPrior = betaPrior,
                          rateEPrior = rateEPrior,
                          rateIPrior = rateIPrior)
        
        initsList <- list(beta = round(betaInit, 8), 
                          rateE = round(rateEInit, 8), 
                          rateI = round(rateIInit, 8))
        
    } else if (infPeriodSpec == 'PS') {
        # path-specific infectious period
        
        psParamsPrior <- function(x) {
            dgamma(x['shape'], 560, 10, log = T) +
                dgamma(x['rate'], 70, 10, log = T)
        }
        psParamsInit <- list(shape = rgamma(1, 560, 10), 
                             rate = rgamma(1, 70, 10))
        
        
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
            
            if (datGen == 'IDD_peak') {
                
                iddParamsPrior <- function(x) {
                    dnorm(x['meanlog'], log(4), 0.2, log = T) +
                        dgamma(x['sdlog'], 1, 1, log = T)
                }
                
                iddParamsInit = list(meanlog = rnorm(1, log(4), 0.2),
                                     sdlog = rgamma(1, 1, 1))
                
            } else if (datGen %in% c('IDD_exp', 'IDD_logit', 'PS')) {
                
                iddParamsPrior <- function(x) {
                    dgamma(x['meanlog'], log(3), 0.2, log = T) +
                        dgamma(x['sdlog'], 1, 1, log = T)
                }
                
                iddParamsInit = list(meanlog = rnorm(1, log(3), 0.2),
                                     sdlog = rgamma(1, 1, 1))
                
            }
            
        } else if (iddFun == 'logitIDD') {
            
            iddParamsPrior <- function(x) {
                dnorm(x['mid'], 8, 1, log = T) +
                    dgamma(x['rate'], 1, 1, log = T)
            }
            
            iddParamsInit = list(mid = rnorm(1, 8, 1),
                                 rate = runif(1, 0.5, 2))
            
        } else if (iddFun == 'splineIDD') {
            
            iddParamsPrior <- function(x) {
                sum(dnorm(x, 0, 4, log = T)) 
            }
            
            iddParamsInit = list(b1 = rnorm(1, 0, 2),
                                 b2 = rnorm(1, 0, 2),
                                 b3 = rnorm(1, 0, 2),
                                 b4 = rnorm(1, 0, 2),
                                 b5 = rnorm(1, 0, 2),
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
