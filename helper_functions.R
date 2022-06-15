################################################################################
# Helper functions
################################################################################

################################################################################
# function to construct design matrix for intervention 

getX <- function(ntpts, tint){
    # Intercept and intervention term
    cbind(1, cumsum(1:ntpts > tint) / 100)
}


################################################################################
# function to run models in simulation until convergence is acheived 

runModels <- function(datList, X, infPeriodSpec_i, iddFun_i, datGen_i, maxInf_i,
                      EType, i, niter) {
    
    assign('infPeriodSpec_i', infPeriodSpec_i, envir = .GlobalEnv)
    assign('EType', EType, envir = .GlobalEnv)
    assign('niter', niter, envir = .GlobalEnv)
    
    print(paste0('Attempting niter = ', niter))
    
    cl <- makeCluster(3, outfile="")
    clusterExport(cl, list('datList',  'X', 'infPeriodSpec_i', 'iddFun_i',
                           'datGen_i', 'maxInf_i', 'EType',  'i', 'niter'))
    
    resThree <- parLapplyLB(cl, X = 1:3, fun = function(clIdx) {
        
        library(BayesSEIR)
        library(splines)
        
        # MCMC specifications
        nburn <- 2e5
        
        # set seed for reproducibility
        set.seed(clIdx + i)
        
        # get priors and initial values based on model/data generating scenario
        source('get_priors_inits.R')
        priorsInits <- get_priors_inits(infPeriodSpec = infPeriodSpec_i, 
                                        iddFun = iddFun_i, 
                                        datGen = datGen_i, 
                                        maxInf = maxInf_i) 
        
        initsList<- priorsInits$initsList 
        priorList<- priorsInits$priorList 
        
        # WAIC only needs to be calculated for estimated exposure times
        if (EType == 'known') {
            waic <- FALSE
        } else if (EType == 'estimated') {
            waic <- TRUE
        }
        
        # start timing for MCMC efficiency
        startTime <- Sys.time()
        
        if (infPeriodSpec_i == 'exp') {
            
            res <-  mcmcSEIR(dat = datList, X = X, 
                             inits = initsList, 
                             niter = niter, nburn = nburn,
                             infPeriodSpec = infPeriodSpec_i,
                             priors = priorList,
                             WAIC = waic, seed = clIdx)
            
        } else if (infPeriodSpec_i == 'PS') {
            
            res <- mcmcSEIR(dat = datList, X = X, 
                            inits = initsList, 
                            niter = niter, nburn = nburn,
                            infPeriodSpec = infPeriodSpec_i,
                            priors = priorList,
                            dist = 'gamma', maxInf = maxInf_i,
                            WAIC = waic, seed = clIdx)
            
        } else if (infPeriodSpec_i == 'IDD') {
            
            res <-  mcmcSEIR(dat = datList, X = X, 
                             inits = initsList, 
                             niter = niter, nburn = nburn,
                             infPeriodSpec = infPeriodSpec_i,
                             priors = priorList,
                             iddFun = iddFun_i, maxInf = maxInf_i,
                             WAIC = waic, seed = clIdx)
        }
        
        endTime <- Sys.time()
        
        res$chainTime <- as.numeric(endTime - startTime, units = 'mins')
        res
        
    })
    stopCluster(cl)
    
    ### post processing to get Gelman-Rubin and other summaries of interest
    postSummaries <- post_processing(modelOutput = resThree, EType = EType_i,
                                     infPeriodSpec = infPeriodSpec_i, 
                                     datGen = datGen_i, iddFun = iddFun_i, 
                                     simNumber = simNumber_i, maxInf = maxInf_i,
                                     X = X, N = N)
    
    
    allConverge <- all(postSummaries$gdiag$gr < 1.1)
    
    if (allConverge) {
        return(postSummaries)
    } else {
        runModels(datList, X, iddFun_i, datGen_i, maxInf_i, i, niter = niter + 2e5)
    }
    
    
}
