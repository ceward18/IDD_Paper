################################################################################
# Function to compute summary metrics from model ouptut
#
# Called from model_fitting_estimatedE.R
#
# For each model/simulation, want:
#   Gelman Rubin to ensure convergence
#   Posterior mean and 95% CI's for each parameter
#   Posterior median estimated IDD curve
#   Posterior mean estimated R0(t)
#   WAIC
#   MCMC Efficiency
################################################################################

post_processing <- function(modelOutput, EType, infPeriodSpec, 
                            datGen, iddFun, simNumber, maxInf,
                            X, N) {
    
    # thin to use only every 100th iteration
    nthin <- 100
    idxKeep <- seq(1, nrow(modelOutput[[1]]$fullPost), nthin)
    
    samples1 <- modelOutput[[1]]$fullPost[idxKeep,]
    samples2 <- modelOutput[[2]]$fullPost[idxKeep,]
    samples3 <- modelOutput[[3]]$fullPost[idxKeep,]
    
    paramsPost <- rbind.data.frame(samples1,
                                   samples2, 
                                   samples3)
    
    paramsPost <- cbind.data.frame(data.frame(datGen = datGen,
                                              infPeriodSpec = infPeriodSpec,
                                              iddFun = iddFun,
                                              simNumber = simNumber,
                                              maxInf = maxInf),
                                   paramsPost)
    
    ############################################################################
    ### gelman-rubin
    res_mcmc <- mcmc.list(mcmc(samples1), 
                          mcmc(samples2), 
                          mcmc(samples3))
    gdiag <- data.frame(gelman.diag(res_mcmc, multivariate = F)$psrf)
    colnames(gdiag) <- c('gr', 'grUpper')
    gdiag$param <- rownames(gdiag)
    rownames(gdiag) <- NULL
    gdiag <- cbind.data.frame(data.frame(datGen = datGen,
                                         infPeriodSpec = infPeriodSpec,
                                         iddFun = iddFun,
                                         simNumber = simNumber,
                                         maxInf = maxInf),
                              gdiag)
    
    ############################################################################
    ### posterior mean and 95% CI for parameters
    
    postMeans <- colMeans(paramsPost)
    postCI <- apply(paramsPost, 2, quantile, probs = c(0.025, 0.975))
    paramsSummary <- data.frame(param = names(postMeans),
                                mean = postMeans,
                                lower = postCI[1,],
                                upper = postCI[2,])
    rownames(paramsSummary) <- NULL
    
    
    ############################################################################
    ### posterior median IDD curve
    # (only calculated if model fit was IDD)
    
    if (infPeriodSpec == 'IDD') {
        
        iddParamsPost <- paramsPost[,-grep('beta1|beta2|rateE', colnames(paramsPost))]
        beta1Post <- paramsPost[,grep('beta1', colnames(paramsPost))]
        
        iddCurveFun <- get(iddFun)
        
        p0SE <- matrix(NA, nrow = maxInf, ncol = nrow(paramsPost))
        for (j in 1:nrow(paramsPost)) {
            
            if (iddFun == 'splineIDD') {
                
                XBasis <- bs(1:maxInf, df = 5)
                iddCurve <- iddCurveFun(1:maxInf, 
                                        params = as.list(iddParamsPost[j,]),
                                        XBasis = XBasis)
                
            } else {
                iddCurve <- iddCurveFun(1:maxInf, 
                                        params = as.list(iddParamsPost[j,]))
            }
            
            p0SE[,j] <- 1 - exp(-exp(beta1Post[j]) * iddCurve /N)
            
        }
        
        curveMedian <- apply(p0SE, 1, median)
        curveCI <- apply(p0SE, 1, quantile, probs = c(0.025, 0.975))
        
        iddSummary <- data.frame(datGen = datGen,
                                 infPeriodSpec = infPeriodSpec,
                                 iddFun = iddFun,
                                 simNumber = simNumber,
                                 maxInf = maxInf,
                                 infDay = 1:maxInf,
                                 median = curveMedian,
                                 lower = curveCI[1,],
                                 upper = curveCI[1,])
    } else {
        iddSummary <- data.frame(datGen = NA,
                                 infPeriodSpec = NA,
                                 iddFun = NA,
                                 simNumber = NA,
                                 maxInf = NA,
                                 infDay = NA,
                                 median = NA,
                                 lower = NA,
                                 upper = NA)
    }
    
    ############################################################################
    ### posterior mean R0(t)
    
    r0time <- matrix(NA, nrow = nrow(X), ncol = nrow(paramsPost))
    for (j in 1:nrow(paramsPost)) {
        
        betaSamp <- unlist(paramsPost[j, c('beta1', 'beta2')])
        if (infPeriodSpec == 'exp') {
            rateISamp <- paramsPost[j, 'rateI']
            r0time[,j] <- getR0(infPeriodSpec = infPeriodSpec, 
                                beta = betaSamp, X = X, N = N,
                                infExpParams = list(rateI = rateISamp))
            
        } else if (infPeriodSpec == 'PS') {
            
            psParamsSamp <- list(shape = paramsPost[j, 'shape'],
                                 rate = paramsPost[j, 'rate'])
            
            r0time[,j] <- getR0(infPeriodSpec = infPeriodSpec,
                                beta = betaSamp, X = X, N = N,
                                infPSParams = list(dist = 'gamma',
                                                   psParams = psParamsSamp,
                                                   maxInf = maxInf))
            
        } else if (infPeriodSpec == 'IDD') {
            
            iddCurveFun <- get(iddFun)
            iddParamsSamp <- as.list(iddParamsPost[j,])
            
            if (iddFun == 'splineIDD') {
                iddParamsSamp$XBasis <- XBasis
            }
            
            r0time[,j] <- getR0(infPeriodSpec = infPeriodSpec, 
                                beta = betaSamp, X = X, N = N,
                                infIDDParams = list(iddFun = iddCurveFun,
                                                    iddParams = iddParamsSamp,
                                                    maxInf = maxInf))
        }
    }
    
    r0Mean <- apply(r0time, 1, mean)
    r0CI <- apply(r0time, 1, quantile, probs = c(0.025, 0.975))
    
    r0Summary <- data.frame(datGen = datGen,
                            infPeriodSpec = infPeriodSpec,
                            iddFun = iddFun,
                            simNumber = simNumber,
                            maxInf = maxInf,
                            time = 1:length(r0Mean),
                            mean = r0Mean,
                            lower = r0CI[1,],
                            upper = r0CI[1,])
    
    
    ############################################################################
    # MCMC Efficiency = ESS / computation time
    if (EType == 'estimated') {
        
        
        # average time in minutes
        avgTime <- mean(c(modelOutput[[1]]$chainTime,
                          modelOutput[[2]]$chainTime,
                          modelOutput[[3]]$chainTime))
        
        
        # one for each parameter
        mcmcEff <- effectiveSize(res_mcmc) /  avgTime
        
        # R0 efficiency
        chainLength <- length(idxKeep)
        r0_mcmc_list <- mcmc.list(mcmc(r0time[1,1:chainLength]),
                                  mcmc(r0time[1,(chainLength + 1):(2 * chainLength)]),
                                  mcmc(r0time[1,(2 * chainLength + 1):(3 * chainLength)]))
        
        mcmcEffR0 <- effectiveSize(r0_mcmc_list) /  avgTime
        
        mcmcEffSummary <- data.frame(datGen = datGen,
                                     infPeriodSpec = infPeriodSpec,
                                     iddFun = iddFun,
                                     simNumber = simNumber,
                                     maxInf = maxInf,
                                     param = c(names(mcmcEff), 'R0'),
                                     eff = c(mcmcEff, mcmcEffR0))
        
        
    } else {
        mcmcEffSummary <- data.frame(datGen = NA,
                                     infPeriodSpec = NA,
                                     iddFun = NA,
                                     simNumber = NA,
                                     maxInf = NA,
                                     param = NA,
                                     eff = NA)
        
    }
    
    
    ############################################################################
    # WAIC
    
    if (EType == 'estimated') {
        
        avgWAIC <- mean(c(modelOutput[[1]]$WAIC,
                          modelOutput[[2]]$WAIC,
                          modelOutput[[3]]$WAIC))
        
        
        waicSummary <- data.frame(datGen = datGen,
                                  infPeriodSpec = infPeriodSpec,
                                  iddFun = iddFun,
                                  simNumber = simNumber,
                                  maxInf = maxInf,
                                  waic = avgWAIC)
        
        
    } else {
        waicSummary <- data.frame(datGen = NA,
                                  infPeriodSpec = NA,
                                  iddFun = NA,
                                  simNumber = NA,
                                  maxInf = NA,
                                  waic = NA)
        
    }
    
    
    list(gdiag = gdiag,
         paramsSummary = paramsSummary,
         iddSummary = iddSummary,
         r0Summary = r0Summary,
         mcmcEffSummary = mcmcEffSummary,
         waicSummary = waicSummary)
    
}