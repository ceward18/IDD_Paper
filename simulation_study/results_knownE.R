################################################################################
# Combine batch output for known exposure period
################################################################################

library(BayesSEIR)
library(ggplot2)

myTheme <-  theme_bw() +
    theme(plot.title = element_text(h = 0.5, size = 16),
          legend.title = element_text(size = 14, h = 0.5),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.background=element_rect(fill="white"),
          strip.text = element_text(size = 14))


batchFiles <- list.files('./batch_output/')
batchFilesKnown <- batchFiles[grep('knownE', batchFiles)]
batchFilesEst <- batchFiles[grep('estimatedE', batchFiles)]



################################################################################

# combine results with known infectious period
iddCurveKnown <- readRDS(paste0('./batch_output/', batchFilesKnown[1]))
iddCurveKnown$infPeriodSpec <- 'IDD'

for (i in 2:length(batchFilesKnown)) {
    knownE_i <- readRDS(paste0('./batch_output/', batchFilesKnown[i]))
    
    if (i %in% c(2, 9:16)) {
        knownE_i$infPeriodSpec <- 'IDD'
    }
    
    iddCurveKnown <-rbind.data.frame(iddCurveKnown, knownE_i)
}

iddCurveKnown$EstarType <- 'Known'

# combine results with est infectious period
estE_i <- readRDS(paste0('./batch_output/', batchFilesEst[1]))

gdiagEst <- estE_i$gdiagAll
postParamsEst <- estE_i$postParamsAll
iddCurveEst <- estE_i$iddCurveAll
r0Est <- estE_i$r0All
mcmcEffEst <- estE_i$mcmcEffAll
waicEst <- estE_i$waicAll

for (i in 2:length(batchFilesEst)) {
    estE_i <- readRDS(paste0('./batch_output/', batchFilesEst[i]))
    
    gdiagEst <-rbind.data.frame(gdiagEst, estE_i$gdiagAll)
    postParamsEst <-rbind.data.frame(postParamsEst, estE_i$postParamsAll)
    iddCurveEst <-rbind.data.frame(iddCurveEst, estE_i$iddCurveAll)
    r0Est <-rbind.data.frame(r0Est, estE_i$r0All)
    mcmcEffEst <-rbind.data.frame(mcmcEffEst, estE_i$mcmcEffAll)
    waicEst <-rbind.data.frame(waicEst, estE_i$waicAll)
}

gdiagEst$EstarType <- 'Estimated'
postParamsEst$EstarType <- 'Estimated'
iddCurveEst$EstarType <- 'Estimated'
r0Est$EstarType <- 'Estimated'
mcmcEffEst$EstarType <- 'Estimated'
waicEst$EstarType <- 'Estimated'

# combine IDD Curves
iddCurveAll <- rbind.data.frame(
    iddCurveKnown[-which(colnames(iddCurveKnown) == 'allConverge')], 
    iddCurveEst)

################################################################################
# Check convergence

# known infectious period
notConvergeKnown <- iddCurveKnown[!iddCurveKnown$allConverge,]
notConvergeKnown <- notConvergeKnown[notConvergeKnown$infDay == 1,]

table(notConvergeKnown$iddFun, notConvergeKnown$datGen, notConvergeKnown$maxInf)

# estimated infectious period
notConvergeEst <- gdiagEst[gdiagEst$gr > 1.1,]
notConvergeModels <-  notConvergeEst[
    !duplicated(notConvergeEst
                [,-which(colnames(notConvergeEst) %in% 
                             c('gr', 'grUpper', 'param'))]),
    c('datGen', 'infPeriodSpec', 'iddFun', 'simNumber', 'maxInf')]


with(subset(notConvergeModels, maxInf == 15 & infPeriodSpec == 'IDD'),
     table(datGen, iddFun))

with(subset(notConvergeModels, maxInf == 20 & infPeriodSpec == 'IDD'),
     table(datGen, iddFun))

with(subset(notConvergeModels, infPeriodSpec == 'Exp'),
     table(datGen, maxInf))

with(subset(notConvergeModels, infPeriodSpec == 'PS'),
     table(datGen, maxInf))

################################################################################
# plot median IDD curve compared to truth

maxInf <- 15
N <- 5363500

### get true curves
iddPeakCurve <- dgammaIDD(1:maxInf, params = list(shape = 4, rate = 1))
iddPeakCurve <-  1 - exp(-exp(0.25) * iddPeakCurve /N)

iddExpCurve <- dgammaIDD(1:maxInf, params = list(shape = 0.9, rate = 0.3))
iddExpCurve <-  1 - exp(-exp(0.4) * iddExpCurve /N)

iddLogitCurve <- logitIDD(1:maxInf, params = list(mid = 8, rate = 1.5))
iddLogitCurve <-  1 - exp(-exp(-1.77) * iddLogitCurve /N)


trueCurves <- data.frame(datGen = rep(c('IDD_peak', 'IDD_exp', 'IDD_logit'), 
                                      each = maxInf),
                         infDay = 1:maxInf,
                         iddCurveTruth = c(iddPeakCurve, iddExpCurve, iddLogitCurve))


# merge
knownEAll <- merge(knownEAll, trueCurves, by = c('datGen', 'infDay'), all.x = T)
knownEAll <- knownEAll[order(knownEAll$datGen, knownEAll$iddFun,
                             knownEAll$simNumber, knownEAll$maxInf, knownEAll$infDay),]

# format for plotting
knownEAll$iddFun <- factor(knownEAll$iddFun,
                           levels = c('dgammaIDD', 'dlnormIDD', 'logitIDD', 'splineIDD'),
                           labels = c('Gamma pdf', 'Log-normal pdf', 'Logistic Decay', 'Basis Spline'))

knownEAll$datGen <- factor(knownEAll$datGen,
                           levels = c('IDD_peak', 'IDD_exp', 'IDD_logit'),
                           labels = c('IDD Peak', 'IDD Exp', 'IDD Logit'))

knownEAll$EstarType <- 'known'

pal <- c('magenta3', 'darkorange', 'green3', 'royalblue2')

# results using true function

ggplot(subset(knownEAll, iddFun != 'Basis Spline' & maxInf == 15)) +
    geom_line(aes(x = infDay, y = iddCurveMedian, group = simNumber),
              color =  adjustcolor('grey', alpha = 0.5)) +
    facet_wrap(~datGen + iddFun, nrow = 1) + 
    geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1) +
    myTheme +
    labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
         title = 'Known Exposure Times')  +
    scale_color_manual(values =pal) +
    theme(legend.position = "none")

ggplot(subset(knownEAll, iddFun != 'Basis Spline' & maxInf == 20)) +
    geom_line(aes(x = infDay, y = iddCurveMedian, group = simNumber),
              color =  adjustcolor('grey', alpha = 0.5)) +
    facet_wrap(~datGen + iddFun, nrow = 1) + 
    geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1) +
    myTheme +
    labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
         title = 'Known Exposure Times')  +
    scale_color_manual(values =pal) +
    theme(legend.position = "none")


## using basis splines


ggplot(subset(knownEAll, iddFun == 'Basis Spline' & maxInf == 15)) +
    geom_line(aes(x = infDay, y = iddCurveMedian, group = simNumber),
              color =  adjustcolor('grey', alpha = 0.5)) +
    facet_wrap(~datGen + iddFun, nrow = 1) + 
    geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1) +
    myTheme +
    labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
         title = 'Known Exposure Times')  +
    scale_color_manual(values =pal) +
    theme(legend.position = "none")


ggplot(subset(knownEAll, iddFun == 'Basis Spline' & maxInf == 20)) +
    geom_line(aes(x = infDay, y = iddCurveMedian, group = simNumber),
              color =  adjustcolor('grey', alpha = 0.5)) +
    facet_wrap(~datGen + iddFun, nrow = 1) + 
    geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1) +
    myTheme +
    labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
         title = 'Known Exposure Times')  +
    scale_color_manual(values =pal) +
    theme(legend.position = "none")





############################################################################
# MSE of R0
# 1/100 sum (mean_i - truth)^2

fullX <- getX(250, 120)

# get true R0's
iddPeakR0 <- getR0(infPeriodSpec = 'IDD', 
                   beta = c(0.25, -7), X = fullX, N = N,
                   infIDDParams = list(iddFun = dgammaIDD,
                                       iddParams = list(shape = 4, rate = 1),
                                       maxInf = 15))

iddExpR0 <- getR0(infPeriodSpec = 'IDD', 
                  beta = c(0.4, -7), X = fullX, N = N,
                  infIDDParams = list(iddFun = dgammaIDD,
                                      iddParams = list(shape = 0.9, rate = 0.3),
                                      maxInf = 15))

iddLogitR0 <- getR0(infPeriodSpec = 'IDD', 
                    beta = c(-1.77, -7), X = fullX, N = N,
                    infIDDParams = list(iddFun = logitIDD,
                                        iddParams = list(mid = 8, rate = 1.5),
                                        maxInf = 15))

iddPSR0 <- getR0(infPeriodSpec = 'PS', 
                 beta = c(-1.77, -7), X = fullX, N = N,
                 infPSParams = list(dist = 'gamma',
                                    psParams = list(shape = 56, rate = 7),
                                    maxInf = 15))

trueR0 <- data.frame(datGen = rep(c('IDD_peak', 'IDD_exp', 'IDD_logit', 'PS'), 
                                  each = 15),
                     infDay = 1:15,
                     truth = c(iddPeakR0, iddExpR0, iddLogitR0, iddPSR0))




