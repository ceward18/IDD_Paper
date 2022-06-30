################################################################################
# Combine batch output for known exposure period
################################################################################

library(BayesSEIR)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)
library(RColorBrewer)
library(scales)

myTheme <-  theme_bw() +
    theme(strip.background=element_rect(fill="white"))

source('../helper_functions.R')

################################################################################
### Supplemental Table 1: Simulation characteristics (median size and length)

# load all epidemics
psDat <- readRDS('./Data/PS_data.rds')
iddPeakDat <- readRDS('./Data/IDD_peak_data.rds')
iddExpDat <- readRDS('./Data/IDD_exp_data.rds')
iddLogitDat <- readRDS('./Data/IDD_logit_data.rds')

# total number infected from each simulation
medianSize <- c(median(colSums(psDat$Istar)),
                median(colSums(iddPeakDat$Istar)),
                median(colSums(iddExpDat$Istar)),
                median(colSums(iddLogitDat$Istar)))

# last day of new case + maximum duration of infectious period
medianLength <- c(median(apply(psDat$Istar, 2, function(x) max(which(x > 1))) + 15),
                  median(apply(iddPeakDat$Istar, 2, function(x) max(which(x > 1))) + 15),
                  median(apply(iddExpDat$Istar, 2, function(x) max(which(x > 1))) + 15),
                  median(apply(iddLogitDat$Istar, 2, function(x) max(which(x > 1))) + 15))

data.frame(Scenarios = c('PS', 'IDD Peak', 'IDD Exp', 'IDD Logit'),
           medianSize = medianSize,
           medianLength = medianLength)

################################################################################
### Combine batch files
# Note: takes about 10 min to combine all batches

batchFiles <- list.files('./batch_output/')
batchFilesKnown <- batchFiles[grep('knownE', batchFiles)]
batchFilesEst <- batchFiles[grep('estimatedE', batchFiles)]

# allBatchFiles <- paste0("estimatedE_batch", 1:4400, ".rds")
# allBatchFiles[!allBatchFiles %in% batchFilesEst]

# combine results with known infectious period
knownE_i <- readRDS(paste0('./batch_output/', batchFilesKnown[1]))

gdiagKnown <- knownE_i$gdiag
postParamsKnown <- knownE_i$paramsSummary
iddCurveKnown <- knownE_i$iddSummary
r0Known <- knownE_i$r0Summary
mcmcEffKnown <- knownE_i$mcmcEffSummary

for (i in 2:length(batchFilesKnown)) {
    knownE_i <- readRDS(paste0('./batch_output/', batchFilesKnown[i]))
    
    gdiagKnown <-rbind.data.frame(gdiagKnown, knownE_i$gdiag)
    postParamsKnown <-rbind.data.frame(postParamsKnown, knownE_i$paramsSummary)
    iddCurveKnown <-rbind.data.frame(iddCurveKnown, knownE_i$iddSummary)
    r0Known <-rbind.data.frame(r0Known, knownE_i$r0Summary)
    mcmcEffKnown <-rbind.data.frame(mcmcEffKnown, knownE_i$mcmcEffSummary)
}


# combine results with estimated infectious period
estE_i <- readRDS(paste0('./batch_output/', batchFilesEst[1]))

gdiagEst <- estE_i$gdiag
postParamsEst <- estE_i$paramsSummary
iddCurveEst <- estE_i$iddSummary
r0Est <- estE_i$r0Summary
mcmcEffEst <- estE_i$mcmcEffSummary
waicEst <- estE_i$waicSummary

for (i in 2:length(batchFilesEst)) {
    estE_i <- readRDS(paste0('./batch_output/', batchFilesEst[i]))
    
    gdiagEst <-rbind.data.frame(gdiagEst, estE_i$gdiag)
    postParamsEst <-rbind.data.frame(postParamsEst, estE_i$paramsSummary)
    iddCurveEst <-rbind.data.frame(iddCurveEst, estE_i$iddSummary)
    r0Est <-rbind.data.frame(r0Est, estE_i$r0Summary)
    mcmcEffEst <-rbind.data.frame(mcmcEffEst, estE_i$mcmcEffSummary)
    waicEst <-rbind.data.frame(waicEst, estE_i$waicSummary)
}


# combine IDD Curves
iddCurveAll <- rbind.data.frame(iddCurveKnown, iddCurveEst)

################################################################################
# Double check convergence 

# known infectious period (should be 0)
sum(gdiagKnown[gdiagKnown$gr > 1.1,])

# estimated infectious period (should be 0)
sum(gdiagEst[gdiagEst$gr > 1.1,])

################################################################################
# Check all models have finished

# should be 16
length(batchFilesKnown)

# should be 4400
length(batchFilesEst)


############################################################################
### Figure 1: IDD transmissibility curves used in the simulation

maxInf <- 15
N <- 5363500

### get true curves
iddPeakCurve <- dgammaIDD(1:maxInf, params = list(shape = 4, rate = 1))
iddPeakCurve <-  1 - exp(-exp(0.25) * iddPeakCurve /N)

iddExpCurve <- dgammaIDD(1:maxInf, params = list(shape = 0.9, rate = 0.3))
iddExpCurve <-  1 - exp(-exp(0.4) * iddExpCurve /N)

iddLogitCurve <- logitIDD(1:maxInf, params = list(mid = 8, rate = 1.5))
iddLogitCurve <-  1 - exp(-exp(-1.77) * iddLogitCurve /N)


pal <- c('magenta3', 'darkorange', 'green3')

# bottom, left, top, and right.
# pdf('../figures/main/sim_iddCurves_fig1.pdf', height = 3.5, width = 9.5)
par(mfrow = c(1, 3), mar = c(5.1, 2.1, 4.1, 2.1))
plot(1:maxInf, iddPeakCurve, ylab = '', yaxt = 'n', type = 'l',
     xlab = 'Days Infectious', cex.axis = 1.5, cex.main = 2, cex.lab = 2,
     main = 'IDD Peak', lwd = 3, col = pal[1])

plot(1:maxInf, iddExpCurve, ylab = '', yaxt = 'n', type = 'l',
     xlab = 'Days Infectious', cex.axis = 1.5, cex.main = 2, cex.lab = 2,
     main = 'IDD Exp', lwd = 3, col = pal[2])

plot(1:maxInf, iddLogitCurve, ylab = '', yaxt = 'n', type = 'l',
     xlab = 'Days Infectious', cex.axis = 1.5, cex.main = 2, cex.lab = 2,
     main = 'IDD Logit', lwd = 3, col = pal[3])
# dev.off()


################################################################################
### Figure 2: Posterior median IDD transmissibility curves compared to truth

trueCurves <- data.frame(datGen = rep(c('IDD_peak', 'IDD_exp', 'IDD_logit'), 
                                      each = maxInf),
                         infDay = 1:maxInf,
                         iddCurveTruth = c(iddPeakCurve, iddExpCurve, iddLogitCurve))

iddCurveAllSub <- merge(subset(iddCurveAll, datGen != 'PS'), 
                        trueCurves, by = c('datGen', 'infDay'), all.x = T)

iddCurveAllSub$iddCurveTruth[iddCurveAllSub$infDay > 15] <- 0


iddCurveAllSub$iddFunLab <- factor(iddCurveAllSub$iddFun,
                                   levels = c('dgammaIDD', 'dlnormIDD',
                                              'logitIDD',  'splineIDD'),
                                   labels = paste0('Model: ',
                                                   c('Gamma pdf', 'Log-normal pdf', 
                                                     'Logistic Decay', 'Basis Spline')))

iddCurveAllSub$datGenLab <- factor(iddCurveAllSub$datGen,
                                   levels = c('IDD_peak', 'IDD_exp',
                                              'IDD_logit', 'PS'),
                                   labels = paste0('Data: ',
                                                   c('IDD Peak', 'IDD Exp',
                                                     'IDD Logit', 'PS')))


# when exposure times are estimated, we only care about certain scenarios for
# estimating IDD curves

iddCurveAllSub <- iddCurveAllSub[-which(iddCurveAllSub$datGen == 'IDD_peak' &
                                            iddCurveAllSub$iddFun == 'logitIDD'),]
iddCurveAllSub <- iddCurveAllSub[-which(iddCurveAllSub$datGen == 'IDD_exp' & 
                                            iddCurveAllSub$iddFun == 'logitIDD'),]
iddCurveAllSub <- iddCurveAllSub[-which(iddCurveAllSub$datGen == 'IDD_logit' & 
                                            iddCurveAllSub$iddFun %in% 
                                            c('dgammaIDD', 'dlnormIDD')),]

pal <- c('magenta3', 'darkorange', 'green3')

p1 <- ggplot(data = subset(iddCurveAllSub, 
                           iddFun != 'splineIDD' & 
                               EType == 'known' & 
                               maxInf == 20)) +
    geom_line(aes(x = infDay, y = median, group = simNumber),
              color =  adjustcolor('grey', alpha = 0.5)) +
    geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1.5) +
    facet_wrap(~datGenLab + iddFunLab, nrow = 1) + 
    myTheme +
    labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
         title = 'Known Exposure Times')  +
    scale_color_manual(values =pal) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(h = 0.5, size = 18),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 15))

p2 <- ggplot(data = subset(iddCurveAllSub, 
                           iddFun != 'splineIDD' & 
                               EType == 'estimated' & 
                               maxInf == 20)) +
    geom_line(aes(x = infDay, y = median, group = simNumber),
              color =  adjustcolor('grey', alpha = 0.5)) +
    geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1.5) +
    facet_wrap(~datGenLab + iddFunLab, nrow = 1) + 
    myTheme +
    labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
         title = 'Estimated Exposure Times')  +
    scale_color_manual(values =pal) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(h = 0.5, size = 18),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 15))

pdf('../figures/main/sim_iddEst_fig2.pdf', height = 7.5, width = 14)
grid.arrange(p1, p2,
             top = textGrob(expression('             Posterior median estimates of '~pi[0]^(SE)),
                            gp = gpar(fontsize = 20, font = 2)))
dev.off()


################################################################################
# Figure 3: Posterior IDD curves using basis splines compared to truth


p1 <- ggplot(data = subset(iddCurveAllSub, 
                           iddFun == 'splineIDD' & 
                               EType == 'known' & 
                               maxInf == 20)) +
    geom_line(aes(x = infDay, y = median, group = simNumber),
              color =  adjustcolor('grey', alpha = 0.5)) +
    geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1.5) +
    facet_wrap(~datGenLab+iddFunLab, nrow = 1) + 
    myTheme +
    labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
         title = 'Known Exposure Times')  +
    scale_color_manual(values =pal) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(h = 0.5, size = 18),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 15))

p2 <- ggplot(data = subset(iddCurveAllSub, 
                           iddFun == 'splineIDD' & 
                               EType == 'estimated' & 
                               maxInf == 20)) +
    geom_line(aes(x = infDay, y = median, group = simNumber),
              color =  adjustcolor('grey', alpha = 0.5)) +
    geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1.5) +
    facet_wrap(~datGenLab+iddFunLab, nrow = 1) + 
    myTheme +
    labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
         title = 'Estimated Exposure Times')  +
    scale_color_manual(values =pal) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(h = 0.5, size = 18),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 15))

pdf('../figures/main/sim_iddEstSpline_fig3.pdf', height = 7.5, width = 14)
grid.arrange(p1, p2,
             top = textGrob(expression('             Posterior median estimates of '~pi[0]^(SE)),
                            gp = gpar(fontsize = 20, font = 2)))
dev.off()



################################################################################
### Supplemental Figure 1: Posterior median IDD transmissibility curves ALL


iddCurveAllSub <- merge(iddCurveAll, 
                        trueCurves, by = c('datGen', 'infDay'), all.x = T)

iddCurveAllSub$iddCurveTruth[iddCurveAllSub$infDay > 15 & 
                                 iddCurveAllSub$datGen != 'PS'] <- 0
iddCurveAllSub <- subset(iddCurveAllSub, EType == 'estimated')
iddCurveAllSub <- subset(iddCurveAllSub, infPeriodSpec == 'IDD')

iddCurveAllSub$iddFun <- factor(iddCurveAllSub$iddFun,
                                levels = c('dgammaIDD', 'dlnormIDD',
                                           'logitIDD',  'splineIDD'),
                                labels = paste0('',
                                                c('Gamma pdf', 'Log-normal pdf', 
                                                  'Logistic Decay', 'Basis Spline')))

iddCurveAllSub$datGen <- factor(iddCurveAllSub$datGen,
                                levels = c('IDD_peak', 'IDD_exp',
                                           'IDD_logit', 'PS'),
                                labels = paste0('Data: ',
                                                c('IDD Peak', 'IDD Exp',
                                                  'IDD Logit', 'PS')))

pal <- c('magenta3', 'darkorange', 'green3', 'black')

maxInfLabs <- c('15-day inf. period',
                '20-day inf. period')
names(maxInfLabs) <- c('15', '20')
iddCurveAllSub$maxInf <- as.character(iddCurveAllSub$maxInf)

pdf('../figures/supplemental/sim_iddEst_suppfig1.pdf', height = 10, width = 20)
ggplot(data = subset(iddCurveAllSub)) +
    geom_line(aes(x = infDay, y = median, group = simNumber),
              color =  adjustcolor('grey', alpha = 0.5)) +
    geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1.8) +
    facet_grid(iddFun~datGen + maxInf, scales = 'free',
               labeller = labeller(maxInf = maxInfLabs)) +  
    myTheme +
    labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
         title = expression('Posterior median estimates of '~pi[0]^(SE)~' for all scenarios'))  +
    scale_color_manual(values =pal) +
    theme(legend.position = "none",
          strip.text = element_text(size = 18),
          plot.title = element_text(size = 22, h = 0.5),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
dev.off()



################################################################################
# Figure 4: MSE of R0

# get true R0
maxInf <- 15
N <- 5363500
fullX <- getX(250, 120)

# get true R0's
iddPSR0 <- getR0(infPeriodSpec = 'PS', 
                 beta = c(-1.77, -7), X = fullX, N = N,
                 infPSParams = list(dist = 'gamma',
                                    psParams = list(shape = 56, rate = 7),
                                    maxInf = maxInf))
iddPeakR0 <- getR0(infPeriodSpec = 'IDD', 
                   beta = c(0.25, -7), X = fullX, N = N,
                   infIDDParams = list(iddFun = dgammaIDD,
                                       iddParams = list(shape = 4, rate = 1),
                                       maxInf = maxInf))

iddExpR0 <- getR0(infPeriodSpec = 'IDD', 
                  beta = c(0.4, -7), X = fullX, N = N,
                  infIDDParams = list(iddFun = dgammaIDD,
                                      iddParams = list(shape = 0.9, rate = 0.3),
                                      maxInf = maxInf))

iddLogitR0 <- getR0(infPeriodSpec = 'IDD', 
                    beta = c(-1.77, -7), X = fullX, N = N,
                    infIDDParams = list(iddFun = logitIDD,
                                        iddParams = list(mid = 8, rate = 1.5),
                                        maxInf = maxInf))


trueR0 <- data.frame(datGen = rep(c('PS', 'IDD_peak', 'IDD_exp', 'IDD_logit'), 
                                  each = 250),
                     time = 1:250,
                     truth = c(iddPSR0, iddPeakR0, iddExpR0, iddLogitR0))

# merge estimates with truth
r0EstTrue <- merge(r0Est, trueR0, by = c('datGen', 'time'), all.x = T)

r0MSE <- ddply(subset(r0EstTrue, time == 1), .(datGen, infPeriodSpec, iddFun, maxInf),
               summarize,
               mse = mean((mean-truth)^2))

# barplot
r0MSE$fitType <- r0MSE$iddFun
r0MSE$fitType[r0MSE$infPeriodSpec == 'exp'] <- 'exp'
r0MSE$fitType[r0MSE$infPeriodSpec == 'PS'] <- 'PS'

# replicate exponential for maxInf of 15 as maxInf does not impact exponential
r0MSEExp <- r0MSE[r0MSE$fitType == 'exp',]
r0MSEExp$maxInf <- 15
r0MSE <- rbind.data.frame(r0MSE, r0MSEExp)

# format for plotting labels
r0MSE$fitType <- factor(r0MSE$fitType, 
                        levels = c('dlnormIDD', 'dgammaIDD',
                                   'logitIDD', 'splineIDD',
                                   'exp',  'PS'),
                        labels = c('IDD - Log-normal Pdf',
                                   'IDD - Gamma Pdf', 
                                   'IDD - Logistic Decay',
                                   'IDD - Basis Spline',
                                   'Exponential', 
                                   'Path-specific'))

r0MSE$datGen <- factor(r0MSE$datGen,
                       levels = c('IDD_peak', 'IDD_exp', 
                                  'IDD_logit', 'PS'),
                       labels = c('IDD Peak', 'IDD Exp', 
                                  'IDD Logit', 'PS'))

maxInfLabs <- c('15-day infectious period',
                '20-day infectious period')
names(maxInfLabs) <- c('15', '20')
r0MSE$maxInf <- as.character(r0MSE$maxInf)


pal <- brewer.pal(6, 'Dark2')
pal <- pal[c(1:3, 6, 4:5)]


pdf('../figures/main/sim_r0_mse_fig4.pdf', height = 4, width = 10)
ggplot(data = r0MSE, 
       aes(x = datGen, y = mse, fill = fitType, group = fitType)) +
    geom_bar(position="dodge", stat="identity", col = 'black', size = 0.2) +
    facet_wrap(~maxInf, nrow =1,
               labeller = labeller(maxInf = maxInfLabs)) +
    scale_fill_manual(values=pal) +
    myTheme +
    labs(title = expression('MSE of'~ R[0]),
         fill = 'Modeling Approach',
         y = 'MSE', x = 'Data Type') +
    theme(legend.text.align = 0,
          plot.title = element_text(h = 0.5, size = 16),
          legend.title = element_text(size = 12, h = 0.5),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 13))
dev.off()


################################################################################
# Supplemental Figure 2, 3: Posterior mean R0(t) compared to truth 

r0EstFull <- r0EstTrue

r0EstFull$fitType <- r0EstFull$iddFun
r0EstFull$fitType[r0EstFull$infPeriodSpec == 'exp'] <- 'exp'
r0EstFull$fitType[r0EstFull$infPeriodSpec == 'PS'] <- 'PS'

# replicate exponential for maxInf of 15 as maxInf does not impact exponential
r0EstFullExp <- r0EstFull[r0EstFull$fitType == 'exp',]
r0EstFullExp$maxInf <- 15
r0EstFull <- rbind.data.frame(r0EstFull, r0EstFullExp)

r0EstFull$datGen <- factor(r0EstFull$datGen, 
                           levels = c('IDD_peak', 'IDD_exp', 'IDD_logit', 'PS'),
                           labels = paste0('Data: ',
                                           c('IDD Peak', 'IDD Exp', 
                                             'IDD Logit', 'PS')))

r0EstFull$fitType <- factor(r0EstFull$fitType, 
                            levels = c('exp', 'PS', 'dgammaIDD', 
                                       'dlnormIDD', 'logitIDD', 'splineIDD'),
                            labels = paste0('Model: ', 
                                            c('Exponential','Path-specific', 
                                              'Gamma Pdf', 'Log-normal Pdf', 
                                              'Logistic Decay', 'Basis Spline')))


pal <- c('magenta3', 'darkorange', 'green3', 'royalblue2')


pdf('../figures/supplemental/sim_r0_est15_suppfig2.pdf', height = 6, width = 11.5)
ggplot(subset(r0EstFull, maxInf == 15), aes(x = time, group = simNumber)) + 
    geom_line(aes(y = mean), color =  adjustcolor('grey', alpha = 0.5)) +
    geom_line(aes(y = truth, col = datGen), size = 1) +
    facet_grid(datGen ~ fitType) +
    myTheme +
    labs(x = 'Epidemic Time', y = expression(R[0](t)),
         col = 'Data Generation\nScenario')  +
    scale_color_manual(values =pal) +
    ggtitle(expression('Posterior mean estimates of '~R[0](t)~' for 15-day infectious period')) +
    theme(legend.position = "none",
          plot.title = element_text(h = 0.5, size = 15),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 9),
          strip.text = element_text(size = 10.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
dev.off()

pdf('../figures/supplemental/sim_r0_est20_suppfig3.pdf', height = 6, width = 11.5)
ggplot(subset(r0EstFull, maxInf == 20), aes(x = time, group = simNumber)) + 
    geom_line(aes(y = mean), color =  adjustcolor('grey', alpha = 0.5)) +
    geom_line(aes(y = truth, col = datGen), size = 0.8) +
    facet_grid(datGen ~ fitType) +
    myTheme +
    labs(x = 'Epidemic Time', y = expression(R[0](t)),
         col = 'Data Generation\nScenario')  +
    scale_color_manual(values =pal) +
    ggtitle(expression('Posterior mean estimates of '~R[0](t)~' for 20-day infectious period')) +
    theme(legend.position = "none",
          plot.title = element_text(h = 0.5, size = 15),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 9),
          strip.text = element_text(size = 10.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
dev.off()


################################################################################
# Supplemental Figure 4: MSE of R0(t) (over time)

r0MSE_time <- ddply(r0EstTrue, .(datGen, infPeriodSpec, iddFun, maxInf, time),
                    summarize,
                    mse = mean((mean-truth)^2))

# barplot
r0MSE_time$fitType <- r0MSE_time$iddFun
r0MSE_time$fitType[r0MSE_time$infPeriodSpec == 'exp'] <- 'exp'
r0MSE_time$fitType[r0MSE_time$infPeriodSpec == 'PS'] <- 'PS'

# replicate exponential for maxInf of 15 as maxInf does not impact exponential
r0MSE_timeExp <- r0MSE_time[r0MSE_time$fitType == 'exp',]
r0MSE_timeExp$maxInf <- 15
r0MSE_time <- rbind.data.frame(r0MSE_time, r0MSE_timeExp)

# format for plotting labels
r0MSE_time$fitType <- factor(r0MSE_time$fitType, 
                             levels = c('dlnormIDD', 'dgammaIDD',
                                        'logitIDD', 'splineIDD',
                                        'exp',  'PS'),
                             labels = c('IDD - Log-normal Pdf',
                                        'IDD - Gamma Pdf', 
                                        'IDD - Logistic Decay',
                                        'IDD - Basis Spline',
                                        'Exponential', 
                                        'Path-specific'))

r0MSE_time$datGen <- factor(r0MSE_time$datGen,
                            levels = c('IDD_peak', 'IDD_exp', 
                                       'IDD_logit', 'PS'),
                            labels = paste0('Data: ', 
                                            c('IDD Peak', 'IDD Exp',
                                              'IDD Logit', 'PS')))


maxInfLabs <- c('15-day infectious period',
                '20-day infectious period')
names(maxInfLabs) <- c('15', '20')
r0MSE_time$maxInf <- as.character(r0MSE_time$maxInf)

pal <- brewer.pal(6, 'Dark2')
pal <- pal[c(1:3, 6, 4:5)]

pdf('../figures/supplemental/sim_r0t_mse_suppfig4.pdf', height = 6, width = 12.5)
ggplot(data = r0MSE_time, 
       aes(x = time, y = mse, col = fitType, linetype = fitType)) +
    geom_line(size = 0.8, alpha = 0.8) +
    facet_wrap(maxInf~datGen, nrow =2,
               labeller = labeller(maxInf = maxInfLabs)) +
    scale_color_manual(values=pal) +
    scale_linetype_manual(values = 1:6) +
    myTheme +
    labs(title = expression('MSE of'~ R[0](t)),
         col = 'Modeling Approach         ', 
         linetype = 'Modeling Approach         ', 
         y = 'MSE', x = 'Epidemic Time') +
    theme(legend.text.align = 0,
          legend.key.width = unit(1.5,"cm"),
          plot.title = element_text(h = 0.5, size = 16),
          legend.title = element_text(size = 14, h = 0.5),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 14))
dev.off()


################################################################################
# Table 1: Mean, SD of WAIC

waicEst$fitType <- waicEst$iddFun
waicEst$fitType[waicEst$infPeriodSpec == 'exp'] <- 'exp'
waicEst$fitType[waicEst$infPeriodSpec == 'PS'] <- 'PS'

# replicate exponential for maxInf of 15 as maxInf does not impact exponential
waicEstExp <- waicEst[waicEst$fitType == 'exp',]
waicEstExp$maxInf <- 15
waicEst <- rbind.data.frame(waicEst, waicEstExp)

ggplot(waicEst, aes(x = datGen, y = waic, fill = fitType)) +
    facet_wrap(~maxInf) +
    geom_boxplot(outlier.size = 0.4) +
    myTheme

waicSummary <- ddply(waicEst, .(maxInf, datGen, fitType), summarize,
                     waic = paste0(format(round(mean(waic), 2), nsmall = 2) ,
                                   ' (', format(round(sd(waic), 1), nsmall = 1),
                                   ')'))


waicSummary$fitType <- factor(waicSummary$fitType, 
                              levels = c('exp', 'PS', 'dgammaIDD', 'dlnormIDD',
                                         'logitIDD', 'splineIDD'),
                              labels = c('Exponential', 'Path-specific', 
                                         'IDD - Gamma pdf', 'IDD - Log-normal pdf',
                                         'IDD - Logistic Decay', 'IDD - Basis Spline'))

waicSummary$datGen <- factor(waicSummary$datGen,
                             levels = c('IDD_peak', 'IDD_exp', 
                                        'IDD_logit', 'PS'),
                             labels = c('IDD Peak', 'IDD Exp',
                                        'IDD Logit', 'PS'))


waicSummary <- waicSummary[order(waicSummary$maxInf, 
                                 waicSummary$datGen,
                                 waicSummary$fitType),]

w15 <- reshape(subset(waicSummary, maxInf == 15), 
               timevar = "datGen",
               idvar = c('fitType', 'maxInf'),
               direction = "wide")
w15 <- w15[,-1]
colnames(w15) <- c('fitType', paste0(levels(waicSummary$datGen), ' 15'))

w20 <- reshape(subset(waicSummary, maxInf == 20), 
               timevar = "datGen",
               idvar = c('fitType', 'maxInf'),
               direction = "wide")
w20 <- w20[,-1]
colnames(w20) <- c('fitType', paste0(levels(waicSummary$datGen), ' 20'))


tab1 <- merge(w15, w20, by = 'fitType')
tab1 <- tab1[,c(1, 2, 6, 3, 7, 4, 8, 5, 9)]
tab1 <- tab1[order(tab1$fitType),]
tab1
# 
# 
# selectedModel <- ddply(waicEst, .(maxInf, datGen, simNumber), summarize,
#                        min = min(waic),
#                        bestFit = fitType[which.min(waic)])
# 
# selectedModel$bestFit <- factor(selectedModel$bestFit, 
#                                 levels = c('exp', 'PS', 'dgammaIDD', 'dlnormIDD',
#                                            'logitIDD', 'splineIDD'),
#                                 labels = c('Exponential', 'Path-specific', 
#                                            'IDD - Gamma pdf', 'IDD - Log-normal pdf',
#                                            'IDD - Logistic Decay', 'IDD - Basis Spline'))
# 
# selectedModel$datGen <- factor(selectedModel$datGen,
#                                levels = c('IDD_peak', 'IDD_exp', 
#                                           'IDD_logit', 'PS'),
#                                labels = c('IDD Peak', 'IDD Exp',
#                                           'IDD Logit', 'PS'))
# 
# with(subset(selectedModel, maxInf == 15), table(datGen, bestFit))
# with(subset(selectedModel, maxInf == 20), table(datGen, bestFit))
# 
# waicTab15 <- with(subset(selectedModel, maxInf == 15), table(datGen, bestFit))
# waicTab20 <- with(subset(selectedModel, maxInf == 20), table(datGen, bestFit))
# 
# round(prop.table(waicTab15, 1) * 100, 2)
# round(prop.table(waicTab20, 1) * 100, 2)



################################################################################
# Table 2: Mean, minimum, R0 MCMC efficiency
# Supplemental table 3

mcmcEffAll <- mcmcEffEst

mcmcEffAll$fitType <- mcmcEffAll$iddFun
mcmcEffAll$fitType[mcmcEffAll$infPeriodSpec == 'exp'] <- 'exp'
mcmcEffAll$fitType[mcmcEffAll$infPeriodSpec == 'PS'] <- 'PS'


# replicate exponential for maxInf of 15 as maxInf does not impact exponential
mcmcEffAllExp <- mcmcEffAll[mcmcEffAll$fitType == 'exp',]
mcmcEffAllExp$maxInf <- 15
mcmcEffAll <- rbind.data.frame(mcmcEffAll, mcmcEffAllExp)


mcmcEffAll$fitType <- factor(mcmcEffAll$fitType, 
                             levels = c('exp', 'PS', 'dgammaIDD', 'dlnormIDD',
                                        'logitIDD', 'splineIDD'),
                             labels = c('Exponential', 'Path-specific', 
                                        'Gamma pdf', 'Log-normal pdf',
                                        'Logistic Decay', 'Basis Spline'))

# get mean and standard deviation for each scenario 
# collapse over data generating scenario and length of infectious period
mcmcEffSummary <- ddply(mcmcEffAll, .(fitType, param), summarize,
                        meanEff = mean(eff),
                        sdEff = sd(eff))

# rename parameters to match table 
mcmcEffSummary$param[mcmcEffSummary$param == 'rateE'] <- 'rho'
mcmcEffSummary$param[mcmcEffSummary$param == 'rateI'] <- 'gamma'

mcmcEffSummary$param[mcmcEffSummary$param == 'shape' & 
                         mcmcEffSummary$fitType == 'Path-specific'] <- 'alpha_I'
mcmcEffSummary$param[mcmcEffSummary$param == 'rate' & 
                         mcmcEffSummary$fitType == 'Path-specific'] <- 'beta_I'

mcmcEffSummary$param[mcmcEffSummary$param == 'shape' & 
                         mcmcEffSummary$fitType == 'Gamma pdf'] <- 'alpha'
mcmcEffSummary$param[mcmcEffSummary$param == 'rate' & 
                         mcmcEffSummary$fitType == 'Gamma pdf'] <- 'beta'

mcmcEffSummary$param[mcmcEffSummary$param == 'mid'] <- 'w_0'
mcmcEffSummary$param[mcmcEffSummary$param == 'rate' & 
                         mcmcEffSummary$fitType == 'Logistic Decay'] <- 'k'

minEff <- ddply(mcmcEffSummary, .( fitType), summarize,
                param = 'Min',
                meanEff = meanEff[which.min(meanEff)],
                sdEff = sdEff[which.min(meanEff)])
meanEff <- ddply(mcmcEffAll, .(fitType), summarize,
                 param = 'Mean',
                 meanEff = mean(eff),
                 sdEff = sd(eff))

mcmcEffSummary <- rbind.data.frame(mcmcEffSummary,
                                   meanEff,
                                   minEff)


mcmcEffSummary$mcmcEff <- paste0(format(round(mcmcEffSummary$meanEff, 1), nsmall = 1),
                                 ' (', trimws(format(round(mcmcEffSummary$sdEff, 1), nsmall = 1)), ')')

mcmcEffSummary$param <- factor(mcmcEffSummary$param,
                               levels = c('beta1', 'beta2', 'rho', 'gamma', 
                                          'alpha_I',  'beta_I', 
                                          'alpha', 'beta',
                                          'meanlog',  'sdlog', 
                                          'w_0', 'k', 
                                          'b1', 'b2', 'b3', 'b4', 'b5',
                                          'Mean', 'Min', 'R0'))
mcmcEffSummary <- mcmcEffSummary[order(mcmcEffSummary$param, 
                                       mcmcEffSummary$fitType),]

supp_tab_3 <- reshape(subset(mcmcEffSummary)[,c('fitType', 'param', 'mcmcEff')], 
                      timevar = "fitType",
                      idvar = c("param"),
                      direction = "wide")
colnames(supp_tab_3) <- c('', levels(mcmcEffSummary$fitType))
supp_tab_3



################################################################################
# Supplemental Figure 5-10: Posterior means and 95% CIs for parameters
# Supp fig 5: exponential
# Supp fig 6: path-specific
# Supp fig 7: gamma pdf
# Supp fig 8: log-normal pdf
# Supp fig 9: logistic decay
# Supp fig 10: basis splines

postParamsAll <- postParamsEst

# Supp fig 6: exponential
postParamsAll$fitType <- postParamsAll$iddFun
postParamsAll$fitType[postParamsAll$infPeriodSpec == 'exp'] <- 'exp'
postParamsAll$fitType[postParamsAll$infPeriodSpec == 'PS'] <- 'PS'

postParamsAll$fitType <- factor(postParamsAll$fitType, 
                                levels = c('exp', 'PS', 'dgammaIDD', 
                                           'dlnormIDD', 'logitIDD', 'splineIDD'),
                                labels = c('Exponential','Path-specific', 
                                           'Gamma pdf', 'Log-normal pdf',
                                           'Logistic Decay', 'Basis Spline'))


# rename parameters to match manuscript 
postParamsAll$param[postParamsAll$param == 'rateE'] <- 'rho'
postParamsAll$param[postParamsAll$param == 'rateI'] <- 'gamma'

postParamsAll$param[postParamsAll$param == 'shape' & 
                        postParamsAll$fitType == 'Path-specific'] <- 'alpha_I'
postParamsAll$param[postParamsAll$param == 'rate' & 
                        postParamsAll$fitType == 'Path-specific'] <- 'beta_I'

postParamsAll$param[postParamsAll$param == 'shape' & 
                        postParamsAll$fitType == 'Gamma pdf'] <- 'alpha'
postParamsAll$param[postParamsAll$param == 'rate' & 
                        postParamsAll$fitType == 'Gamma pdf'] <- 'beta'

postParamsAll$param[postParamsAll$param == 'mid'] <- 'w_0'
postParamsAll$param[postParamsAll$param == 'rate' & 
                        postParamsAll$fitType == 'Logistic Decay'] <- 'k'

iddPeakParam <- data.frame(datGen = 'IDD_peak',
                           param = c('beta1', 'beta2', 'rho', 'alpha', 'beta'),
                           truth = c(0.25, -7, 0.2, 4, 1))
iddExpParam <- data.frame(datGen = 'IDD_exp',
                          param = c('beta1', 'beta2', 'rho', 'alpha', 'beta'),
                          truth = c(0.4, -7, 0.2, 0.9, 0.3))
iddLogitParam <- data.frame(datGen = 'IDD_logit',
                            param = c('beta1', 'beta2', 'rho', 'w_0', 'k'),
                            truth = c(-1.77, -7, 0.2, 8, 1.5))
psParam <- data.frame(datGen = 'PS',
                      param = c('beta1', 'beta2', 'rho', 'alpha_I', 'beta_I'),
                      truth = c(-1.77, -7, 0.2, 56, 7))

trueParams <- rbind.data.frame(iddPeakParam, 
                               iddExpParam, 
                               iddLogitParam, 
                               psParam)


postParamsAll <- merge(postParamsAll, trueParams, 
                       by = c('datGen', 'param'),
                       all.x = T)

postParamsAll$datGen <- factor(postParamsAll$datGen, 
                               levels = c('IDD_peak', 'IDD_exp', 'IDD_logit', 'PS'),
                               labels = paste0('Data: ',
                                               c('IDD Peak', 'IDD Exp',
                                                 'IDD Logit', 'PS')))

postParamsAll$param <- factor(postParamsAll$param,
                              levels = c('beta1', 'beta2', 'rho', 'gamma',
                                         'alpha_I', 'beta_I', 'alpha', 'beta',
                                         'meanlog', 'sdlog', 'k', 'w_0', 
                                         'b1', 'b2', 'b3', 'b4', 'b5'),
                              labels = c('beta[1]', 'beta[2]', 'rho', 'gamma',
                                         'alpha[I]', 'beta[I]', 'alpha', 'beta',
                                         expression("log("~mu~")"),
                                         expression("log("~sigma~")"), 
                                         'k', 'w[0]', 
                                         'b[1]', 'b[2]', 'b[3]', 'b[4]', 'b[5]'))


maxInfLabs <- c('15-day inf. period',
                '20-day inf. period')
names(maxInfLabs) <- c('15', '20')
postParamsAll$maxInf <- as.character(postParamsAll$maxInf)

# Exponential
pdf('../figures/supplemental/sim_paramEst_Exp_suppfig5.pdf', 
    width = 14, height = 8)
ggplot(subset(postParamsAll, infPeriodSpec == 'exp'), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper)) +
    geom_errorbar(width = 2, col = 'grey25') +
    geom_point(size= 1.5, color = 'black') + 
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, size = 0.8) +
    facet_grid(param~datGen, scales = 'free',
               labeller = labeller(param = label_parsed)) + 
    myTheme + 
    guides(col = 'none') +
    ggtitle('Posterior means and 95% credible intervals\nModel: Exponential') +
    labs(x = 'Simulation Number', y = 'Value') +
    theme(plot.title = element_text(size = 18, h = 0.5),
          strip.text = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
dev.off()

# Path-specific
pdf('../figures/supplemental/sim_paramEst_PS_suppfig6.pdf', 
    width = 24, height = 16)
ggplot(subset(postParamsAll, infPeriodSpec == 'PS'), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper,
           col = maxInf)) +
    geom_errorbar(width = 3) +
    geom_point(size= 1.5, color = 'black') + 
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, size = 0.8) +
    facet_grid(param~datGen + maxInf, scales = 'free',
               labeller = labeller(maxInf = maxInfLabs,
                                   param = label_parsed)) + 
    myTheme + 
    scale_color_manual(values = c('grey25', 'royalblue1')) +
    guides(col = 'none') +
    ggtitle('Posterior means and 95% credible intervals\nModel: Path-specific') +
    labs(x = 'Simulation Number', y = 'Value') +
    theme(plot.title = element_text(size = 28, h = 0.5),
          strip.text = element_text(size = 24),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
dev.off()

# IDD Gamma
pdf('../figures/supplemental/sim_paramEst_IDDGamma_suppfig7.pdf', 
    width = 24, height = 16)
ggplot(subset(postParamsAll, iddFun == 'dgammaIDD'), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper,
           col = maxInf)) +
    geom_errorbar(width = 3) +
    geom_point(size= 1.5, color = 'black') + 
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, size = 0.8) +
    facet_grid(param~datGen + maxInf, scales = 'free',
               labeller = labeller(maxInf = maxInfLabs,
                                   param = label_parsed)) + 
    myTheme + 
    scale_color_manual(values = c('grey25', 'royalblue1')) +
    guides(col = 'none') +
    ggtitle('Posterior means and 95% credible intervals\nModel: Gamma pdf') +
    labs(x = 'Simulation Number', y = 'Value')  +
    theme(plot.title = element_text(size = 28, h = 0.5),
          strip.text = element_text(size = 24),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
dev.off()

# IDD Log-normal
pdf('../figures/supplemental/sim_paramEst_IDDLnorm_suppfig8.pdf', 
    width = 24, height = 16)
ggplot(subset(postParamsAll, iddFun == 'dlnormIDD'), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper,
           col = maxInf)) +
    geom_errorbar(width = 3) +
    geom_point(size= 1.5, color = 'black') + 
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, size = 0.8) +
    facet_grid(param~datGen + maxInf, scales = 'free',
               labeller = labeller(maxInf = maxInfLabs,
                                   param = label_parsed)) + 
    myTheme + 
    scale_color_manual(values = c('grey25', 'royalblue1')) +
    guides(col = 'none') +
    ggtitle('Posterior means and 95% credible intervals\nModel: Log-normal pdf') +
    labs(x = 'Simulation Number', y = 'Value') +
    theme(plot.title = element_text(size = 28, h = 0.5),
          strip.text = element_text(size = 24),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
dev.off()

# IDD Logit
pdf('../figures/supplemental/sim_paramEst_IDDLogit_suppfig9.pdf', 
    width = 24, height = 16)
ggplot(subset(postParamsAll, iddFun == 'logitIDD'), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper,
           col = maxInf)) +
    geom_errorbar(width = 3) +
    geom_point(size= 1.5, color = 'black') + 
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, size = 0.8) +
    facet_grid(param~datGen + maxInf, scales = 'free',
               labeller = labeller(maxInf = maxInfLabs,
                                   param = label_parsed)) + 
    myTheme + 
    scale_color_manual(values = c('grey25', 'royalblue1')) +
    guides(col = 'none') +
    ggtitle('Posterior means and 95% credible intervals\nModel: Logistic Decay') +
    labs(x = 'Simulation Number', y = 'Value') +
    theme(plot.title = element_text(size = 28, h = 0.5),
          strip.text = element_text(size = 24),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
dev.off()


# IDD Spline
pdf('../figures/supplemental/sim_paramEst_IDDSpline_suppfig10.pdf', 
    width = 24, height = 22)
ggplot(subset(postParamsAll, iddFun == 'splineIDD'), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper,
           col = maxInf)) +
    geom_errorbar(width = 3) +
    geom_point(size= 1.5, color = 'black') + 
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, size = 0.8) +
    facet_grid(param~datGen + maxInf, scales = 'free',
               labeller = labeller(maxInf = maxInfLabs,
                                   param = label_parsed)) + 
    myTheme + 
    scale_color_manual(values = c('grey25', 'royalblue1')) +
    guides(col = 'none') +
    ggtitle('Posterior means and 95% credible intervals\nModel: Basis Splines') +
    labs(x = 'Simulation Number', y = 'Value')  +
    theme(plot.title = element_text(size = 28, h = 0.5),
          strip.text = element_text(size = 24),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18))
dev.off()

