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
  theme(plot.title = element_text(h = 0.5, size = 16),
        legend.title = element_text(size = 14, h = 0.5),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background=element_rect(fill="white"),
        strip.text = element_text(size = 14))

source('../helper_functions.R')


################################################################################
### Combine batch files

batchFiles <- list.files('./batch_output/')
batchFilesKnown <- batchFiles[grep('knownE', batchFiles)]
batchFilesEst <- batchFiles[grep('estimatedE', batchFiles)]


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
pdf('./figures/Figure1.pdf', height = 3.5, width = 9.5)
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
dev.off()


################################################################################
### Figure 2: Posterior median IDD transmissibility curves compared to truth

trueCurves <- data.frame(datGen = rep(c('IDD_peak', 'IDD_exp', 'IDD_logit'), 
                                      each = maxInf),
                         infDay = 1:maxInf,
                         iddCurveTruth = c(iddPeakCurve, iddExpCurve, iddLogitCurve))

iddCurveAll <- merge(iddCurveAll, trueCurves, by = c('datGen', 'infDay'), all.x = T)


iddCurveAll$iddFun <- factor(iddCurveAll$iddFun,
                             levels = c('dgammaIDD', 'dlnormIDD',
                                        'logitIDD',  'splineIDD'),
                             labels = c('Gamma pdf', 'Log-normal pdf', 
                                        'Logistic Decay', 'Basis Spline'))

iddCurveAll$datGen <- factor(iddCurveAll$datGen,
                             levels = c('IDD_peak', 'IDD_exp', 'IDD_logit'),
                             labels = c('IDD Peak', 'IDD Exp', 'IDD Logit'))

# when exposure times are estimated, we only care about certain scenarios for
# estimating IDD curves

# iddCurveAll <- iddCurveAll[-which(iddCurveAll$datGen == 'IDD Peak' & 
#                                       iddCurveAll$iddFun == 'Logistic Decay'),]
iddCurveAll <- iddCurveAll[-which(iddCurveAll$datGen == 'IDD Exp' & 
                                    iddCurveAll$iddFun == 'Logistic Decay'),]
iddCurveAll <- iddCurveAll[-which(iddCurveAll$datGen == 'IDD Logit' & 
                                    iddCurveAll$iddFun %in% 
                                    c('Gamma pdf', 'Log-normal pdf')),]

pal <- c('magenta3', 'darkorange', 'green3')

p1 <- ggplot(data = subset(iddCurveAll, 
                           iddFun != 'Basis Spline' & 
                             EstarType == 'Known' & 
                             maxInf == 20)) +
  geom_line(aes(x = infDay, y = median, group = simNumber),
            color =  adjustcolor('grey', alpha = 0.5)) +
  geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1) +
  facet_wrap(~datGen+iddFun, nrow = 1) + 
  myTheme +
  labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
       title = 'Known Exposure Times')  +
  scale_color_manual(values =pal) +
  theme(legend.position = "none")

p2 <- ggplot(data = subset(iddCurveAll, 
                           iddFun != 'Basis Spline' & 
                             EstarType == 'Estimated' & 
                             maxInf == 20)) +
  geom_line(aes(x = infDay, y = median, group = simNumber),
            color =  adjustcolor('grey', alpha = 0.5)) +
  geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1) +
  facet_wrap(~datGen+iddFun, nrow = 1) + 
  myTheme +
  labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
       title = 'Estimated Exposure Times')  +
  scale_color_manual(values =pal) +
  theme(legend.position = "none")

pdf('./figures/Figure3.pdf', height = 7, width = 10)
grid.arrange(p1, p2,
             top = textGrob(expression('Posterior median estimates of '~pi[0]^(SE)),
                            gp = gpar(fontsize = 18, font = 2)))
dev.off()


################################################################################
# Figure 3: Posterior IDD curves using basis splines compared to truth


p1 <- ggplot(data = subset(iddCurveAll, 
                           iddFun == 'Basis Spline' & 
                             EstarType == 'Known' & 
                             maxInf == 20)) +
  geom_line(aes(x = infDay, y = median, group = simNumber),
            color =  adjustcolor('grey', alpha = 0.5)) +
  geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1) +
  facet_wrap(~datGen+iddFun, nrow = 1) + 
  myTheme +
  labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
       title = 'Known Exposure Times')  +
  scale_color_manual(values =pal) +
  theme(legend.position = "none")

p2 <- ggplot(data = subset(iddCurveAll, 
                           iddFun == 'Basis Spline' & 
                             EstarType == 'Estimated' & 
                             maxInf == 20)) +
  geom_line(aes(x = infDay, y = median, group = simNumber),
            color =  adjustcolor('grey', alpha = 0.5)) +
  geom_line(aes(x = infDay, y = iddCurveTruth, col = datGen), size = 1) +
  facet_wrap(~datGen+iddFun, nrow = 1) + 
  myTheme +
  labs(x = 'Days Infectious', y = expression(pi[0]^"(SE)"),
       title = 'Estimated Exposure Times')  +
  scale_color_manual(values =pal) +
  theme(legend.position = "none")

pdf('./figures/Figure3.pdf', height = 7, width = 10)
grid.arrange(p1, p2,
             top = textGrob(expression('Posterior median estimates of '~pi[0]^(SE)),
                            gp = gpar(fontsize = 18, font = 2)))
dev.off()

################################################################################
# Figure 4: MSE of R0

# get true R0
maxInf <- 15
N <- 5363500
fullX <- getX(250, 120)

# get true R0's
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

iddPSR0 <- getR0(infPeriodSpec = 'PS', 
                 beta = c(-1.77, -7), X = fullX, N = N,
                 infPSParams = list(dist = 'gamma',
                                    psParams = list(shape = 56, rate = 7),
                                    maxInf = maxInf))

trueR0 <- data.frame(datGen = rep(c('IDD_peak', 'IDD_exp', 'IDD_logit', 'PS'), 
                                  each = 250),
                     time = 1:250,
                     truth = c(iddPeakR0, iddExpR0, iddLogitR0, iddPSR0))

# merge estimates with truth
r0Est <- merge(r0Est, trueR0, by = c('datGen', 'time'), all.x = T)

r0MSE <- ddply(subset(r0Est, time == 1), .(datGen, infPeriodSpec, iddFun, maxInf),
               summarize,
               mse = mean((mean-truth)^2))

# barplot
r0MSE$fitType <- r0MSE$iddFun
r0MSE$fitType[r0MSE$infPeriodSpec == 'exp'] <- 'exp'
r0MSE$fitType[r0MSE$infPeriodSpec == 'PS'] <- 'PS'

r0MSE$fitType <- factor(r0MSE$fitType, 
                        levels = c('dlnormIDD', 'dgammaIDD',
                                   'logitIDD', 'exp', 
                                   'PS', 'splineIDD'),
                        labels = c(expression(`Log-normal`~Pdf^"\206"),
                                   expression("Gamma"~Pdf^"\206"), 
                                   expression(Logistic~Decay^"\206"),
                                   expression(Exponential),
                                   expression(`Path-specific`), 
                                   expression(Basis~Spline^"\206")))


pal <- brewer.pal(6, 'Dark2')
pal <- pal[c(1:3, 6, 4:5)]

makeFootnote <- function(footnoteText,
                         size = .7, color = grey(.5),
                         xunit=0.85, yunit = 35) {
  require(grid)
  pushViewport(viewport())
  grid.text(label = footnoteText ,
            x = unit(xunit,"npc"),
            y = unit(yunit, "mm"),
            just = c("left", "bottom"),
            gp = gpar(cex = size, col = color))
  
  popViewport()
}

r0MSE$datGen <- factor(r0MSE$datGen,
                       levels = c('IDD_peak', 'IDD_exp', 
                                  'IDD_logit', 'PS'),
                       labels = c('IDD Peak', 'IDD Exp', 
                                  'IDD Logit', 'Path-specific'))

maxInfLabs <- c('15 day infectious period',
                '20 day infectious period')
names(maxInfLabs) <- c('15', '20')
r0MSE$maxInf <- as.character(r0MSE$maxInf)

pdf('./figures/Figure4.pdf', height = 6, width = 12)
ggplot(data = r0MSE, 
       aes(x = datGen, y = mse, fill = fitType, group = fitType)) +
  geom_bar(position="dodge", stat="identity", col = 'black', size = 0.2) +
  facet_wrap(~maxInf, nrow = 2,
             labeller = labeller(maxInf = maxInfLabs)) +
  scale_fill_manual(values=pal, labels = parse_format()) +
  myTheme +
  labs(title = expression('MSE of'~ R[0]),
       fill = 'Modeling Approach', y = 'MSE', x = 'Data Type') +
  theme(legend.text.align = 0)
makeFootnote(expression("\206"), 
             color = "black", size=1, xunit = 0.777, yunit = 35)
makeFootnote("models fit using the proposed\nIDD transmissibility approach", 
             color = "black", size=1.2, xunit = 0.787, yunit = 27)
dev.off()



################################################################################
# Table 1: Mean, minimum, R0 MCMC efficiency

mcmcEffEst$fitType <- mcmcEffEst$iddFun
mcmcEffEst$fitType[mcmcEffEst$infPeriodSpec == 'exp'] <- 'exp'
mcmcEffEst$fitType[mcmcEffEst$infPeriodSpec == 'PS'] <- 'PS'

ddply(mcmcEffEst, .(maxInf, fitType, param), summarize,
      meanEff = mean(eff),
      sdEff = sd(eff))

################################################################################
# Table 2: Mean, SD of WAIC

waicEst$fitType <- waicEst$iddFun
waicEst$fitType[waicEst$infPeriodSpec == 'exp'] <- 'exp'
waicEst$fitType[waicEst$infPeriodSpec == 'PS'] <- 'PS'

ggplot(waicEst, aes(x = datGen, y = waic, fill = fitType)) +
  geom_boxplot() +
  myTheme

waicSummary <- ddply(waicEst, .(datGen, fitType), summarize,
                     mean = mean(waic),
                     sd = sd(waic))


waicSummary <- waicSummary[order(waicSummary$datGen, waicSummary$mean),]


################################################################################
# Supplemental Figure 1: Posterior mean R0(t) compared to truth

r0Est$fitType <- r0Est$iddFun
r0Est$fitType[r0Est$infPeriodSpec == 'exp'] <- 'exp'
r0Est$fitType[r0Est$infPeriodSpec == 'PS'] <- 'PS'

r0Est$datGen <- factor(r0Est$datGen, 
                       levels = c('IDD_peak', 'IDD_exp', 'IDD_logit', 'PS'),
                       labels = c('IDD Peak', 'IDD Exp', 'IDD Logit', 'PS'))

r0Est$fitType <- factor(r0Est$fitType, 
                        levels = c('exp', 'PS', 'dgammaIDD', 'dlnormIDD', 'logitIDD', 'splineIDD'),
                        labels = c('Exponential','Path-specific', 
                                   'Gamma Pdf', 'Log-normal Pdf', 'Logistic Decay', 'Basis Spline'))


pal <- c('magenta3', 'darkorange', 'green3', 'royalblue2')

ggplot(subset(r0Est, maxInf == 15), aes(x = time)) + 
  geom_line(aes(y = mean), color =  adjustcolor('grey', alpha = 0.5)) +
  geom_line(aes(y = truth, col = datGen), size = 1) +
  facet_wrap(~datGen + fitType, nrow = 4) +
  myTheme +
  labs(x = 'Epidemic Time', y = expression(R[0](t)),
       col = 'Data Generation\nScenario')  +
  scale_color_manual(values =pal) +
  ggtitle(expression('Posterior mean estimates of '~R[0](t))) +
  theme(legend.position = "none")

ggplot(subset(r0Est, maxInf == 20), aes(x = time)) + 
  geom_line(aes(y = mean), color =  adjustcolor('grey', alpha = 0.5)) +
  geom_line(aes(y = truth, col = datGen), size = 1) +
  facet_wrap(~datGen + fitType, nrow = 4) +
  myTheme +
  labs(x = 'Epidemic Time', y = expression(R[0](t)),
       col = 'Data Generation\nScenario')  +
  scale_color_manual(values =pal) +
  ggtitle(expression('Posterior mean estimates of '~R[0](t))) +
  theme(legend.position = "none")



################################################################################
# Supplemental Figure 2: Posterior means and 95% CIs for parameters

iddPeakParam <- data.frame(datGen = 'IDD_peak',
                           param = c('beta1', 'beta2', 'rateE', 'shape', 'rate'),
                           truth = c(0.25, -7, 0.2, 4, 1))
iddExpParam <- data.frame(datGen = 'IDD_exp',
                          param = c('beta1', 'beta2', 'rateE', 'shape', 'rate'),
                          truth = c(0.4, -7, 0.2, 0.9, 0.3))
iddLogitParam <- data.frame(datGen = 'IDD_logit',
                            param = c('beta1', 'beta2', 'rateE', 'mid', 'rate'),
                            truth = c(-1.77, -7, 0.2, 8, 1.5))
psParam <- data.frame(datGen = 'PS',
                      param = c('beta1', 'beta2', 'rateE', 'shape', 'rate'),
                      truth = c(-1.77, -7, 0.2, 56, 7))

trueParams <- rbind.data.frame(iddPeakParam, 
                               iddExpParam, 
                               iddLogitParam, 
                               psParam)

postParamsEst <- merge(postParamsEst, trueParams, 
                       by = c('datGen', 'param'),
                       all.x = T)


# IDD Gamma
ggplot(subset(postParamsEst, iddFun == 'dgammaIDD' & maxInf == 15), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
  facet_wrap(~ datGen + param, ncol = 5, scales = 'free') + 
  myTheme

ggplot(subset(postParamsEst, iddFun == 'dgammaIDD'& maxInf == 20), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
  facet_wrap(~ datGen + param, ncol = 5, scales = 'free') + 
  myTheme

# IDD log-normal
ggplot(subset(postParamsEst, iddFun == 'dlnormIDD' & maxInf == 15), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
  facet_wrap(~ datGen + param, ncol = 5, scales = 'free') + 
  myTheme

ggplot(subset(postParamsEst, iddFun == 'dlnormIDD'& maxInf == 20), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
  facet_wrap(~ datGen + param, ncol = 5, scales = 'free') + 
  myTheme

# IDD logit
ggplot(subset(postParamsEst, iddFun == 'logitIDD' & maxInf == 15), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
  facet_wrap(~ datGen + param, ncol = 5, scales = 'free') + 
  myTheme

ggplot(subset(postParamsEst, iddFun == 'logitIDD'& maxInf == 20), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
  facet_wrap(~ datGen + param, ncol = 5, scales = 'free') + 
  myTheme

# IDD spline
ggplot(subset(postParamsEst, iddFun == 'dlnormIDD' & maxInf == 15), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
  facet_wrap(~ datGen + param, ncol = 5, scales = 'free') + 
  myTheme

ggplot(subset(postParamsEst, iddFun == 'dlnormIDD'& maxInf == 20), 
       aes(x = simNumber, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  geom_hline(aes(yintercept = truth), col = 'red', linetype = 2) +
  facet_wrap(~ datGen + param, ncol = 5, scales = 'free') + 
  myTheme

