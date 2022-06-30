################################################################################
# Results of the Ebola Analysis
################################################################################

library(BayesSEIR)
library(ABSEIR)
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
### Figure 5: 1995 DRC Ebola epidemic

ebola <- ABSEIR::Kikwit1995

ebola <- ebola[1:135,]

# pdf('../figures/ebola_data_fig5.pdf', width = 6, height = 4)
plot(as.Date(ebola$Date), ebola$Count, type='h', 
     main='Daily Counts of Ebola by Symptom Onset Date', ylab='New Cases',
     xlab = 'Symptom Onset Date', lwd = 2, col = 'grey30', 
     cex.main = 1.3, cex.axis = 1.1, cex.lab = 1.2)
abline(v = as.Date('1995-05-09'), col = 'blue', lty =2, lwd = 3)
legend('topright', 'Date of\nIntervention', col = 'blue', lty =2, lwd = 3, 
       cex = 1.1, bty = 'n')
# dev.off()



################################################################################
### Combine batch files

batchFiles <- list.files('./output/')

# combine results with est infectious period
batch_i <- readRDS(paste0('./output/', batchFiles[1]))

gdiag <- batch_i$gdiag
paramsSummary <- batch_i$paramsSummary
iddSummary <- batch_i$iddSummary
r0Summary <- batch_i$r0Summary
waicSummary <- batch_i$waicSummary

for (i in 2:length(batchFiles)) {
    batch_i <- readRDS(paste0('./output/', batchFiles[i]))
    
    gdiag <-rbind.data.frame(gdiag, batch_i$gdiag)
    paramsSummary <-rbind.data.frame(paramsSummary, batch_i$paramsSummary)
    iddSummary <-rbind.data.frame(iddSummary, batch_i$iddSummary)
    r0Summary <-rbind.data.frame(r0Summary, batch_i$r0Summary)
    waicSummary <-rbind.data.frame(waicSummary, batch_i$waicSummary)
}


################################################################################
### Check convergence

# should be 0 if everything converged
sum(gdiag$gr > 1.1)

################################################################################
### Figure 6: Posterior mean and 95% CI for R0(t)

r0Summary$fitType <- r0Summary$iddFun
r0Summary$fitType[r0Summary$infPeriodSpec == 'exp'] <- 'exp'
r0Summary$fitType[r0Summary$infPeriodSpec == 'PS'] <- 'PS'

r0Summary$fitType <- factor(r0Summary$fitType, 
                            levels = c('exp', 'PS', 'dgammaIDD', 'dlnormIDD',
                                       'logitIDD', 'splineIDD'),
                            labels = c('Exponential', 'Path-specific', 
                                       'IDD - Gamma pdf', 'IDD - Log-normal pdf',
                                       'IDD - Logistic Decay', 'IDD - Basis Spline'))


# merge in dates
ebola$time <- 1:nrow(ebola) - 3

r0Summary <- merge(r0Summary, ebola, by = 'time', all.x = T)
r0Summary <- r0Summary[-which(is.na(r0Summary$Date)),]

pal <- brewer.pal(6, 'Set2')

# pdf('../figures/main/ebola_r0_fig6.pdf', height = 6, width = 10)
ggplot(r0Summary, aes(x = Date, y = mean, 
                      ymin = lower, ymax = upper,
                      col = fitType)) +
    geom_line(size = 1.5) + 
    geom_ribbon(aes(fill = fitType), alpha = 0.1) + 
    geom_hline(yintercept = 1, size = 1, linetype = 2) +
    facet_wrap(~fitType) +
    myTheme +  
    ggtitle(expression('Posterior mean and 95% CI of' ~R[0](t))) +
    labs(x = 'Epidemic Time', y = expression(R[0](t))) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    guides(fill = 'none', col = 'none')
# dev.off()

# which day does R0 cross 1
expR0 <- r0Summary[r0Summary$infPeriodSpec == 'exp',]
expR0[c(max(which(expR0$mean > 1)), max(which(expR0$mean > 1)) + 1),]

psR0 <- r0Summary[r0Summary$infPeriodSpec == 'PS',]
psR0[c(max(which(psR0$mean > 1)), max(which(psR0$mean > 1)) + 1),]

gammaR0 <- r0Summary[which(r0Summary$iddFun == 'dgammaIDD'),]
gammaR0[c(max(which(gammaR0$mean > 1)), max(which(gammaR0$mean > 1)) + 1),]

lnormR0 <- r0Summary[which(r0Summary$iddFun == 'dlnormIDD'),]
lnormR0[c(max(which(lnormR0$mean > 1)), max(which(lnormR0$mean > 1)) + 1),]

logitR0 <- r0Summary[which(r0Summary$iddFun == 'logitIDD'),]
logitR0[c(max(which(logitR0$mean > 1)), max(which(logitR0$mean > 1)) + 1),]

splinesR0 <- r0Summary[which(r0Summary$iddFun == 'splineIDD'),]
splinesR0[c(max(which(splinesR0$mean > 1)), max(which(splinesR0$mean > 1)) + 1),]

################################################################################
### Figure 7: Posterior mean and 95% CI for IDD curves

iddSummary <- iddSummary[!is.na(iddSummary$iddFun),]

iddSummary$iddFun <- factor(iddSummary$iddFun, 
                            levels = c('dgammaIDD', 'dlnormIDD',
                                       'logitIDD', 'splineIDD'),
                            labels = c('Gamma pdf', 'Log-normal pdf',
                                       'Logistic Decay', 'Basis Spline'))

iddCurvePal <- pal[3:6]

# pdf('../figures/ebola_iddCurves_fig7.pdf', height = 4, width = 10)
ggplot(iddSummary, aes(x = infDay, y = median, 
                       ymin = lower, ymax = upper,
                       col = iddFun)) +
    geom_line(size= 1.3) + 
    geom_ribbon(aes(fill = iddFun), alpha = 0.3) + 
    facet_wrap(~iddFun, nrow = 1) +
    myTheme +  
    labs(x = 'Days Infectious', y = expression(pi[0]^(SE))) +
    scale_color_manual(values = iddCurvePal) +
    scale_fill_manual(values = iddCurvePal) +
    guides(fill = 'none', col = 'none') +
    ggtitle(expression('Posterior median and 95% CI of '~pi[0]^(SE)))
# dev.off()

################################################################################
### Table 2: WAIC by model

waicSummary$fitType <- waicSummary$iddFun
waicSummary$fitType[waicSummary$infPeriodSpec == 'exp'] <- 'exp'
waicSummary$fitType[waicSummary$infPeriodSpec == 'PS'] <- 'PS'

waicSummary$fitType <- factor(waicSummary$fitType, 
                              levels = c('exp', 'PS', 'dgammaIDD', 'dlnormIDD',
                                         'logitIDD', 'splineIDD'),
                              labels = c('Exponential', 'Path-specific', 
                                         'Gamma pdf', 'Log-normal pdf',
                                         'Logistic Decay', 'Basis Spline'))

tab5 <- t(trimws(format(round(waicSummary$waic, 2), nsmall = 2)))
row.names(tab5) <- 'WAIC'
colnames(tab5) <- waicSummary$fitType

tab5

################################################################################
### Supplemental Table: Posterior mean and 95% CIs for each parameter

paramsSummary$fitType <- paramsSummary$iddFun
paramsSummary$fitType[paramsSummary$infPeriodSpec == 'exp'] <- 'exp'
paramsSummary$fitType[paramsSummary$infPeriodSpec == 'PS'] <- 'PS'

paramsSummary$fitType <- factor(paramsSummary$fitType, 
                                levels = c('exp', 'PS', 'dgammaIDD', 'dlnormIDD',
                                           'logitIDD', 'splineIDD'),
                                labels = c('Exponential', 'Path-specific', 
                                           'Gamma pdf', 'Log-normal pdf',
                                           'Logistic Decay', 'Basis Spline'))

paramsSummary$mean <- trimws(format(round(paramsSummary$mean, 2), nsmall = 2))
paramsSummary$ci <- paste0('(', 
                           trimws(format(round(paramsSummary$lower, 2), nsmall = 2)),
                           ', ', 
                           trimws(format(round(paramsSummary$upper, 2), nsmall = 2)),
                           ')')

paramsSummary$postSum <- paste0(paramsSummary$mean, ' ', paramsSummary$ci)

# rename parameters to match table 
paramsSummary$param[paramsSummary$param == 'rateE'] <- 'rho'
paramsSummary$param[paramsSummary$param == 'rateI'] <- 'gamma'

paramsSummary$param[paramsSummary$param == 'shape' & 
                        paramsSummary$fitType == 'Path-specific'] <- 'alpha_I'
paramsSummary$param[paramsSummary$param == 'rate' & 
                        paramsSummary$fitType == 'Path-specific'] <- 'beta_I'

paramsSummary$param[paramsSummary$param == 'shape' & 
                        paramsSummary$fitType == 'Gamma pdf'] <- 'alpha'
paramsSummary$param[paramsSummary$param == 'rate' & 
                        paramsSummary$fitType == 'Gamma pdf'] <- 'beta'

paramsSummary$param[paramsSummary$param == 'mid'] <- 'w_0'
paramsSummary$param[paramsSummary$param == 'rate' & 
                        paramsSummary$fitType == 'Logistic Decay'] <- 'k'


tabLong <- paramsSummary[,c('fitType', 'param', 'postSum')]

supp_tab_4 <- reshape(tabLong, 
             timevar = "fitType",
             idvar = c("param"),
             direction = "wide")

# add in WAIC and posterior R0 for each model

postR0 <- r0Summary[r0Summary$time == 1,]
postR0$mean <- format(round(postR0$mean, 2), nsmall = 2)
postR0$ci <- paste0('(', trimws(format(round(postR0$lower, 2), nsmall = 2)), ', ', 
                    trimws(format(round(postR0$upper, 2), nsmall = 2)), ')')

postR0$postSum <- paste0(postR0$mean, ' ', postR0$ci)

supp_tab_4 <- rbind(c('WAIC', format(round(waicSummary$waic, 2), nsmall = 2)),
                    c('R0', postR0$postSum),
                    supp_tab_4)
colnames(supp_tab_4) <- c('', levels(paramsSummary$fitType))

supp_tab_4
