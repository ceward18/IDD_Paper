

library(ggplot2)

modTime <- readRDS('timePerIter.rds')

# format for plot legend
modTime$model <- factor(modTime$model,
                        levels = c('PS', 'exp', 'IDD'),
                       levels = c('Path-specific', 'Exponential', 'IDD'))

# color palette for plot
myPal <- rev(c('royalblue2', 'goldenrod3', 'orangered2'))

pdf("../figures/figure8.pdf", height = 4, width = 8)
ggplot(modTime, aes(x = nInf, y = sec, group = model, 
                   color = model, linetype = model)) +
    geom_point(size = 2.5) + geom_line(size = 1.2) +
    theme_bw() +
    labs(title = 'Iteration time by epidemic size',
         x = 'Epidemic size',
         y = 'Time per iteration (sec)',
         color = 'Modeling Approach',
         linetype = 'Modeling Approach') +
    theme(plot.title = element_text(h = 0.5, size = 16),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.key.width = unit(1.5,"cm")) +
    guides(linetype = guide_legend(override.aes = list(width = 4))) + 
    scale_x_continuous(labels = scales::comma) +
    scale_linetype_manual(values = c(3,5,1)) +
    scale_color_manual(values = myPal)
dev.off()
