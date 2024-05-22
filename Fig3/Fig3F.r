info = read.table('meta_spatial.csv', row.names = 1, sep = '\t', header = T)
df = as.data.frame(table(info$sample, info$AF))

library(ggplot2)
library(ggpubr)

colnames(df) = c('slide', 'AF', 'value')
# df$slide = factor(df$slide, c(levels = 'SNU38A', 
#                   'SNU19A', 
#                   'SNU22A', 
#                   'SNU40A', 
#                   'SNU16A', 
#                   'SNU16B', 
#                   'SNU25A', 
#                   'SNU25B', 
#                   'SNU27A', 
#                   'SNU27B', 
#                   'SNU17A', 
#                   'SNU17B', 
#                   'SNU18A', 
#                   'SNU18B', 
#                   'SNU18Afrozen', 
#                   'SNU46A', 
#                   'SNU46B', 
#                   'SNU23A', 
#                   'SNU51A', 
#                   'SNU51B', 
#                   'SNU21A', 
#                   'SNU21B', 
#                   'SNU21Afrozen', 
#                   'SNU21Bfrozen', 
#                   'SNU21Cfrozen', 
#                   'SNU21Dfrozen', 
#                   'SNU43A', 
#                   'SNU24A', 
#                   'SNU24B', 
#                   'SNU33A', 
#                   'SNU33B', 
#                   'SNU34A'))
df$slide = factor(df$slide, levels = c('SNU38',
                                       'SNU19',
                                       'SNU22',
                                       'SNU40',
                                       'SNU16',
                                       'SNU25',
                                       'SNU27',
                                       'SNU17',
                                       'SNU18',
                                       'SNU46',
                                       'SNU23',
                                       'SNU51',
                                       'SNU21',
                                       'SNU43',
                                       'SNU24',
                                       'SNU33',
                                       'SNU34'))

meta = read.csv('infoAnno2.csv', row.names = 1)
meta = meta[levels(df$slide), ]
#df$slide = factor(df$slide, levels = meta[,1])

identical(rownames(meta), levels(df[,1]))

df$AF = factor(df$AF, levels = c('LE_GM', 'LE_WGM', 'LE_WM', 'CT_control',
                                 'CT', 'PNZ_ccGSC', 'PAN', 'HBV_ccGSC',
                                 'MVP'))

p1 = ggplot(df, aes(x = slide, y = value, fill = AF)) +
  geom_bar(position = 'fill', stat = 'identity') + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c('LE_GM' = 'darkolivegreen4',
                               'LE_WGM' = 'chartreuse',
                               'LE_WM' = 'goldenrod2',
                               'CT_control' = 'yellow',
                               'CT' = 'salmon',
                               'PNZ_ccGSC' = 'navy',
                               'PAN' = 'dodgerblue',
                               'HBV_ccGSC' = 'magenta',
                               'MVP' = 'darkmagenta'))
