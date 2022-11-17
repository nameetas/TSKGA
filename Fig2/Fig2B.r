library(Seurat)
library(ggplot2)

meta = read.csv('meta_16Sept.csv', row.names = 1) #Spatial meta data
binary = readRDS('binary_15Sept.rds') #Binary cell type data
identical(rownames(binary), rownames(meta))

meta$numCelltype = apply(binary, 1, sum)
m1 = as.data.frame(table(meta$AF, meta$numCelltype))
colnames(m1) = c('AF', 'numCelltype', 'Freq')

m1$AF = factor(m1$AF, levels = c('LE_GM', 'LE_WGM', 'LE_WM',
                                 'CT_control', 'CT', 'PNZ_ccGSC',
                                 'PAN', 'HBV_ccGSC', 'MVP'))

m1$AF = factor(m1$AF, levels = c('LE_GM', 'LE_WGM', 'LE_WM',
                                 'CT_control', 'CT', 'PNZ_ccGSC',
                                 'PAN', 'HBV_ccGSC', 'MVP'))

m1$numCelltype = factor(m1$numCelltype, 
                        levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 
                                   "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                   "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", 
                                   "30", "31", "32", "33", "34", ">=35"))

pdf('Fig2B.pdf', width = 10, height = 5)
ggplot(m1, aes(x = numCelltype, y = Freq, fill = AF)) + 
  geom_bar(position = 'fill', stat = 'identity') + 
  theme_classic() + xlab('#celltypes in each spot') + 
  ylab('Proportions') +
  scale_fill_manual(values = c('LE_GM' = 'darkolivegreen4',
                               'LE_WGM' = 'chartreuse',
                               'LE_WM' = 'goldenrod2',
                               'CT_control' = 'yellow',
                               'CT' = 'salmon',
                               'PNZ_ccGSC' = 'navy',
                               'PAN' = 'dodgerblue',
                               'HBV_ccGSC' = 'magenta',
                               'MVP' = 'darkmagenta')) + 
  RotatedAxis() + 
  theme(axis.text = element_text(size = 10, color = 'black')) + 
  NoLegend()
dev.off()
