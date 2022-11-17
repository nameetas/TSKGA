library(ggplot2)
library(ggpubr)

umap = read.csv('scRef_1Sep_UMAP.csv', row.names = 1) #UMAP coordinates from scanpy for Reference data
info = read.csv('1Sep/finalCells/info2.csv', row.names = 1) #Reference meta data 

df = cbind.data.frame(umap, info)

p = ggplot(df, aes(x = X0, y = X1, color = riverCelltype)) + 
  geom_point(size = 0.1) + theme_classic() + 
  scale_color_manual(values = c('Neuron' = 'aquamarine4',
                                'Oligo_OPC' = 'plum4',
                                'Astrocyte' = 'orange',
                                'Neoplastic' = 'dodgerblue',
                                'Microglia' = 'purple',
                                'Myeloid' = 'hotpink',
                                'Lymphoid' = 'navy',
                                'Fibroblast' = 'burlywood',
                                'Pericyte' = 'maroon4',
                                'Endothelial' = 'chartreuse',
                                'RBC' = 'olivedrab4')) +
  xlab('UMAP_1') + ylab('UMAP_2') + 
  guides(colour = guide_legend(override.aes = list(size=10)))
p

q = ggplot(df, aes(x = X0, y = X1, color = label1)) + 
  geom_point(size = 0.1) + theme_classic() + 
  xlab('UMAP_1') + ylab('UMAP_2') + 
  guides(colour = guide_legend(override.aes = list(size=10)))
q


pdf('Fig2A.pdf', width = 10, height = 10)
  print(p)
  print(q)
dev.off()
