library(Seurat)
obj = readRDS('/media/user/disk21/completeAnalysis/visium_15Sept/visiumObj_20Sept.rds')
meta = read.table('/media/user/disk21/completeAnalysis/cellbrowser/Spatial-Data/meta_visium.csv', row.names = 1,
                  sep = '\t', header = T)

obj@meta.data = meta

obj@meta.data$AF = factor(obj@meta.data$AF,
                          levels = c('LE_GM', 'LE_WGM', 'LE_WM',
                                     'CT_control', 'CT', 'PNZ_ccGSC',
                                     'PAN', 'HBV_ccGSC', 'MVP'))

pdf('visiumDotplot_AF.pdf', width = 10, height = 7.5)

p = DotPlot(obj, features = c('SLC17A7', 'SNAP25', 'NRGN', 'MOBP', 'PLP1',
                              'MBP', 'GALNT13', 'ERMN', 'SNX22',
                              'TOP2A', 'DLL3', 'SOX11', 'CHI3L1',
                              'TIMP1', 'NAMPT', 'VEGFA', 'ADM',
                              'HILPDA', 'CD163', 'IFI30', 'LYZ',
                              'COL1A2', 'ACTA2', 'PLVAP'), group.by = 'AF') +
  RotatedAxis() + NoLegend()
p = p +
  theme(axis.text = element_text(color = 'black', size = 20))
print(p)
dev.off()  

library(ggplot2)

ggplot(meta, aes(x = AF_ivy, y = inferCNV, fill = AF_ivy)) + geom_boxplot() + 
  RotatedAxis() + NoLegend() + 
  geom_hline(yintercept = c(100, 150), linetype = 'dashed')