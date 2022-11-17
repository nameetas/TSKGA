obj = readRDS('finalObj_citeseq.rds') #Cite-seq seurat object

library(ggplot2)
library(ggpubr)
library(Seurat)

meta = read.csv('meta_17Sept.csv', row.names = 1) #Cite-seq meta data
identical(rownames(obj@meta.data), rownames(meta))
obj@meta.data = meta
unique(obj@meta.data$Major.celltype)

obj1 = subset(obj, subset = Major.celltype != 'misc')

obj1@meta.data$Major.celltype = factor(obj1@meta.data$Major.celltype, 
                                       levels = c('DC',
                                                  'endothelial',
                                                  'granulocyte',
                                                  'lymphocyte',
                                                  'macrophage',
                                                  'microglia',
                                                  'monocyte',
                                                  'oligodendrocyte',
                                                  'stromal cell',
                                                  'tumor cell'))
pdf('Fig1C.pdf', width = 10, height = 10)
  DimPlot(obj1, group.by = 'Major.celltype', cols = c('#FA8072',
                                                    '#4EEE94', 
                                                    '#EEB422',
                                                    '#B22222', 
                                                    '#FF69B4',
                                                    '#A020F0',
                                                    '#000080',
                                                    '#EEAEEE',
                                                    '#DEB887',
                                                    '#1E90FF'))
dev.off()

rownames(obj1@assays$ADT@data) = gsub('-TotalSeqC', '', rownames(obj1@assays$ADT@data))
rownames(obj1@assays$ADT@scale.data) = gsub('-TotalSeqC', '', rownames(obj1@assays$ADT@scale.data))

pdf('Fig1E_citeseqDotplot_majorCelltype.pdf', width = 32, height = 5.5)

p = DotPlot(obj1, features = c('AREG',
                               'PLAC8',
                               'FCER1A',
                               'CLDN5',
                               'VWF',
                               'EGFL7',
                               'CMTM2',
                               'FCGR3B',
                               'IFITM2',
                               'CD3E',
                               'CCL5',
                               'CD52',
                               'GPNMB',
                               'CTSL',
                               'MRC1',
                               'C3',
                               'P2RY12',
                               'CX3CR1',
                               'FCN1',
                               'VCAN',
                               'THBS1',
                               'TF',
                               'PLP1',
                               'MAG',
                               'MGP',
                               'ISLR',
                               'IGFBP6',
                               'PTPRZ1',
                               'BCAN',
                               'FABP7'), group.by = 'Major.celltype') + 
  RotatedAxis()
colfunc = colorRampPalette(c('lightgrey', 'blue'))(100)
p = p +
  theme(axis.text = element_text(color = 'black', size = 20)) + 
  scale_color_gradientn(colors = colfunc,
                        limits = c(-1.2, 2.5)) 

q = DotPlot(obj1, features = c('FceRIa',
                               'CD101',
                               'CD49f',
                               'CD29',
                               'CD16',
                               'CD352 ',
                               'CD278',
                               'CD14',
                               'CX3CR1',
                               'CD36',
                               'CD57',
                               'CD82',
                               'CD73',
                               'CD224',
                               'GPR56',
                               'CD24'), 
            group.by = 'Major.celltype',
            assay = 'ADT', cols = c('lightgrey', 'maroon4')) + 
  RotatedAxis()
colfunc = colorRampPalette(c('lightgrey', 'maroon4'))(100)
q = q +
  theme(axis.text = element_text(color = 'black', size = 20)) + 
  scale_color_gradientn(colors = colfunc,
                        limits = c(-1.2, 2.5)) 
ggarrange(p,q,ncol = 2, nrow = 1)

dev.off()
