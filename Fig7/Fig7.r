library(Seurat)
library(ggplot2)
library(ggpubr)
library(inlmisc)

obj = readRDS('visiumObj_20Sept.rds') #Spatial seurat object
binary = readRDS('binary_15Sept.rds') #Binary data

obj1 = obj

cols = colnames(binary)

totalpixels = nrow(obj@meta.data)

thres = read.csv('celltypeThreshold_15Sept.csv', row.names = 1) #Threshold values per celltype

obj1@active.ident = factor(obj1@meta.data$AF)
names(obj1@active.ident) = rownames(obj1@meta.data)

p = SpatialDimPlot(obj1, pt.size.factor = 1.4, alpha = 0.9, stroke = 0)


img1[[1]] = p[[22]] + scale_fill_manual(values = c('CT' = 'salmon',
                                              'CT_control' = 'yellow',
                                              'CT_ccGSC' = 'firebrick', 
                                              'LE_GM' = 'darkolivegreen4',
                                              'LE_WGM' = 'chartreuse',
                                              'LE_WM' = 'goldenrod2',
                                              'MVP' = 'darkmagenta',
                                              'HBV_ccGSC' = 'magenta',
                                              'PAN' = 'dodgerblue',
                                              'PNZ_ccGSC' = 'navy')) + 
  NoLegend() + theme(plot.title = element_text(size = 30))

cols = c('gbmEndo_9_1', 'gbmEndo_9_3', 'gbmEndo_9_4', 'gbmEndoPeri_2',
         'gbmFib_4', 'gbmFib_5', 'gbmFib_6')

plist = list()
plist[[1]] = p[[22]] + scale_fill_manual(values = c('CT' = 'salmon',
                                                    'CT_control' = 'yellow',
                                                    'CT_ccGSC' = 'firebrick', 
                                                    'LE_GM' = 'darkolivegreen4',
                                                    'LE_WGM' = 'chartreuse',
                                                    'LE_WM' = 'goldenrod2',
                                                    'MVP' = 'darkmagenta',
                                                    'HBV_ccGSC' = 'magenta',
                                                    'PAN' = 'dodgerblue',
                                                    'PNZ_ccGSC' = 'navy')) + 
  NoLegend() + theme(plot.title = element_text(size = 30))

for(i in 1:length(cols)){
  j = i + 1
  obj1@meta.data$Class = obj1@meta.data[, cols[i]]
  obj2 = subset(obj1, subset = Class == 1)
  obj2@active.ident = factor(obj2@meta.data[, cols[i]])
  names(obj2@active.ident) = rownames(obj2@meta.data)
  
  q = SpatialDimPlot(obj2, pt.size.factor = 1.4, alpha = 1,
                     image.alpha = 0.4, crop = F) + NoLegend()
    
  plist[[j]] = q[[22]] + scale_fill_manual(values = c('0' = adjustcolor('white', alpha.f = 0),
                                                      '1' = 'firebrick')) +
    theme(plot.title = element_text(size = 30)) + NoLegend() +
    ggtitle(cols[i])
  
}

pdf('Fig7A_gbmBinary.pdf', width = 40, height = 20)
  ggarrange(plotlist = plist, ncol = 4, nrow = 2)
dev.off()


data = readRDS('/media/user/disk3/singleCellRef/1Sep/finalCells/log2cpm.rds') #Normalized reference data
meta = read.csv('/media/user/disk3/singleCellRef/1Sep/finalCells/info2.csv', row.names = 1) #Reference meta data

identical(colnames(data), rownames(meta))

genes = c('JAG1',
          'CD46',
          'NOTCH1',
          'NOTCH2',
          'NOTCH3',
          'NOTCH4',
          'VASN',
          'THY1',
          'ITGAM',
          'ITGAV',
          'ITGAX',
          'ITGB2',
          'ITGB3',
          'FYN')

data = data[genes, ]

mat = aggregate(t(data), list(meta$label1), mean)
rownames(mat) = mat[,1]
mat = mat[,-1]
mat = t(mat)

library(matrixStats)
library(ComplexHeatmap)
library(circlize)

mat1 = (mat - rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]

df1 = unique(meta[, c('label1', 'riverCelltype')])
df1 = df1[order(df1$label1), ]

mat1 = mat1[, order(colnames(mat1))]
identical(df1$label1, colnames(mat1))
df1 = as.data.frame(df1$riverCelltype)
colnames(df1) = 'celltype'

ha = HeatmapAnnotation(df = df1,
                       col = list(celltype = c('Neuron' = 'aquamarine4',
                                               'Oligo_OPC' = 'plum4',
                                               'Astrocyte' ='orange',
                                               'Neoplastic' = 'dodgerblue',
                                               'Microglia' = 'purple',
                                               'Myeloid' = 'hotpink',
                                               'Lymphoid' = 'navy',
                                               'Fibroblast' = 'burlywood',
                                               'Pericyte' = 'maroon4',
                                               'Endothelial' = 'chartreuse',
                                               'RBC' = 'olivedrab')))

scAvg = read.csv('/media/user/disk2/completeAnalysis/visium_15Sept/LRpairs/scAVG_outertrack1.csv', row.names = 1)
sc = scAvg[genes, 1:9]

identical(rownames(sc), rownames(mat1))

h1 = Heatmap(mat1, 
        top_annotation = ha,
        right_annotation = rowAnnotation(df = as.data.frame(sc),
                                         col = list(LE_GM = colorRamp2(c(0, 2), c('white', 'darkolivegreen4')),
                                                    LE_WGM = colorRamp2(c(0, 2), c('white', 'chartreuse')),
                                                    LE_WM = colorRamp2(c(0, 2), c('white', 'goldenrod2')),
                                                    CT_control = colorRamp2(c(0, 2), c('white', 'yellow')),
                                                    CT = colorRamp2(c(0, 2), c('white', 'salmon')),
                                                    PNZ_ccGSC = colorRamp2(c(0, 2), c('white', 'navy')),
                                                    PAN = colorRamp2(c(0, 2), c('white', 'dodgerblue')),
                                                    HBV_ccGSC = colorRamp2(c(0, 2), c('white', 'magenta')),
                                                    MVP = colorRamp2(c(0, 2), c('white', 'darkmagenta')))))
pdf('Fig7B_refLabel1_heatmap.pdf', width = 15, height = 6)
  draw(h1)
dev.off()
