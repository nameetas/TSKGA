data = readRDS('log2cpm.rds') #normalized single-cell data
meta = read.csv('meta_28Sept.csv', row.names = 1) #single-cell meta data

identical(colnames(data), rownames(meta))

num = which(meta$X28Sept %in% c('Treg', 
                                'Tcell_prolif',
                                'CD8', 'CD4',
                                'CD8_TRM', 'NKT',
                                'Bcell_plasma', 'Bcell'))
meta = meta[num, ]
data = data[, rownames(meta)]

identical(colnames(data), rownames(meta))

mat = aggregate(t(data), list(meta$X28Sept), mean)
rownames(mat) = mat[,1]
mat = mat[,-1]
mat = t(mat)

genes = c('B2M','SLAMF6','TNFSF14','FASLG','CRTAM','GZMK','KLRC1',
          'KLRB1','CD226','KLRF1','KLRD1','GZMB','FAS','SIRPG','CD28',
          'TIGIT','SLAMF1','CEACAM1','CTLA4','CLEC2D','FLT3LG','CD2','CD96',
          'CD6','CD40LG','SLAMF7','TNFRSF17','SDC1','CD38','ITGB7',
          'TNFRSF13B','LY9','SELL')
mat = mat[genes, ]

library(matrixStats)
mat1 = (mat - rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]

library(ComplexHeatmap)
library(circlize)

Heatmap(mat1, col = colorRamp2(c(-2, 0, 2), c('orangered', 'white', 'purple')))
