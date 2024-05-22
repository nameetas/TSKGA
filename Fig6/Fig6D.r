library(Seurat)
library(ggplot2)
library(ggpubr)
library(Nebulosa)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(patchwork)
library(monocle)


IOU = function(x,y) {
  U = sum(x | y)
  if(U < 100){
    iou = 0
  } else {
    I = sum(x & y)
    iou = (I/U) * 100
  }
  return(iou)
}


refinfo = read.csv('/media/user/disk31/singleCellRef/1Sep/finalCells/info2.csv',
                   row.names = 1)

vismeta = read.csv('visium_15Sept/meta_16Sept.csv', row.names = 1)
data = readRDS('visium_15Sept/binary_15Sept.rds')
identical(rownames(vismeta), rownames(data))

num = grep('GBM', vismeta$short_histology)

vismeta = vismeta[num, ]
data = data[num, ]
identical(rownames(vismeta), rownames(data))

labels = c('Oligo_2_1', 'Oligo_2_2', 'Oligo_2_3_1', 'Oligo_2_3_2')

binary = data

mat = matrix(NA, ncol = ncol(binary), nrow = ncol(binary))
rownames(mat) = colnames(binary)
colnames(mat) = colnames(binary)

for(j in 1:ncol(binary)){
  vec1 = binary[, j]
  for(k in 1:ncol(binary)){
    vec2 = binary[, k]
    mat[colnames(binary)[j], colnames(binary)[k]] = IOU(vec1, vec2)
  }
}
mat1 = mat[, labels]

cells = unique(refinfo[, c('riverCelltype', 'label1')])
cells = as.data.frame(cells)
cells = cells[order(cells$label1),]

mat1 = mat1[order(rownames(mat1)), ]
rownames(mat1) = gsub('CD8_6_1_3', 'CD8_TRM', rownames(mat1))
rownames(mat1) = gsub(' ', '', rownames(mat1))
identical(rownames(mat1), cells$label1)

df = as.data.frame(cells[, 1])
colnames(df) = 'celltype'

heatmapMat = mat1

num = which(rownames(mat1) %in% c('RBC', 'OPC', 'TC_MTC'))
heatmapMat = mat1[-num, ]

df = as.data.frame(df[-num, ])
colnames(df) = 'celltype'

h1 = Heatmap(t(heatmapMat), 
             col = colorRamp2(c(0, 10, 40), 
                              c('orange', 'white', 'purple')),
             column_names_gp = gpar(fontsize = 25),
             top_annotation = HeatmapAnnotation(df = df, 
                                                col = list(celltype = c('Astrocyte' = '#FFA500',
                                                                        'Endothelial' = '#7FFF00',
                                                                        'Fibroblast' = '#DEB887',
                                                                        'Lymphoid' = '#000080',
                                                                        'Microglia' = '#A020F0',
                                                                        'Myeloid' = '#FF69B4',
                                                                        'Neoplastic' = '#1E90FF',
                                                                        'Neuron' = '#458B74',
                                                                        'Oligodendrocyte' = '#8B668B',
                                                                        'Pericyte' = '#8B1C62',
                                                                        'RBC' = '#698B22'))),
             row_names_gp = gpar(fontsize = 25), name = 'Exp zscore')
