#IOU calculation function
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

library(ComplexHeatmap)
library(circlize)

info = read.csv('/media/user/disk3/singleCellRef/1Sep/finalCells/info2.csv',
                row.names = 1)  #Reference meta data
cells = info[, c('riverCelltype', 'label1')]
cells = unique(cells)
rownames(cells) = cells[,2]
cells = cells[order(rownames(cells)), ]
df = as.data.frame(cells$riverCelltype)
colnames(df) = 'celltype'

meta = read.csv('meta_16Sept.csv', row.names = 1) #Spatial meta data
data = readRDS('binary_15Sept.rds') #Binary cell type data
identical(rownames(meta), rownames(data))

num = grep('GBM', meta$short_histology)
meta = meta[num, ]
data = data[num, ]
identical(rownames(meta), rownames(data))

af = c('all', 'LE_GM', 'LE_WGM', 'LE_WM',
       'CT_control', 'CT', 'PNZ_ccGSC', 'PAN', 'HBV_ccGSC',
       'MVP')

p = list()

for(a in 1:length(af)){
  
  if(af[a] == 'all'){
    binary = data
  } else {
    num = which(meta$AF == af[a])
    binary = data[num, ]
  }
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
  mat = mat[order(rownames(mat)), order(rownames(mat))]
  mat1 = mat
  diag(mat1) = NA
  print('checking!!!!!!!!!!')
  print(identical(colnames(mat1), rownames(cells)))
  
  p[[a]] = Heatmap(mat1, 
                   top_annotation = HeatmapAnnotation(df = df, 
                                                      col = list('celltype' = c('Neuron' = 'aquamarine4',
                                                                                'Oligo_OPC' = 'plum4',
                                                                                'Astrocyte' = 'orange',
                                                                                'Neoplastic' = 'dodgerblue',
                                                                                'Microglia' = 'purple',
                                                                                'Myeloid' = 'hotpink',
                                                                                'Lymphoid' = 'navy',
                                                                                'Fibroblast' = 'burlywood',
                                                                                'Pericyte' = 'maroon4',
                                                                                'Endothelial' = 'chartreuse',
                                                                                'RBC' = 'olivedrab4')),
                                                      show_legend = F,
                                                      show_annotation_name = F),
                   col = colorRamp2(c(0, 10, 40), 
                                    c('orange', 'white', 'purple')),
                   right_annotation = rowAnnotation(df = df,
                                                    col = list('celltype' = c('Neuron' = 'aquamarine4',
                                                                              'Oligo_OPC' = 'plum4',
                                                                              'Astrocyte' = 'orange',
                                                                              'Neoplastic' = 'dodgerblue',
                                                                              'Microglia' = 'purple',
                                                                              'Myeloid' = 'hotpink',
                                                                              'Lymphoid' = 'navy',
                                                                              'Fibroblast' = 'burlywood',
                                                                              'Pericyte' = 'maroon4',
                                                                              'Endothelial' = 'chartreuse',
                                                                              'RBC' = 'olivedrab4')),
                                                    show_legend = F,
                                                    show_annotation_name = F),
                   column_title = paste0('IOU: ', af[a], 'GBM'),
                   show_heatmap_legend = FALSE,
                   column_names_gp = gpar(fontsize = 20),
                   row_names_gp = gpar(fontsize = 20))
}

ht_list = p[[10]] + p[[1]] + p[[2]] + p[[3]] + p[[4]] +
  p[[5]] + p[[7]] + p[[8]] + p[[9]]

pdf('Fig2D.pdf', width = 150, height = 15)
  draw(ht_list, auto_adjust = F)
dev.off()
