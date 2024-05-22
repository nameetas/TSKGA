setwd('/media/user/disk21/completeAnalysis/')

data = read.csv('scRef_xeniumPanel.csv', row.names = 1)

df = data[, 1:3]

rowdf = data[, 1:2]

coldf = as.data.frame(c('Neurons', 'Neurons',
          'Oligodendrocyte/OPC', 'Oligodendrocyte/OPC', 'Oligodendrocyte/OPC', 'Oligodendrocyte/OPC', 'Oligodendrocyte/OPC',
          'Astrocytes', 'Astrocytes', 'Astrocytes', 
          'Neoplastic', 'Neoplastic', 'Neoplastic', 'Neoplastic', 'Neoplastic', 'Neoplastic', 'Neoplastic',
          'Microglia', 'Microglia', 'Microglia',
          'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid',
          'Lymphoid', 'Lymphoid', 'Lymphoid', 'Lymphoid', 'Lymphoid', 'Lymphoid', 'Lymphoid', 'Lymphoid',
          'Fibroblasts', 'Fibroblasts', 'Fibroblasts', 'Fibroblasts', 'Fibroblasts',
          'Pericytes', 'Pericytes', 'Pericytes', 'Pericytes', 'Pericytes',
          'Endothelial', 'Endothelial', 'Endothelial', 'Endothelial', 'Endothelial', 'Endothelial', 'Endothelial', 'Endothelial', 'Endothelial', 'Endothelial'))
colnames(coldf) = 'majorCelltype'
coldf$celltype = colnames(data)[-c(1:3)]

mat = as.matrix(data[, -c(1:3)])

library(ComplexHeatmap)
library(circlize)

cols1 = read.csv('celltypeColors_11Dec.csv')

cols = cols1$colors
names(cols) = cols1$names

idx = 1:nrow(df)
num = which(df$label == '')
idx = idx[-num]

labs = unique(df$label)
labs = labs[-1]

ha1 = rowAnnotation(df = rowdf, 
                    col = list(celltype = cols,
                               major_celltype = c('Neurons' = '#458B73',
                                                  'Oligodendrocyte/OPC' = '#89638A',
                                                  'Astrocytes' = '#FFA300',
                                                  'Neoplastic' = '#1B8FFE',
                                                  'Microglia' = '#9E1EEF',
                                                  'Myeloid' = '#FF69B3',
                                                  'Lymphoid' = '#00007E',
                                                  'Fibroblasts' = '#DEB985',
                                                  'Pericytes' = '#8B1961',
                                                  'Endothelial' = '#7DFF00')),
                    mark = anno_mark(at = idx, 
                              labels = labs, 
                              which = "row", 
                              labels_gp = gpar(fontsize = 5, 
                                               fontfamily = 'sans')),
                    annotation_name_gp = gpar(fontsize = 5, fontfamily = 'sans'))

rowdf$major_celltype = factor(rowdf$major_celltype, levels = unique(rowdf$major_celltype))

ha2 = HeatmapAnnotation(df = coldf, 
                        col = list(majorCelltype = c('Neurons' = '#458B73',
                                                     'Oligodendrocyte/OPC' = '#89638A',
                                                     'Astrocytes' = '#FFA300',
                                                     'Neoplastic' = '#1B8FFE',
                                                     'Microglia' = '#9E1EEF',
                                                     'Myeloid' = '#FF69B3',
                                                     'Lymphoid' = '#00007E',
                                                     'Fibroblasts' = '#DEB985',
                                                     'Pericytes' = '#8B1961',
                                                     'Endothelial' = '#7DFF00'),
                                   celltype = cols),
                        annotation_name_gp = gpar(fontsize = 5, fontfamily = 'sans'))

coldf$majorCelltype = factor(coldf$majorCelltype, levels = c('Neurons', 'Oligodendrocyte/OPC',
                                                             'Astrocytes', 'Neoplastic',
                                                             'Microglia', 'Myeloid',
                                                             'Lymphoid', 'Fibroblasts',
                                                             'Pericytes', 'Endothelial'))
# 
# pdf('scRef_xeniumPanel.pdf', width = 8, height = 30)
# Heatmap(mat, col = colorRamp2(c(-2, 0, 2), c('orange', 'white', 'purple')),
#         cluster_rows = F, cluster_columns = F, 
#         right_annotation = ha1, top_annotation = ha2, 
#         column_split = coldf$majorCelltype, column_title = NULL,
#         row_split = rowdf$major_celltype, row_title = NULL,
#         row_names_gp = gpar(fontsize = 0))
# dev.off()
# num = which(is.na(df$label))
# df[num, 'label'] = ' '
# 
# text_list = split(df$label, df$major_celltype)
# 
# text_list = list(text_list$Neurons, text_list$Astrocytes, text_list$`Oligodendrocyte/OPC`,
#                  text_list$Neoplastic, text_list$Microglia, text_list$Myeloid,
#                  text_list$Lymphoid, text_list$Fibroblasts, text_list$Pericytes,
#                  text_list$Endothelial)
# 
# names(text_list) = unique(rowdf$major_celltype)

majorCols = c('Neurons' = '#458B73',
              'Oligodendrocyte/OPC' = '#89638A',
              'Astrocytes' = '#FFA300',
              'Neoplastic' = '#1B8FFE',
              'Microglia' = '#9E1EEF',
              'Myeloid' = '#FF69B3',
              'Lymphoid' = '#00007E',
              'Fibroblasts' = '#DEB985',
              'Pericytes' = '#8B1961',
              'Endothelial' = '#7DFF00')

pdf('scRef_xeniumPanel.pdf', width = 3, height = 9)
h1 = Heatmap(mat, col = colorRamp2(c(-2, 0, 2), c('orange', 'white', 'purple')),
        cluster_rows = F, cluster_columns = F, 
        right_annotation = ha1, top_annotation = ha2, 
        column_split = coldf$majorCelltype, column_title = NULL,
        row_split = rowdf$major_celltype, row_title = NULL,
        row_names_gp = gpar(fontsize = 0), name = 'Log2CPM',
        column_names_gp = gpar(fontsize = 0))

draw(h1, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
# for(i in 1:10) {
#   decorate_annotation("major_celltype", slice = i, {
#     grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(0.5, "mm"),
#               gp=gpar(fontsize = 8),
#               just = "left", draw = T)
#   })
# }

for(i in 1:10){
  for(j in 1:10){
    if(i == j){
      decorate_heatmap_body("Log2CPM", row_slice = i, column_slice = j, {
        
        grid.rect(gp = gpar(fill = 'transparent', col = "black", lwd = 2))
        
      })
    }
  }
}
dev.off()


ha = Heatmap(t(heatmapMat), 
             col = colorRamp2(c(0, 10, 40), 
                              c('orange', 'white', 'purple')),
             column_names_gp = gpar(fontsize = 6, fontfamily = 'sans'),
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
                                                                        'RBC' = '#698B22')),
                                                annotation_name_gp = gpar(fontsize = 5, fontfamily = 'sans')),
             row_names_gp = gpar(fontsize = 6, fontfamily = 'sans'), name = 'Exp zscore')

pdf('Fig6C.pdf', width = 6, height = 2)
draw(ha, show_annotation_legend = FALSE)
dev.off()
