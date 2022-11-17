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


scObj = readRDS('scRNA_25Aug/finalObj_13Sept.rds')  #single-cell seurat object
scMeta = read.csv('scRNA_25Aug/meta_28Sept.csv', row.names = 1) #single-cell meta data

identical(rownames(scObj@meta.data), rownames(scMeta))

scObj@meta.data = scMeta

Tumap = read.csv('scRNA_25Aug/RCC1/Tumap.csv', row.names = 1) #UMAP coordinates calculated by RCC

identical(rownames(Tumap), rownames(scObj@reductions$umap@cell.embeddings))

scObj@reductions$umap@cell.embeddings[,1] = Tumap$Lev3_UMAP1
scObj@reductions$umap@cell.embeddings[,2] = Tumap$Lev3_UMAP2

obj1 = subset(scObj, subset = Major.celltype == 'oligodendrocyte')
unique(obj1@meta.data$tissue)

obj2 = subset(obj1, subset = tissue != 'other')

obj2@meta.data$tissue = factor(obj2@meta.data$tissue, 
                               levels = c('control', 'peri', 'core'))

p1 = DimPlot(obj2, group.by = 'X28Sept') + ggtitle(' ')
p2 = DimPlot(obj2, group.by = 'X28Sept', split.by = 'tissue', ncol = 1) + 
  ggtitle(' ')

img1 = ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = T)

meta = obj2@meta.data
num = which(meta$groups == '')
meta = meta[-num, ]
df = as.data.frame(table(meta[, 'groups'],
                         meta[, 'X28Sept']))
colnames(df) = c('tissue', 'celltype', 'Freq')
df$cells = nrow(meta)
df$percent = df$Freq/df$cells
df$tissue = as.character(df$tissue)
df$label = sapply(df$tissue, 
                  function(x) strsplit(x, ":")[[1]][1], 
                  USE.NAMES=FALSE)
df$label = gsub('core', 'GBMcore', df$label)
df$label = gsub('peri', 'GBMperi', df$label)

df$label = factor(df$label, levels = c('control', 'IDHmut_G3',
                                       'GBMperi', 'GBMcore'))
subcells = unique(df$celltype)

my_comparisons = list(c('control', 'IDHmut_G3'), 
                      c('control', 'GBMcore'), 
                      c('control', 'GBMperi'),
                      c('IDHmut_G3', 'GBMcore'),
                      c('IDHmut_G3', 'GBMperi'),
                      c('GBMperi', 'GBMcore'))

plist = list()
for(j in 1:length(subcells)){
  num1 = which(df$celltype == subcells[j])
  df1 = df[num1, ]
  
  plist[[j]] = ggplot(df1, aes(x = label, y = percent, fill = label)) + 
    geom_boxplot() + theme_classic() +
    stat_compare_means(method = 't.test', size = 15,
                       comparisons = my_comparisons) +
    NoLegend() + RotatedAxis() + xlab('') +
    theme(axis.text = element_text(color = 'black',
                                   size = 30),
          strip.text = element_text(size = 35))
  plist[[j]] = annotate_figure(plist[[j]], 
                               top = text_grob(subcells[j], 
                                               color = "black", 
                                               size = 35))
}

img2 = ggarrange(plotlist = plist, ncol = 2, nrow = 2)

deg = read.csv('scRNA_25Aug/finalAnnot_1Sep/oligodendrocyteDEG_GBM.csv', 
               row.names = 1) #Differential gene expression STable 6
num = which(is.na(deg$pvalue))
deg = deg[-num, ]
deg = as.data.frame(deg)

img3 = EnhancedVolcano(deg,
                     lab = rownames(deg),
                     x = 'log2fc',
                     y = 'pvalue',
                     selectLab = c('OPALIN', 'HLA-A', 'LGALS1', 
                                   'HLA-E', 'MIF', 'IFITM3',
                                   'LGALS3BP', 'CLU'),
                     labSize = 10)

o1 = plot_density(obj2, c('NGFR'), size = 0.5) + 
  theme(plot.title = element_text(size = 30))
o2 = plot_density(obj2, c('FGFR2'), size = 0.5) + 
  theme(plot.title = element_text(size = 30))
o3 = plot_density(obj2, c('S100A1'), size = 0.5) + 
  theme(plot.title = element_text(size = 30))
o4 = plot_density(obj2, c('LGALS1'), size = 0.5) + 
  theme(plot.title = element_text(size = 30))
o5 = plot_density(obj2, c('OPALIN'), size = 0.5) + 
  theme(plot.title = element_text(size = 30))

img4 = ggarrange(o1, o2, o3, o4, o5, ncol = 5, nrow = 1)

refinfo = read.csv('/media/user/disk3/singleCellRef/1Sep/finalCells/info2.csv',
                row.names = 1)  #Reference meta data

vismeta = read.csv('visium_15Sept/meta_16Sept.csv', row.names = 1)  #Spatial meta data
data = readRDS('visium_15Sept/binary_15Sept.rds') #Binary data per spot
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
identical(rownames(mat1), cells$label1)

df = as.data.frame(cells[, 1])
colnames(df) = 'celltype'

heatmapMat = mat1

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
                                                                   'Oligo_OPC' = '#8B668B',
                                                                   'Pericyte' = '#8B1C62',
                                                                   'RBC' = '#698B22'))),
        row_names_gp = gpar(fontsize = 25))

img5 = grid.grabExpr(draw(h1))

visobj = readRDS('/media/user/disk2/completeAnalysis/visium_15Sept/visiumObj_20Sept.rds') #Spatial seurat object
mat = readRDS('visium_15Sept/LRpairs_spatial/spots_pairs.rds')  #Ligand receptor pairs

pairs = c('FGF11:FGFR2', 'OMG:NGFR', 'S100A1:TLR4')
num = which(colnames(mat) %in% pairs)

mat1 = mat[, num]
identical(rownames(mat1), rownames(visobj@meta.data))

img6 = list()

visobj@active.ident = factor(visobj@meta.data$AF, levels = c('LE_GM', 
                                                             'LE_WGM',
                                                             'LE_WM',
                                                             'CT_control',
                                                             'CT',
                                                             'PNZ_ccGSC',
                                                             'PAN',
                                                             'HBV_ccGSC', 
                                                             'MVP'))
names(visobj@active.ident) = rownames(visobj@meta.data)

p = SpatialDimPlot(visobj, pt.size.factor = 1.5, alpha = 0.9) + NoLegend()

p1 = list()
for(j in 1:32){
  p1[[j]] = p[[j]] + scale_fill_manual(values = c('CT' = 'salmon',
                                         'CT_control' = 'yellow',
                                         'CT_ccGSC' = 'firebrick', 
                                         'LE_GM' = 'darkolivegreen4',
                                         'LE_WGM' = 'chartreuse',
                                         'LE_WM' = 'goldenrod2',
                                         'MVP' = 'darkmagenta',
                                         'HBV_ccGSC' = 'magenta',
                                         'PAN' = 'dodgerblue',
                                         'PNZ_ccGSC' = 'navy')) + 
    NoLegend() + 
    theme(plot.title = element_text(size = 35))
}

s1 = p1[[4]]
s2 = p1[[18]]

img6[[1]] = ggarrange(s1, s2, ncol = 2, nrow = 1, widths = c(1, 1))

for(i in 1:ncol(mat1)){
  print(i)
  
  g1 = strsplit(colnames(mat1)[i], ':')[[1]][1]
  g2 = strsplit(colnames(mat1)[i], ':')[[1]][2]
  
  if(g1 == g2){
    next
  }
  
  if(length(unique(mat1[,i])) == 1){
    next
  }
  
  visobj@active.ident = factor(mat1[, i], levels = c('none', g1, g2, 'both'))
  names(visobj@active.ident) = rownames(mat1)
  visobj@meta.data$Class = mat1[,i]
  
  visobj1 = subset(visobj, subset = Class != 'none')
  
  p = SpatialDimPlot(visobj1, pt.size.factor = 1.5,
                     alpha = 0.7, stroke = 0) + NoLegend()
  
  p1 = list()
  for(j in 1:32){
    p1[[j]] = p[[j]] + scale_fill_manual(name = c(g1, g2, 'both'), 
                                         values = c('red', 'green', 'yellow')) + 
      NoLegend() + ggtitle('')
  }
  
  s1 = p1[[4]]
  s2 = p1[[18]]
  
  img6[[i + 1]] = ggarrange(s1, s2, ncol = 2, nrow = 1, widths = c(1, 1))
}

img7 = ggplot() + theme_void()

pdf('Fig5.pdf', width = 40, height = 50)
  print(ggarrange(img1, img2, 
                  img3, 
                  ggarrange(img4, img5, ncol = 1, nrow = 2),
                  ggarrange(plotlist = img6,
                                            ncol = 1, nrow = 4),
                  img7, 
                  ncol = 2, nrow = 3,
                  widths = c(0.8, 1.2),
                  heights = c(1.2, 0.8, 1.2)))
dev.off()
