library(Seurat)
library(ggplot2)
library(ggpubr)


visobj = readRDS('visiumObj_20Sept.rds')  #spatial seurat object
mat = readRDS('spots_pairs.rds') #matrix showing gene expression data per spot. Subset interactions provided as RDS object

pairs = c('CD86:CD28',
          'CD86:CTLA4',
          'CADM1:CRTAM')
num = which(colnames(mat) %in% pairs)

mat1 = mat[, num]
mat1 = mat1[, 1:3]
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

#select specific slides
s1 = p1[[7]]
s2 = p1[[28]]

img6[[1]] = ggarrange(s1, s2, ncol = 1, nrow = 2, widths = c(1, 1))

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
  cols = c('red', 'green', 'yellow')
  names(cols) = c(g1, g2, 'both')
  
  p1 = list()
  for(j in 1:32){
    p1[[j]] = p[[j]] + scale_fill_manual(values = cols)
  }
  
  png(paste0(colnames(mat1)[i], '.png'), width = 8000, height = 4000)
  print(ggarrange(plotlist = p1, ncol = 8, nrow = 4))
  dev.off()
  
  s1 = p1[[7]]
  s2 = p1[[28]]
  
  img6[[i + 1]] = ggarrange(s1, s2, ncol = 1, nrow = 2, widths = c(1, 1))
}


binary = readRDS('binary_15Sept.rds')
identical(rownames(binary), rownames(visobj@meta.data))

visobj@active.ident = factor(binary[, 'Treg'])
names(visobj@active.ident) = rownames(binary)
visobj@meta.data$Class = binary[, 'Treg']
visobj1 = subset(visobj, subset = Class == 1)

p = SpatialDimPlot(visobj1, pt.size.factor = 1.4,
                   alpha = 0.9, crop = F)

for(j in 1:length(names(visobj1@images))){
  p[[j]] = p[[j]] + scale_fill_manual(values = c('1' = 'red4')) +
    theme(plot.title = element_text(size=50))
}

img6[[length(img6) + 1]] = ggarrange(p[[7]], p[[28]], ncol = 1, nrow = 2)

colnames(binary) = gsub('CD8_6_1_3', 'CD8_TRM', colnames(binary))
visobj@active.ident = factor(binary[, 'CD8_TRM'])
names(visobj@active.ident) = rownames(binary)
visobj@meta.data$Class = binary[, 'CD8_TRM']
visobj1 = subset(visobj, subset = Class == 1)

p = SpatialDimPlot(visobj1, pt.size.factor = 1.4,
                   alpha = 0.9, crop = F)

for(j in 1:length(names(visobj1@images))){
  p[[j]] = p[[j]] + scale_fill_manual(values = c('1' = 'red4')) +
    theme(plot.title = element_text(size=50))
}

img6[[length(img6) + 1]] = ggarrange(p[[7]], p[[28]], ncol = 1, nrow = 2)

length(img6)

#Figure 3D
ggarrange(plotlist = img6, ncol = length(img6), nrow = 1)
