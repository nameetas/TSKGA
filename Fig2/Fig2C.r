library(matrixStats)
library(ggridges)
library(reshape2)

avg = readRDS('l2cpmavg_15Sept.rds')  #log2cpm avg cell type data
#To generate this object, use gene sets for cell types and average it across all the spots in spatial data
colnames(avg) = gsub('CD8_6_1_3', 'CD8_TRM', colnames(avg))
meta = read.csv('meta_16Sept.csv', row.names = 1) #spatial meta data
identical(rownames(meta), rownames(avg))

num1 = which(meta$source == 'Frozen')
num2 = which(meta$source == 'FFPE')
mat = avg

avg1 = mat[num1, ]
avg1 = t(avg1)
m1 = (avg1 - rowMeans(avg1, na.rm = T))/(rowSds(as.matrix(avg1), na.rm = T))[row(avg1)]
m1 = t(m1)

avg2 = mat[num2, ]
avg2 = t(avg2)
m2 = (avg2 - rowMeans(avg2, na.rm = T))/(rowSds(as.matrix(avg2), na.rm = T))[row(avg2)]
m2 = t(m2)

mat1 = rbind.data.frame(m1, m2)
mat1 = mat1[rownames(meta), ]

mat2 = cbind.data.frame(meta$AF, mat1)
colnames(mat2)[1] = 'AF' 

df = melt(mat2, id.vars = 'AF')
colnames(df) = c('AF', 'celltype', 'value')

anno = read.csv('../survival_BulkData/celtypeAnno.csv') #cell type annotation data
colnames(anno)[1] = 'celltype'

df1 = merge(df, anno, by = 'celltype', all.x = T)

cells = unique(df1$celltype)

pdf('Fig2C.pdf', width = 200, height = 120)
p1 = list()
for(i in 1:length(cells)){
  num = which(df1$celltype == cells[i])
  df2 = df1[num, ]
  
  
  df2$celltype = 'label'
  p1[[i]] = ggplot(df2, aes(x = value, y = celltype, fill = AF, color = AF)) + 
    geom_density_ridges(alpha = 0.2) +
    theme_ridges() + theme(text = element_text(size = 80)) +
    theme(axis.text = element_text(size = 50)) +
    ggtitle(cells[i]) + 
    theme(legend.position = "none") + xlim(c(-2.5, 5)) + 
    theme(panel.grid.major = element_line(color = "black",
                                          size = 0.15,
                                          linetype = 2)) +
    theme(plot.title = element_text(size = 80)) + 
    scale_fill_manual(values = c('CT' = 'salmon',
                                 'CT_control' = 'yellow',
                                 'LE_GM' = 'darkolivegreen4',
                                 'LE_WGM' = 'chartreuse',
                                 'LE_WM' = 'goldenrod2',
                                 'MVP' = 'darkmagenta',
                                 'HBV_ccGSC' = 'magenta',
                                 'PAN' = 'dodgerblue',
                                 'PNZ_ccGSC' = 'navy')) + 
    scale_color_manual(values = c('CT' = 'salmon',
                                  'CT_control' = 'yellow',
                                  'LE_GM' = 'darkolivegreen4',
                                  'LE_WGM' = 'chartreuse',
                                  'LE_WM' = 'goldenrod2',
                                  'MVP' = 'darkmagenta',
                                  'HBV_ccGSC' = 'magenta',
                                  'PAN' = 'dodgerblue',
                                  'PNZ_ccGSC' = 'navy'))
}

print(ggarrange(plotlist = p1, ncol = 10, nrow = 6))

dev.off()
