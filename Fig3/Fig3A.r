library(ggplot2)
library(ggpubr)
library(Seurat)

meta = read.csv('meta_28Sept.csv', row.names = 1) #single-cell meta data
cells = unique(meta$Major.celltype)
cells = cells[-6] #remove cells with no Major.celltype assigned

num = which(meta$groups == '')
meta = meta[-num, ]

pdf('Fig3A.pdf', width = 40, height = 25)

my_comparisons = list(c('control', 'IDHmut_G3'), 
                       c('control', 'GBMcore'), 
                       c('control', 'GBMperi'),
                       c('IDHmut_G3', 'GBMcore'),
                       c('IDHmut_G3', 'GBMperi'),
                       c('GBMperi', 'GBMcore'))

for(i in 1:length(cells)){
  num = which(meta$Major.celltype == cells[i])
  df = as.data.frame(table(meta[num, 'groups'],
                           meta[num, 'X28Sept']))
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
  
  print(ggarrange(plotlist = plist, nrow = 2, 
                  ncol = ceiling(length(subcells)/2)))
}
dev.off()
