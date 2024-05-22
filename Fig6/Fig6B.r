setwd('/media/user/disk21/completeAnalysis/')

data = read.csv('Oligo_2_1vsOligo_2_3_2xenium.csv')
num = which(data$select == 1)
data = data[num, ]

library(ggplot2)
data$log2FCsc = as.numeric(data$log2FCsc)
rownames(data) = data$gene

data$labels = ''
sel =  c('OPALIN', 'MOBP', 'APOD', 
         'KIF19', 'HLA-A', 'HLA-DRA',
         'GSN', 'SERPINA3', 'NGFR',
         'TIMP1', 'TUBB2B',
         'S100A1', 'CLU')
data[sel, 'labels'] = c('OPALIN', 'MOBP', 'APOD', 
                        'KIF19', 'HLA-A', 'HLA-DRA', 
                        'GSN', 'SERPINA3', 'NGFR',
                        'TIMP1', 'TUBB2B',
                        'S100A1', 'CLU')

pdf('Oligo_2_1vsOligo_2_3_2xenium_scatterplot.pdf', width = 20, height = 20)
ggplot(data, aes(x = log2FCsc, y = log2fcsqr, label = labels)) + geom_point() + 
  geom_hline(yintercept = 0, col = 'burlywood4', linetype = 'dashed') + 
  geom_vline(xintercept = 0, col = 'burlywood4', linetype = 'dashed') + 
  geom_text(nudge_x = 0, nudge_y = 0.1,
            check_overlap = F, size = 7, family = 'sans') + theme_classic() + 
  xlab('Log2FC in scRNA-seq') + ylab('Log2FC in Xenium') + 
  ggtitle('Oligo_2_1 vs Oligo_2_3_2 DEG') + 
  theme(text = element_text(size = 6, family = 'sans'),
        title = element_text(size = 7, family = 'sans'))
dev.off()
