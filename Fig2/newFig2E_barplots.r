setwd('/media/user/disk21/completeAnalysis/')

meta = read.table('cellbrowser/scRNA-Seq/meta1.tsv', row.names = 1, header = T, sep = '\t')

num = which(meta$tissue_histology %in% c('core:GBM', 'peri:GBM'))
m1 = meta[num, ]

meta$orig.ident = gsub('SNU34citeseq', 'SNU34', meta$orig.ident)

num = which(m1$orig.ident %in% c('core:SNU27', 'core:SNU32', 'core:SNU33', 'core:SNU34',
                                 'core:SNU35', 'core:SNU46', 'core:SNU51',
                                 'peri:SNU27', 'peri:SNU32', 'peri:SNU33', 'peri:SNU34',
                                 'peri:SNU35', 'peri:SNU46', 'peri:SNU51'))
m2 = m1[num, ]

res = as.data.frame(table(m2$orig.ident, m2$Major.celltype, m2$celltype))
colnames(res) = c('origIdent', 'majorCelltype', 'celltype', 'freq')
res$origIdent = as.character(res$origIdent)
res$tissue = sapply(res$origIdent, function(x) strsplit(x, ":")[[1]][1], USE.NAMES=FALSE)
res$sample = sapply(res$ori, function(x) strsplit(x, ":")[[1]][2], USE.NAMES=FALSE)

res$celltype = gsub(' ', '', res$celltype)

res$major = NA

num = which(res$celltype %in% c('Bcell', 'Bcell_plasma', 'CD4', 'CD8', 'CD8_TRM', 'Tcell_prolif', 
                            'Treg', 'NKT'))
res[num, 'major'] = 'lymphocytes'

num = which(res$celltype %in% c('DC', 'Granulocyte', 'Mac_5_1_2', 'Mac_5_1_3', 'Mac_prolif', 'Mg_1_1', 
                            'Mg_1_2', 'Mg_1_3', 'Mg_4_1', 'Mg_prolif', 'Monocyte'))
res[num, 'major'] = 'myeloid'

num = which(res$celltype %in% c('TC_AC', 'TC_mesh', 'TC_mesnh', 'TC_NPC', 'TC_oligo', 'TC_OPC', 'TC_prolif'))
res[num, 'major'] = 'neoplastic'


num = which(res$celltype %in% c('Oligo_2_1', 'Oligo_2_2', 'Oligo_2_3_1', 
                                'Oligo_2_3_2'))
res[num, 'major'] = 'oligodendrocytes'


num = which(is.na(res$major))
res1 = res[-num, ]

cols1 = read.csv('celltypeColors_11Dec.csv')

cols = cols1$colors
names(cols) = cols1$names

library(ggplot2)
library(Seurat)

# pdf('scRNAseq_barplots.pdf', width = 6, height = 10)
# ggplot(res, aes(x = sample, y = Freq, fill = Var3)) + 
#   geom_bar(position = 'fill', stat = 'identity') + facet_wrap(.~ Var2 + tissue, ncol = 2) + 
#   scale_fill_manual(values = cols) + theme_classic() +
#   RotatedAxis()
# dev.off()

res$tissue = gsub('core', 'Core', res$tissue)
res$tissue = gsub('peri', 'Peri', res$tissue)

prop2 = as.data.frame(table(res$tissue, res$major, res$sample))
prop2$id = paste0(prop2$Var1, '_', prop2$Var2, '_', prop2$Var3)
prop2 = prop2[, c('Freq', 'id')]
colnames(prop2)[1] = 'prop2'

res$id = paste0(res$tissue, '_', res$major, '_', res$sample)

sumPatient = as.data.frame(table(m2$orig.ident))
colnames(sumPatient) = c('origIdent', 'sumPatient')

df1 = subset(res, res$major == 'myeloid')
df1 = merge(df1, sumPatient, by = 'origIdent', all.x = T)
df1$prop = df1$freq/df1$sumPatient
df1 = merge(df1, prop2, by = 'id', all.x = T)
df1$finalProp = df1$freq/df1$prop2

library(ggpubr)

p1 = ggplot(df1, aes(x = sample, y = prop, fill = celltype)) + 
  geom_bar(position = 'stack', stat = 'identity') + facet_wrap(.~ tissue, ncol = 2) + 
  scale_fill_manual(values = cols) + theme_classic() +
  RotatedAxis() + NoLegend() + xlab('') + ylab('Proportions') + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        strip.background = element_blank(),
        strip.text.x = element_text(family = 'sans', 
                                    colour = 'black', 
                                    size = 6),
        axis.title = element_text(family = 'sans', 
                                  colour = 'black', 
                                  size = 5),
        plot.margin = unit(c(0, 0, 0, 0), 'lines'))


df2 = subset(res, res$major == 'lymphocytes')
df2 = merge(df2, sumPatient, by = 'origIdent', all.x = T)
df2$prop = df2$freq/df2$sumPatient
df2 = merge(df2, prop2, by = 'id', all.x = T)
df2$finalProp = df2$freq/df2$prop2

p2 = ggplot(df2, aes(x = sample, y = prop, fill = celltype)) + 
  geom_bar(position = 'stack', stat = 'identity') + facet_wrap(.~ tissue , ncol = 2) + 
  scale_fill_manual(values = cols) + theme_classic() +
  RotatedAxis() + NoLegend() + xlab('') + ylab('Proportions') + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        strip.background = element_blank(),
        strip.text.x = element_text(family = 'sans', 
                                    colour = 'black', 
                                    size = 6),
        axis.title = element_text(family = 'sans', 
                                  colour = 'black', 
                                  size = 5),
        plot.margin = unit(c(0, 0, 0, 0), 'lines'))

df3 = subset(res, res$major == 'neoplastic')
df3 = merge(df3, sumPatient, by = 'origIdent', all.x = T)
df3$prop = df3$freq/df3$sumPatient
df3 = merge(df3, prop2, by = 'id', all.x = T)
df3$finalProp = df3$freq/df3$prop2

p3 = ggplot(df3, aes(x = sample, y = prop, fill = celltype)) + 
  geom_bar(position = 'stack', stat = 'identity') + facet_wrap(.~ tissue, ncol = 2) + 
  scale_fill_manual(values = cols) + theme_classic() +
  RotatedAxis() + NoLegend() + xlab('') + ylab('Proportions') + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        strip.background = element_blank(),
        strip.text.x = element_text(family = 'sans', 
                                    colour = 'black', 
                                    size = 6),
        axis.title = element_text(family = 'sans', 
                                  colour = 'black', 
                                  size = 5),
        plot.margin = unit(c(0, 0, 0, 0), 'lines'))


q1 = ggplot(df1, aes(x = sample, y = prop, fill = celltype)) + 
  geom_bar(position = 'fill', stat = 'identity') + facet_wrap(.~ tissue, ncol = 2) + 
  scale_fill_manual(values = cols) + theme_classic() +
  RotatedAxis() + NoLegend() + xlab('') + ylab('Relative proportions') + 
  theme(axis.text.x = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        axis.text.y = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_text(family = 'sans', 
                                  colour = 'black', 
                                  size = 5),
        plot.margin = unit(c(0, 0, 0, 0), 'lines'))

q2 = ggplot(df2, aes(x = sample, y = prop, fill = celltype)) + 
  geom_bar(position = 'fill', stat = 'identity') + facet_wrap(.~ tissue, ncol = 2) + 
  scale_fill_manual(values = cols) + theme_classic() +
  RotatedAxis() + NoLegend() + xlab('') + ylab('Relative proportions') + 
  theme(axis.text.x = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        axis.text.y = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_text(family = 'sans', 
                                  colour = 'black', 
                                  size = 5),
        plot.margin = unit(c(0, 0, 0, 0), 'lines'))

q3 = ggplot(df3, aes(x = sample, y = prop, fill = celltype)) + 
  geom_bar(position = 'fill', stat = 'identity') + facet_wrap(.~ tissue, ncol = 2) + 
  scale_fill_manual(values = cols) + theme_classic() +
  RotatedAxis() + NoLegend() + xlab('') + ylab('Relative proportions') + 
  theme(axis.text.x = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        axis.text.y = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_text(family = 'sans', 
                                  colour = 'black', 
                                  size = 5),
        plot.margin = unit(c(0, 0, 0, 0), 'lines'))


df4 = subset(res, res$major == 'oligodendrocytes')
df4 = merge(df4, sumPatient, by = 'origIdent', all.x = T)
df4$prop = df4$freq/df4$sumPatient
df4 = merge(df4, prop2, by = 'id', all.x = T)
df4$finalProp = df4$freq/df4$prop2

library(ggpubr)

p4 = ggplot(df4, aes(x = sample, y = prop, fill = celltype)) + 
  geom_bar(position = 'stack', stat = 'identity') + facet_wrap(.~ tissue , ncol = 2) + 
  scale_fill_manual(values = cols) + theme_classic() +
  RotatedAxis() + NoLegend() + xlab('') + ylab('Proportions') + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        strip.background = element_blank(),
        strip.text.x = element_text(family = 'sans', 
                                    colour = 'black', 
                                    size = 6),
        axis.title = element_text(family = 'sans', 
                                  colour = 'black', 
                                  size = 5),
        plot.margin = unit(c(0, 0, 0, 0), 'lines'))

q4 = ggplot(df4, aes(x = sample, y = prop, fill = celltype)) + 
  geom_bar(position = 'fill', stat = 'identity') + facet_wrap(.~ tissue, ncol = 2) + 
  scale_fill_manual(values = cols) + theme_classic() +
  RotatedAxis() + NoLegend() + xlab('') + ylab('Relative proportions') + 
  theme(axis.text.x = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        axis.text.y = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_text(family = 'sans', 
                                  colour = 'black', 
                                  size = 5),
        plot.margin = unit(c(0, 0, 0, 0), 'lines'))



# 
# 
# 
# pdf('newFig2E_barplots.pdf', width = 10, height = 5)
# ggarrange(p1, p2, p3, q1, q2, q3, ncol = 3, nrow = 2)
# dev.off()
# 
# 
# 
# 
# 
# 
# cells = unique(df1$celltype)
# propLympho = matrix(NA, nrow = length(cells), ncol = 4)
# propLympho = as.data.frame(propLympho)
# colnames(propLympho) = c('core', 'peri', 'pvalue', 'pvalue2')
# for(i in 1:length(cells)){
#   num = which(df1$celltype == cells[i])
#   data = df1[num, ]
#   data = subset(data, data$freq > 0)
#   
#   out = t.test(prop ~ tissue, data = data)
#   propLympho[i, 'pvalue1'] = out$p.value
#   
#   out = t.test(finalProp ~ tissue, data = data)
#   propLympho[i, 'pvalue2'] = out$p.value
# }
# propLympho$celltype = cells
# write.csv(propLympho, file = 'propmyeloid_ttest.csv')
# 
# 
# cells = unique(df2$celltype)
# propLympho = matrix(NA, nrow = length(cells), ncol = 4)
# propLympho = as.data.frame(propLympho)
# colnames(propLympho) = c('core', 'peri', 'pvalue', 'pvalue2')
# for(i in 1:length(cells)){
#   num = which(df2$celltype == cells[i])
#   data = df2[num, ]
#   data = subset(data, data$freq > 0)
#   
#   out = t.test(prop ~ tissue, data = data)
#   propLympho[i, 'pvalue1'] = out$p.value
#   
#   out = t.test(finalProp ~ tissue, data = data)
#   propLympho[i, 'pvalue2'] = out$p.value
# }
# propLympho$celltype = cells
# write.csv(propLympho, file = 'propLympho_ttest.csv')
# 
# 
# # cells = unique(df3$celltype)
# # propLympho = matrix(NA, nrow = length(cells), ncol = 4)
# # propLympho = as.data.frame(propLympho)
# # colnames(propLympho) = c('core', 'peri', 'pvalue', 'pvalue2')
# # for(i in 1:length(cells)){
# #   if
# #   num = which(df3$celltype == cells[i])
# #   data = df3[num, ]
# #   data = subset(data, data$freq > 0)
# #   
# #   out = t.test(prop ~ tissue, data = data)
# #   propLympho[i, 'pvalue1'] = out$p.value
# #   
# #   out = t.test(finalProp ~ tissue, data = data)
# #   propLympho[i, 'pvalue2'] = out$p.value
# # }
# # propLympho$celltype = cells
# # write.csv(propLympho, file = 'propneoplastic_ttest.csv')



pdf('Fig2C-J.pdf', width = 2.15, height = 7.87)
  ggarrange(p3, q3, p4, q4, p1, q1, p2, q2,
            ncol = 1, nrow = 8)
dev.off()
