info = read.csv('meta_16Sept.csv', row.names = 1) #visium meta data
df = as.data.frame(table(info$sample, info$AF))

library(ggplot2)
library(ggpubr)

colnames(df) = c('slide', 'AF', 'value')
df$slide = factor(df$slide, levels = c('SNU38',
                                       'SNU19',
                                       'SNU22',
                                       'SNU40',
                                       'SNU27',
                                       'SNU17',
                                       'SNU24',
                                       'SNU25',
                                       'SNU21',
                                       'SNU23',
                                       'SNU43',
                                       'SNU18',
                                       'SNU51',
                                       'SNU46',
                                       'SNU16',
                                       'SNU34',
                                       'SNU33'))

meta = read.csv('infoAnno2.csv', row.names = 1) #sample info file
meta = meta[levels(df$slide), ]
#df$slide = factor(df$slide, levels = meta[,1])

identical(rownames(meta), levels(df[,1]))

df$AF = factor(df$AF, levels = c('LE_GM', 'LE_WGM', 'LE_WM', 'CT_control',
                                 'CT', 'PNZ_ccGSC', 'PAN', 'HBV_ccGSC',
                                 'MVP'))

p1 = ggplot(df, aes(x = slide, y = value, fill = AF)) +
  geom_bar(position = 'fill', stat = 'identity') + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c('LE_GM' = 'darkolivegreen4',
                               'LE_WGM' = 'chartreuse',
                               'LE_WM' = 'goldenrod2',
                               'CT_control' = 'yellow',
                               'CT' = 'salmon',
                               'PNZ_ccGSC' = 'navy',
                               'PAN' = 'dodgerblue',
                               'HBV_ccGSC' = 'magenta',
                               'MVP' = 'darkmagenta'))

mat = matrix(1:17, nrow = 2, ncol = nrow(meta))
library(ComplexHeatmap)
library(circlize)

p2 = Heatmap(mat, 
             top_annotation = HeatmapAnnotation(df = as.data.frame(meta[, -c(1, 2)]),
                                                col = list('slides' = colorRamp2(c(0, 3, 6), c('white', 'deeppink2', 'deeppink4')),
                                                           'Gender' = c('Female' = 'plum2',
                                                                        'Male' = 'lightblue'),
                                                           'Age' = colorRamp2(c(27, 55, 82), c('darkslategray1', 'darkslategray3', 'darkslategrey')),
                                                           'Diagnosis' = c('Central neurocytoma' = 'yellow',
                                                                           'Odg, IDHm,1p/19q-codel, G3' = 'violet',
                                                                           'Meningioma' = 'pink', 'PLNTY, G1' = 'dodgerblue',
                                                                           'Medulloblastoma, G4' = 'lightgreen',
                                                                           'Ast, IDHm, G4' = 'deeppink4', 
                                                                           'Glioblastoma' = 'goldenrod2', 
                                                                           'Ast, IDHm, G3' = 'darkslategray4',
                                                                           'Ast, IDHwt, G3' = 'indianred1'),
                                                           'IDH1_mutation' = c('mutant' = 'springgreen2', 'wildtype' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'IDH2_mutation' = c('wildtype' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'X1p_19q' = c('co-deletion' = 'salmon', 'no-deletion' = 'seagreen',
                                                                         '19q only deleted' = 'violet'),
                                                           'EGFR_amplification' = c('absent' = adjustcolor('burlywood4', alpha.f = 0.6), 'amplification' = 'salmon2'),
                                                           'X9p21_CDKN2A_CDKN2B_deletion' = c('absent' = adjustcolor('burlywood4', alpha.f = 0.6),
                                                                                              'hemizygous deletion' = 'blueviolet',
                                                                                              'homozygous deletion' = 'firebrick4'),
                                                           'X10q23_PTEN_locus_deletion' = c('absent' = adjustcolor('burlywood4', alpha.f = 0.6),
                                                                                            'hemizygous deletion' = 'blueviolet'),
                                                           'ATRX_mutation' = c('mutant' = 'springgreen2', 'wildtype' = adjustcolor('burlywood4', alpha.f = 0.6),
                                                                               'VUS' = 'plum2'),
                                                           'TP53_mutation' = c('mutant' = 'springgreen2', 'wildtype' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'TERTp_mutation' = c('mutant' = 'springgreen2', 'wildtype' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'BRAF_mutation' = c('mutant' = 'springgreen2', 'wildtype' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'BRAF_rearrangement' = c('absent' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'H3F3A_mutation' = c('wildtype' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'HIST1H3B_mutation' = c('wildtype' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'FGFR1_mutation' = c('wildtype' = adjustcolor('burlywood4', alpha.f = 0.6), 'mutant' = 'springgreen2'),
                                                           'PDGFRA_amplification' = c('absent' = adjustcolor('burlywood4', alpha.f = 0.6), 'amplification' = 'salmon2'),
                                                           'RB1_mutation' = c('wildtype' = adjustcolor('burlywood4', alpha.f = 0.6), 'mutant' = 'springgreen2'),
                                                           'RB1_deletion' = c('absent' = adjustcolor('burlywood4', alpha.f = 0.6), 'hemizygous deletion' = 'blueviolet'),
                                                           'NF1_mutation' = c('mutant' = 'springgreen2', 'wildtype' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'NF1_deletion' = c('absent' = adjustcolor('burlywood4', alpha.f = 0.6), 'hemizygous deletion' = 'blueviolet'),
                                                           'NF2_mutation' = c('wildtype' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'NF2_deletion' = c('hemizygous deletion' = 'blueviolet', 'absent' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'MYCN_amplification' = c('absent' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'PTEN_mutation' = c('mutant' = 'springgreen2', 'absent' = adjustcolor('burlywood4', alpha.f = 0.6)),
                                                           'Ki67.' = colorRamp2(c(5.4, 48, 79.8), c('darkslategray1', 'darkslategray3', 'darkslategrey')),
                                                           'MGMT_methylation_status' = c('missing' = 'lightgrey', 'methylated' = 'darkolivegreen1',
                                                                                         'unmethylated' = adjustcolor('hotpink4', alpha.f = 0.7))
                                                ), border = T, gp = gpar(col = "black"),
                                                na_col = 'gray88'),
             cluster_columns = F, col = colorRamp2(c(0, 15, 30), c(adjustcolor('orange', alpha.f = 0),
                                                                   adjustcolor('white', alpha.f = 0),
                                                                   adjustcolor('purple', alpha.f = 0))),
             cluster_rows = F)

p3 = grid.grabExpr(draw(p2))

df1 = read.csv('scFreq.csv')  #sample cell type frequency distribution
df1$sample = factor(df1$sample, levels = c('SNU38',
                                           'SNU19',
                                           'SNU22',
                                           'SNU40',
                                           'SNU27',
                                           'SNU17',
                                           'SNU24',
                                           'SNU25',
                                           'SNU21',
                                           'SNU23',
                                           'SNU43',
                                           'SNU18',
                                           'SNU51',
                                           'SNU46',
                                           'SNU16',
                                           'SNU34',
                                           'SNU33'))
p4 = ggplot(df1, aes(x = sample, y = value, fill = celltype)) +
  geom_bar(position = 'fill', stat = 'identity') + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(color = 'black')) +
  scale_fill_manual(values = c('DC' = 'salmon',
                               'endothelial' = 'seagreen2',
                               'granulocyte' = 'goldenrod2', 
                               'lymphocyte' = 'firebrick',
                               'macrophage' = 'hotpink',
                               'microglia' = 'purple',
                               'monocyte' = 'navy',
                               'oligodendrocyte' = 'plum2',
                               'stromal cell' = 'burlywood',
                               'tumor cell' = 'dodgerblue'))

meta1 = read.csv('scInfoAnno.csv', row.names = 1) #scRNA-seq patient data

mat = matrix(1:17, nrow = 2, ncol = nrow(meta1))
library(ComplexHeatmap)
library(circlize)

p5 = Heatmap(mat, 
             top_annotation = HeatmapAnnotation(df = as.data.frame(meta1[, -c(1)]),
                                                col = list('control' = c('1' = 'black'),
                                                           'peri' = c('1' = 'black'),
                                                           'core' = c('1' = 'black')), na_col = 'white', 
                                                gp = gpar(col = "black", size = 1),
                                                border = F),
             cluster_columns = F, col = colorRamp2(c(0, 15, 30), c(adjustcolor('orange', alpha.f = 0),
                                                                   adjustcolor('white', alpha.f = 0),
                                                                   adjustcolor('purple', alpha.f = 0))),
             cluster_rows = F)

p6 = grid.grabExpr(draw(p5))


pdf('fig1A.pdf', width = 9, height = 20)
ggarrange(p1, p3, p4, p6, ncol = 1, nrow = 4)
dev.off()
#User will have to adjust the width and size of plots in image editor.
