library(ComplexHeatmap)
library(circlize)
library(survival)
library(survminer)

data = read.csv('mergedSurvival_NS.csv', row.names = 1) #STable 11

df = as.data.frame(data[, 1:8])
df1 = df[, -5]
mat = as.data.frame(data[, 9:ncol(data)])

info = read.csv('celtypeAnno.csv')  

info = info[order(info[,1]), ]
mat = mat[, order(colnames(mat))]

rdf = as.data.frame(info[, 2:3])
colnames(rdf) = c('AF', 'celltype')
r1 = rowAnnotation(df = as.data.frame(rdf), col = list(celltype = c('Neuron' = 'aquamarine4',
                                                                         'Oligo_OPC' = 'plum4',
                                                                         'Astrocyte' = 'orange',
                                                                         'Neoplastic' = 'dodgerblue',
                                                                         'Microglia' = 'purple',
                                                                         'Myeloid' = 'hotpink',
                                                                         'Lymphoid' = 'navy',
                                                                         'Fibroblast' = 'burlywood',
                                                                         'Pericyte' = 'maroon4',
                                                                         'Endothelial' = 'chartreuse',
                                                                         'RBC' = 'olivedrab4', 
                                                                         'af' = 'white'),
                                                       AF = c('CT' = 'salmon',
                                                              'CT_control' = 'yellow',
                                                              'CT_ccGSC' = 'firebrick', 
                                                              'LE_GM' = 'darkolivegreen4',
                                                              'LE_WM' = 'goldenrod2',
                                                              'MVP' = 'darkmagenta',
                                                              'HBV_ccGSC' = 'magenta',
                                                              'PAN' = 'dodgerblue',
                                                              'PNZ_ccGSC' = 'navy',
                                                              'ct' = 'white')))


h1 = Heatmap(t(mat), column_names_gp = gpar(fontsize = 0),
             row_names_gp = gpar(fontsize = 8),
        col = colorRamp2(c(-2,0,2), c('orange', 'white', 'purple')),
        top_annotation = HeatmapAnnotation(df = as.data.frame(df1),
                                           col = list(Dataset = c('SNUH' = 'salmon2',
                                                                  'CPTAC' = 'purple'),
                                                      vital_status = c('1' = 'black',
                                                                       '0' = 'white'),
                                                      Ecotype = c('good' = 'dodgerblue',
                                                                  'poor' = 'indianred1'),
                                                      Ecotype_MGMT = c('good' = 'dodgerblue',
                                                                       'intermediate' = 'goldenrod2',
                                                                       'poor' = 'indianred1'),
                                                      cluster = c('1' = 'salmon',
                                                                  '2' = 'seagreen2',
                                                                  '3' = 'hotpink',
                                                                  '4' = 'lightblue',
                                                                  '4a' = 'dodgerblue',
                                                                  '5' = 'goldenrod2',
                                                                  '6' = 'navy',
                                                                  '7' = 'orange'),
                                                      MGMT = c('methylated' = 'seagreen2',
                                                               'unmethylated' = 'plum2'),
                                                      OS = colorRamp2(c(0, 27, 55),
                                                                      c('yellow', 'orange', 'red')))),
        cluster_columns = F,
        right_annotation = r1)

pdf('Fig8A_combHeatmapSSGSEA.pdf', width = 15, height = 10)
  h1
dev.off()


fit = survfit(Surv(OS, vital_status) ~ Ecotype_MGMT, data = df)
pdf('Fig8C_ecotypeMGMT1.pdf', width = 8, height = 7)
  ggsurvplot(fit, data = df, risk.table = T, pval = T, 
           palette = c('dodgerblue',
                       'goldenrod2',
                       'indianred1'))
dev.off()

df1 = subset(df, df$MGMT == 'methylated')
fit1 = survfit(Surv(OS, vital_status) ~ Ecotype, data = df1)
df2 = subset(df, df$MGMT == 'unmethylated')
fit2 = survfit(Surv(OS, vital_status) ~ Ecotype, data = df2)

pdf('Fig8B_ecotypeMGMT1.pdf', width = 10, height = 7)
p1 = ggsurvplot(fit1, data = df1, risk.table = T, pval = T, 
           palette = c('dodgerblue',
                       'indianred1'))

p2 = ggsurvplot(fit2, data = df2, risk.table = T, pval = T, 
                palette = c('dodgerblue',
                            'indianred1'))

arrange_ggsurvplots(x = list(p1, p2), ncol = 2, nrow = 1)
dev.off()
