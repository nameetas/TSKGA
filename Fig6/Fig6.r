library(Seurat)
library(ggplot2)
library(ggpubr)
library(randomcoloR)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)
library(EnhancedVolcano)

query = readRDS('vasc_27Jul_GBM.harmony.rds') #Vascular data seurat obj
meta = read.csv('queryMetaData_21Sept_riverplot.csv', row.names = 1)  #vascular meta data
meta = meta[rownames(query@meta.data), ]
identical(rownames(meta), rownames(query@meta.data))

query@meta.data = meta

palette = distinctColorPalette(20)

query = subset(query, subset = label != 'misc')
p1 = DimPlot(query, group.by = 'label', cols = palette)
img1 = p1

counts = as.matrix(query@assays$RNA@counts)
cpm = apply(counts, 2, function(x)(x/sum(x))*1000000)
cpm = log2(cpm + 1)

genes = read.csv('fig6b_genes.csv')
cpm1 = cpm[genes$genes, ]

mat = aggregate(t(cpm1), list(query@meta.data$label), mean)
rownames(mat) = mat[,1]
mat = mat[,-1]
mat = t(mat)

df = as.data.frame(genes$label)
colnames(df) = 'celltype'
zs = (mat - rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]

h1 = Heatmap(t(zs), cluster_columns = F, 
             col = colorRamp2(c(-2,0,2), c('orange', 'white', 'purple')),
             top_annotation = HeatmapAnnotation(df = df, 
                                                col = list(celltype = c('Endo' = '#CEE47B',
                                                                        'FBMC' = '#EE82EE',
                                                                        'Fib' = '#000080', 
                                                                        'Peri' = '#DA9371',
                                                                        'SMC' = '#EEB422', 
                                                                        'Proliferating' = '#D85C7D'))))
h1
img2 = grid.grabExpr(draw(h1))

deg = read.csv('Endothelial_filteredMetaData_catgeryDEG.csv', row.names = 1)  #STable 8
num = which(is.na(deg$pval))
deg = deg[-num, ]
deg = as.data.frame(deg)

pc = read.csv("/media/user/disk2/references/x1", header = F)
sel = intersect(pc[,1], rownames(deg))
deg = deg[sel, ] 

img3 = EnhancedVolcano(deg,
                       lab = rownames(deg),
                       x = 'log2FC',
                       y = 'pval',
                       labSize = 10,
                       selectLab = c('CTSB',
                                     'CTSD',
                                     'FCGRT',
                                     'FCGBP',
                                     'CXCL8',
                                     'GDI2',
                                     'PLVAP'))

ref = readRDS('/media/user/disk3/singleCellRef/1Sep/finalCells/log2cpm.rds')  #Normalized reference data
refInfo = read.csv('/media/user/disk3/singleCellRef/1Sep/finalCells/info2.csv', row.names = 1)  #Reference meta data

labels = c('Art1', 'Art2', 'Art3', 'Endo_cap', 'gbmEndo_9_1', 
           'gbmEndo_9_4', 'Venous', 'Venule')

num = which(refInfo$label1 %in% labels)
refI = refInfo[num, ]
refCPM = ref[, rownames(refI)]

genes = read.csv('Fig6E.csv') #genes of interest

refCPM = refCPM[genes$Gene, ]

mat = aggregate(t(refCPM), list(refI$label1), mean)
rownames(mat) = mat[,1]
mat = mat[,-1]
mat = t(mat)

df = as.data.frame(genes$Class)
colnames(df) = 'Class'
zs = (mat - rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]

df2 = as.data.frame(colnames(zs))
colnames(df2) = 'celltype'

h1 = Heatmap(t(zs), 
        top_annotation = HeatmapAnnotation(df = df,
                                           col = list(Class = c('actin remodeling and cytoskeleton' = 'salmon',
                                                                'collagen and associated factors' = 'seagreen2',
                                                                'ECM interaction' = 'hotpink',
                                                                'immature cell markers' = 'yellow2',
                                                                'matricellular proteins' = 'dodgerblue',
                                                                'scavenging capillary ' = 'purple',
                                                                'scavenging capillary_R ' = 'orange',
                                                                'tip markers' = 'firebrick1'))), 
        cluster_rows = F, 
        col = colorRamp2(c(-2,0,2), c('orange', 'white', 'purple')), 
        cluster_columns = F,
        right_annotation = rowAnnotation(df = df2, 
                                         col = list(celltype = c('Art1' = 'plum2',
                                                                 'Art2' = 'khaki2',
                                                                 'Art3' = 'greenyellow',
                                                                 'Endo_cap' = 'dodgerblue',
                                                                 'gbmEndo_9_1' = 'mediumorchid1',
                                                                 'gbmEndo_9_4' = 'darkmagenta',
                                                                 'Venous' = 'olivedrab4',
                                                                 'Venule' = 'red'))))
img4 = grid.grabExpr(draw(h1))

endo = readRDS('JuneData/Endothelial_filtered.RDS') #Reference vascular object
endoMeta = read.csv('JuneData/Endothelial_filteredMetaData.csv', row.names = 1) #Reference vascular meta data
identical(rownames(endo@meta.data), rownames(endoMeta))

endo@meta.data = endoMeta

endo1 = subset(endo, subset = merged.cell.type %in% c('Art1', 'Art2', 'Art3',
                                                      'Endo_cap', 'gbmEndo_9_1',
                                                      'gbmEndo_9_4', 'Venous', 'Venule'))
p1 = DimPlot(endo1, group.by = 'merged.cell.type', cols = c('plum2', 'khaki2', 'greenyellow', 'dodgerblue',
                                                            'mediumorchid1', 'darkmagenta', 'olivedrab4', 'red'))
endo1@meta.data$category = gsub('AVM_control', 'Normal', endo1@meta.data$category)
endo1@meta.data$category = factor(endo1@meta.data$category, levels = c('Normal', 'GBM'))
p2 = DimPlot(endo1, group.by = 'merged.cell.type', split.by = 'category',
             ncol = 1,
             cols = c('plum2', 'khaki2', 'greenyellow', 'dodgerblue',
                                                                'mediumorchid1', 'darkmagenta', 'olivedrab4', 'red'))
img5 = ggarrange(p1, p2, ncol = 2, nrow = 1)

query = readRDS('vasc_27Jul_GBM.harmony.rds')
qInfo = read.csv('queryMetaData_21Sept_riverplot.csv', row.names = 1)
qInfo = qInfo[rownames(query@meta.data), ]

endoRef = readRDS('vascObj/ref.rds')
rInfo = read.csv('vascObj/refMeta.csv', row.names = 1)
rInfo = rInfo[rownames(endoRef@meta.data), ]

identical(rownames(query@meta.data), rownames(qInfo))
identical(rownames(endoRef@meta.data), rownames(rInfo))

query@meta.data = qInfo
endoRef@meta.data = rInfo

query1 = subset(query, subset = clusters %in% c(4,5,6))
endoRef1 = subset(endoRef, subset = celltype.l1 %in% c('Fib', 'FBMC'))

c1 = as.matrix(query1@assays$RNA@counts)
qcpm = apply(c1, 2, function(x)(x/sum(x))*1000000)
qcpm = log2(qcpm + 1)

c2 = as.matrix(endoRef1@assays$RNA@counts)
rcpm = apply(c2, 2, function(x)(x/sum(x))*1000000)
rcpm = log2(rcpm + 1)

genes = read.csv('fig6F.csv')

gIntersect = intersect(rownames(qcpm), rownames(rcpm))

qcpm = qcpm[gIntersect, ]
rcpm = rcpm[gIntersect, ]


rmat = aggregate(t(rcpm), list(endoRef1@meta.data$celltype.l1), mean)
rownames(rmat) = rmat[,1]
rmat = rmat[,-1]
rmat = t(rmat)

qmat = aggregate(t(qcpm), list(query1@meta.data$clusters), mean)
rownames(qmat) = qmat[,1]
qmat = qmat[,-1]
qmat = t(qmat)

comb = cbind.data.frame(rmat[genes$genes,], qmat[genes$genes, ])
