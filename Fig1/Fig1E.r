library(matrixStats)
library(ComplexHeatmap)
library(circlize)

obj = readRDS('finalObj_citeseq.rds') #Cite-Seq suerat object
meta = read.csv('meta_13Oct.csv', row.names = 1)  #Cite-seq meta data
identical(rownames(obj@meta.data), rownames(meta))

adt = obj@assays$ADT@data
identical(colnames(adt), rownames(meta))

num = which(is.na(meta$celltypes))
meta = meta[-num, ]
adt = adt[, rownames(meta)]
identical(colnames(adt), rownames(meta))

mat = aggregate(t(adt), list(meta$celltypes), mean)
rownames(mat) = mat[,1]
mat = mat[,-1]

mat = t(mat)
rownames(mat) = gsub('-TotalSeqC', '', rownames(mat))
rownames(mat) = gsub(' ', '', rownames(mat))
mat1 = (mat - rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]

df = read.csv('ab_heatmap_celltypes_anno1.csv', row.names = 1)  #column annotation file
rdf = read.csv('ab_heatmap_celltypes_anno.csv', row.names = 1)  #row annotation file
colnames(mat1) = gsub(' ', '', colnames(mat1))
colnames(mat1) = gsub('stromalcell', 'stromal cell', colnames(mat1))
mat2 = mat1[rownames(df), rownames(rdf)]

df = as.data.frame(df[, 1:2])
col_fun = colorRamp2(c(-0.5, 0, 1), c("orangered", "white", "purple"))

colnames(mat2) = gsub('CD8_6_1_3', 'CD8_TRM', colnames(mat2))

pdf('Fig1E.pdf', width = 23, height = 8)
Heatmap(t(mat2), 
        top_annotation = HeatmapAnnotation(df = df, 
                                         col = list('MajorCelltype' = c('DC' = 'salmon',
                                                                        'endothelial' = 'seagreen2',
                                                                        'granulocyte' = 'goldenrod2', 
                                                                        'lymphocyte' = 'firebrick',
                                                                        'macrophage' = 'hotpink',
                                                                        'microglia' = 'purple',
                                                                        'monocyte' = 'navy',
                                                                        'oligodendrocyte' = 'plum2',
                                                                        'stromal cell' = 'burlywood',
                                                                        'tumor cell' = 'dodgerblue',
                                                                        'misc' = 'grey80'),
                                                    'correlation' = col_fun)),
        cluster_columns = F, cluster_rows = F,
        right_annotation = rowAnnotation(df = rdf,
                                         col = list('MajorCelltype' = c('DC' = 'salmon',
                                                                        'endothelial' = 'seagreen2',
                                                                        'granulocyte' = 'goldenrod2', 
                                                                        'lymphocyte' = 'firebrick',
                                                                        'macrophage' = 'hotpink',
                                                                        'microglia' = 'purple',
                                                                        'monocyte' = 'navy',
                                                                        'oligodendrocyte' = 'plum2',
                                                                        'stromal cell' = 'burlywood',
                                                                        'tumor cell' = 'dodgerblue',
                                                                        'misc' = 'grey80'))))
dev.off()
