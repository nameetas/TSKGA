#Trajectory analysis
library(monocle)
library(Seurat)

oligo_1Sep_GBM_control.harmony <- readRDS("/media/user/disk21/completeAnalysis/scRNA_25Aug/oligo_1Sep.rds")
head(oligo_1Sep_GBM_control.harmony)

unique(oligo_1Sep_GBM_control.harmony$tissue)
table(oligo_1Sep_GBM_control.harmony$tissue)

Idents(oligo_1Sep_GBM_control.harmony) <- oligo_1Sep_GBM_control.harmony$tissue
oligo_1Sep_GBM_control.harmony.filtered <- subset(oligo_1Sep_GBM_control.harmony,ident = "control", invert = TRUE)
unique(oligo_1Sep_GBM_control.harmony.filtered$tissue)


subset <- oligo_1Sep_GBM_control.harmony.filtered

meta = read.table('cellbrowser/scRNA-Seq/meta1.tsv', row.names = 1, header = T,
                  sep = '\t')
meta = meta[rownames(subset@meta.data), ]
subset@meta.data = meta

subset = subset(subset,orig.ident %in% c('core:SNU27', 'core:SNU32', 'core:SNU33', 
                                         'core:SNU34citeseq',
                                         'core:SNU35', 'core:SNU46', 'core:SNU51',
                                         'peri:SNU27', 'peri:SNU32', 'peri:SNU33', 
                                         'peri:SNU34citeseq',
                                         'peri:SNU35', 'peri:SNU46', 'peri:SNU51', 
                                         'NormalBrain:SNU52', 'NormalBrain:SNU57'), invert = FALSE)

unique(subset@meta.data$orig.ident)

#subset = FindVariableFeatures(subset, nfeatures = 2000)
data <- as(as.matrix(subset@assays$RNA@data), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", data = subset@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds<- newCellDataSet(data,
                             phenoData = pd,
                             featureData = fd,
                             lowerDetectionLimit = 0.5,
                             expressionFamily = uninormal()) #because already normalized, scaled in Seurat before


var_genes <- subset[["RNA"]]@var.features
ordering_genes <- var_genes
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)

# monocle_cds <- reduceDimension(monocle_cds, norm_method = "none", 
#                                max_components = 5, scaling = TRUE, 
#                                reduction_method = "DDRTree", pseudo_expr = 0)
# monocle_cds <- orderCells(monocle_cds, root_state = 3)
# saveRDS(monocle_cds, "~/Dropbox/GBM_comnbined/Publication/oligoSep1_monocle_cds_2.RDS")
# 
# #saveRDS(monocle_cds, "~/Dropbox/GBM_comnbined/Publication/oligoSep1_monocle_cds.RDS")
# 
# plot_cell_trajectory(monocle_cds, color_by = "cell.subtype", )+ facet_wrap(~tissue)
# 
# my_genes <- c("MOBP","OPALIN", "HLA-A", "HLA-B","HLA-E","LGALS1")
# cds_subset <- monocle_cds[my_genes,]
# plot_genes_in_pseudotime(cds_subset, color_by = "cell.subtype")


 ##  remote duck :  run 2022 nov 18 1107
monocle_cds <- reduceDimension(monocle_cds, norm_method = "none", max_components = 2, 
                               scaling = TRUE, reduction_method = "DDRTree", pseudo_expr = 0)
monocle_cds <- orderCells(monocle_cds, root_state = 1)
# saveRDS(monocle_cds, "~/Dropbox/GBM_comnbined/Publication/oligoSep1_monocle_cds_2.RDS")
# ##  remote duck :  end

cols1 = read.csv('/media/user/disk21/completeAnalysis/celltypeColors_11Dec.csv')
cols = cols1$colors
names(cols) = cols1$names

### re-level phenodata (tissue): start
library(dplyr)
monocle_cds@phenoData@data %>% head()
monocle_cds@phenoData@data$tissue %>% head()


library(ggplot2)

monocle_cds@phenoData@data$tissue = gsub('core', 'Core', monocle_cds@phenoData@data$tissue)
monocle_cds@phenoData@data$tissue = gsub('peri', 'Peri', monocle_cds@phenoData@data$tissue)
monocle_cds@phenoData@data$tissue_factor <- factor(monocle_cds@phenoData@data$tissue, 
                                                   levels = c("NormalBrain", "Peri", "Core"))


pdf('Fig6E_2Apr24.pdf', width = 6.77, height = 2.5)
plot_cell_trajectory(monocle_cds, color_by = "celltype", cell_size = 1)+ facet_wrap(~ tissue_factor) + 
  scale_color_manual(values = cols) + 
  theme(axis.text.x = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
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
        legend.text = element_text(size = 5, family = 'sans', color = 'black'), 
        legend.title = element_text(size = 6, family = 'sans', color = 'black'))

plot_cell_trajectory(monocle_cds, color_by = "celltype", cell_size = 0.5, 
                     state_number_size = 0.1)+ facet_wrap(~ tissue_factor) + 
  scale_color_manual(values = cols) + 
  theme(axis.text.x = element_text(family = 'sans', 
                                   colour = 'black', 
                                   size = 5),
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
        legend.text = element_text(size = 5, family = 'sans', color = 'black'), 
        legend.title = element_text(size = 6, family = 'sans', color = 'black')) + NoLegend()

dev.off()

genes = c('OPALIN', 'MOBP', 'APOD', 
  'KIF19', 'HLA_A', 'HLA_DRA', 
  'GSN', 'SERPINA3', 'NGFR',
  'TIMP1', 'TUBB2B',
  'S100A1', 'CLU')
monocle_cds@phenoData@data$OPALIN = as.numeric(data['OPALIN', ])
monocle_cds@phenoData@data$MOBP = as.numeric(data['MOBP', ])
monocle_cds@phenoData@data$APOD = as.numeric(data['APOD', ])
monocle_cds@phenoData@data$KIF19 = as.numeric(data['KIF19', ])
monocle_cds@phenoData@data$HLA_A = as.numeric(data['HLA-A', ])
monocle_cds@phenoData@data$HLA_DRA = as.numeric(data['HLA-DRA', ])
monocle_cds@phenoData@data$GSN = as.numeric(data['GSN', ])
monocle_cds@phenoData@data$SERPINA3 = as.numeric(data['SERPINA3', ])
monocle_cds@phenoData@data$NGFR = as.numeric(data['NGFR', ])
monocle_cds@phenoData@data$TIMP1 = as.numeric(data['TIMP1', ])
monocle_cds@phenoData@data$TUBB2B = as.numeric(data['TUBB2B', ])
monocle_cds@phenoData@data$S100A1 = as.numeric(data['S100A1', ])
monocle_cds@phenoData@data$CLU = as.numeric(data['CLU', ])

monocle_cds@phenoData@data$OPALIN = as.numeric(data['OPALIN', ])


rng = colorRampPalette(c('grey90', 'purple'))
rng(10)

p1 = list()

for(i in 1:length(genes)){
  p1[[i]] = plot_cell_trajectory(monocle_cds, color_by = genes[i], cell_size = 0.5)+ 
    facet_wrap(~ tissue_factor) + 
    scale_color_gradientn(colors = rng(10)) + 
    theme(axis.text.x = element_text(family = 'sans', 
                                     colour = 'black', 
                                     size = 5),
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
          legend.text = element_text(size = 5, family = 'sans', color = 'black'), 
          legend.title = element_text(size = 6, family = 'sans', color = 'black'), 
          legend.key.height = unit(5, 'mm'), #change legend key height
          legend.key.width = unit(10, 'mm'))
}

library(ggpubr)

pdf('monocleSupplFig_oligoGenes.pdf', width = 6.77, height = 15)
  ggarrange(plotlist = p1, ncol = 2, nrow = 7)
dev.off()
### re-level phenodata (tissue): end
