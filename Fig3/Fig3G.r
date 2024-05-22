data = read.csv('/media/user/disk21/completeAnalysis/visium_2Jun/AFanno_coordinates_CE_counts_Other.csv')

library(ggplot2)
library(Seurat)

ggplot(data, aes(x = AF, y = CE_counts)) + geom_boxplot() + RotatedAxis()
