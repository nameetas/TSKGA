data = read.csv('/media/user/disk31/visiumNewSubmission/CNV_AFreplacementordered.csv')

library(reshape2)
library(Seurat)
library(ggplot2)

df = melt(data, id.vars = 'X')
df$X = factor(df$X, levels = data$X)
ggplot(df, aes(x = X, y = variable, color = value)) + geom_point() + 
  RotatedAxis() + xlab('IvyGAP assigned AF') + ylab('CNV value assigned AF')
