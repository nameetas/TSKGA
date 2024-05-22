cols = read.csv('/media/user/disk21/completeAnalysis/celltypeColors_11Dec.csv')
cols1 = as.vector(cols$colors)
names(cols1) = cols$names

setwd('/media/user/disk21/completeAnalysis/visium_2Jun/')

library(ggplot2)

mat = read.csv('/media/user/disk21/completeAnalysis/visium_2Jun/SNUHAUC_AFvsCelltype_27Jul.csv', row.names = 1)
mat = as.data.frame(mat)

library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)

data = mat

# num = which(rownames(data) == 'TC_MTC')
# data = data[-num, ]

data_long <- data %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)
colnames(data_long) <- c("target", "source", "value")
#cells = c('Bcell', 'Bcell_plasma', 'CD4', 'CD8', 'CD8_TRM', 'Treg', 'Tcell_prolif', 'NKT')

# num = which(data_long$source %in% cells)
# data_long = data_long[num, ]

cells = rownames(data)

afs = c('LE_GM', 'LE_WM', 'IT', 'CT', 'PAN2',
        'PAN1', 'PAN', 'PNZ', 'PNZ1', 'PNZ2',
        'PVN2', 'PVN1', 'BV')
afs = rev(afs)

num = setdiff(afs, unique(data_long$source))
if(length(num) > 0){
  num1 = which(afs %in% num)
  afs = afs[-num1]
}

colors = paste(cols1[c(afs, cells)], collapse = '", "')
colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')


data_long$source = factor(data_long$source, levels = afs) 
data_long$target = factor(data_long$target, levels = cells)

nodes <- data.frame(name=c(levels(data_long$source), levels(data_long$target)))
# num = setdiff(levels(data_long$target), unique(data_long$target))
# if(length(num) > 1){
#   num1 = which(nodes$name %in% num)
#   nodes = as.data.frame(nodes[-num1, ])
#   colnames(nodes) = 'name'
# }


# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_long$IDsource = match(data_long$source, nodes$name) - 1 
data_long$IDtarget = match(data_long$target, nodes$name) - 1

#ColourScal = 'd3.scaleOrdinal() .range(["#FFC802", "#CD3278", "#44D1F8", "#4FFFFF", "#0FFF9F", "#3CB371", "#548B54", "#04D002", "#D105D0", "#208FA4", "#203CA4"])'

p = sankeyNetwork(Links = data_long, Nodes = nodes,
                  Source = "IDsource", Target = "IDtarget",
                  Value = "value", NodeID = "name", iterations = 0,
                  sinksRight = F, nodeWidth=40, fontSize=13, nodePadding=20, colourScale=colorJS,
                  height = 2500, width = 900)
p
require(htmlwidgets)
saveWidget(p, file = "allCelltypes_sankeyPlot.html")
