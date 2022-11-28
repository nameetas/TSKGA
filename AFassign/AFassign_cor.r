cpmCalc = function(counts){
  i = 1
  while (i <= ncol(counts)) {
    j = i + 10000
    if(j > ncol(counts)){
      j = ncol(counts)
    }
    counts1 = as.matrix(counts[, i:j])
    mat = apply(counts1, 2, function(x)(x/sum(x))*1000000)
    mat = log2(mat + 1)
    if(i == 1){
      cpm = mat
    } else {
      cpm = cbind.data.frame(cpm, mat)
    }
    i = j + 1
  }
  return(cpm)
}

library(matrixStats)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)

AFassign = function(obj, method, thres = 0.15){
  #obj: Seurat Object
  #method: FFPE/Frozen
  #thres: 0.15
  if(method == 'FFPE'){
    mat = readRDS('AFavg_ffpe.rds')
  } else if(method == 'Frozen'){
    mat = readRDS('AFavg_frozen.rds')
  } else {
    print('Incorrect input for method')
    return(0)
  }
  counts = obj@assays$Spatial@counts
  cpm = cpmCalc(counts)
  print('CPM calculation done')

  genes = intersect(rownames(mat), rownames(cpm))
  cpm1 = cpm[genes, ]
  mat = mat[genes, ]

  identical(rownames(cpm1), rownames(mat))
  mat = (mat -rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]
  M = cor(cpm1, mat)

  meta = obj@meta.data
  identical(rownames(meta), rownames(M))

  info = meta[,1:2]
  info[, 1] = NA
  info[, 2] = NA
  colnames(info) = c('AF', 'maxAF_cor')

  for(i in 1:nrow(M)){
    if(length(which(is.na(M[i, ]))) > 7){
      next
    }
    num = which(M[i,] == max(M[i,], na.rm = T))
    info[i, 'maxAF_cor'] = M[i, num]
    info[i, 'AF'] = colnames(M)[num]
  }
  num = which(info$maxAF_cor < thres)
  info[num, 'AF'] = NA
  write.csv(info, file = 'AFassign.csv')
}
