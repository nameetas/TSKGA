af1 = readRDS('visium_2Jun/SNUH_Fig2C_input.rds')
af2 = readRDS('visium_2Jun/SNUHshuffle_Fig2C_input.rds')

af1$dataset = 'orig'
af2$dataset = 'shuffle'

af2$AF = paste0(af2$AF, '_shuffle')

pdf('SNUH_AF_ridgeplot.pdf', width = 8, height = 3)
for(i in 1:length(cells)){
  num = which(af1$celltype == cells[i])
  df1 = af1[num, ]
  num = which(is.na(af1$AF))
  if(length(num) > 0){
    df1 = df1[-num,  ]
  }
  #df1$celltype = 'label'
  
  num = which(af2$celltype == cells[i])
  df2 = af2[num, ]
  num = which(is.na(af2$AF))
  if(length(num) > 0){
    df2 = df2[-num,  ]
  }
  #df2$celltype = 'label'
  
  df3 = rbind.data.frame(df1, df2)
  df3$celltype = ''

  jpeg(paste0(cells[i], '_AFridgeplots.jpeg'), width = 800, height = 300)
  p = ggplot(df3, aes(x = value, y = celltype, fill = AF, color = AF)) + 
    geom_density_ridges(aes(linetype = dataset), alpha = 0) +
    theme_classic() + theme(text = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10, color = 'black')) +
    geom_hline(yintercept = 0, col = 'black', linetype = 'dashed') + 
    geom_vline(xintercept = 2, col = 'black', linetype = 'dashed') +
    theme(legend.position = "none") + xlim(c(-2.5, 5)) +
    theme(plot.title = element_text(size = 80)) + 
    scale_y_discrete(expand = c(0, 0)) + 
    scale_linetype_manual(values = c('orig' = 'solid', 'shuffle' = 'dotted')) +
    scale_color_manual(values = c('BV' = '#FFC802',
                                    'PVN1' = 'violetred3',
                                    'PVN2' = '#941751',
                                    'PNZ1' = '#44D1F8',
                                    'PNZ2' = '#44A1F8',
                                    'PNZ' = '#4FFFFF',
                                    'PAN' = '#0FFF9F',
                                    'PAN1' = 'mediumseagreen',
                                    'PAN2' = 'palegreen4',
                                    'CT' = '#04D002',
                                    'IT' = '#D105D0',
                                    'LE_WM' = '#208FA4',
                                    'LE_GM' = '#203CA4',
                                    'BV_shuffle' = '#FFC802',
                                    'PVN1_shuffle' = 'violetred3',
                                    'PVN2_shuffle' = '#941751',
                                    'PNZ1_shuffle' = '#44D1F8',
                                    'PNZ2_shuffle' = '#44A1F8',
                                    'PNZ_shuffle' = '#4FFFFF',
                                    'PAN_shuffle' = '#0FFF9F',
                                    'PAN1_shuffle' = 'mediumseagreen',
                                    'PAN2_shuffle' = 'palegreen4',
                                    'CT_shuffle' = '#04D002',
                                    'IT_shuffle' = '#D105D0',
                                    'LE_WM_shuffle' = '#208FA4',
                                    'LE_GM_shuffle' = '#203CA4')) + 
    ylab('') + xlab('')
    p = annotate_figure(p, left = text_grob(cells[i], 
                                            color = 'black', rot = 90, vjust = 3))
  print(p)
  
  dev.off()
}
dev.off()
