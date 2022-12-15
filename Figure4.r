##figure4 d bubbleplot showing the representitve markers of identified lienages
if(TRUE) {
library(readxl);library(readr);library(dplyr)
    colors = c('#fcfbfd', '#efedf5', '#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#4a1486')
    colors1 = colorRampPalette(colors)(50)
data= readRDS(file = 'Input_Data/Fig4d_bubblePlotData.rds')
data$features.plot = factor(data$features.plot,
                            levels = rev(c(
'Cav1', 'Pdpn', 'Ager', 'Fabp5', 'Lamp3', 'Sftpc', 'Lyz2', 'Etv5', 'Kng2', 'Lcn2', 'Meg3', 'Ereg', 'Cldn4', 'Krt8', 'Cavin3', 'Krt18', 'Clu', 'Cdkn2a'
                            )))
    fig4d_bubblePlot <-ggplot(data, aes(y = features.plot,x = id)) +        ## global aes
        geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  + 
    scale_fill_gradientn(colors = colors1) +
        labs(x='',y='')+
        scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+             ## to tune the size of circles
        theme(panel.grid.major = element_line(colour = "grey90",size=0.2), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )
print(fig4d_bubblePlot)
}

##figure4 g barplots showing the fraction Krasg12d mutation carrying cells among cell lineages
if(TRUE) {
             data = readr::read_tsv('Input_Data/Fig4g_FractionBarPlotData.txt')
             mycolor = c(AT1 = '#b56cad',
                         AT2 = '#0d7072',
                         'Early tumour/AT2' = '#a59c35',
                         'KAC/KAC-like' = '#de1e26')
             data$celltype = factor(data$celltype,
                                    levels = c('AT1', 'AT2',
                                        'KAC/KAC-like',
                                        'Early tumour/AT2'
                                    ))
        fig4g = ggplot(data = data,aes(y = `celltype`,x=Fraction))+
    geom_bar(width = .3,stat = 'identity',aes(fill = celltype)) +
    xlab('') +
        scale_fill_manual(values = mycolor) +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
        print(fig4g)
}

