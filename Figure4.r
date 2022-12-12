##figure4 a violinplot results generation of monocle2 trajectory inference
if(TRUE) {
    mc2data = read_tsv('')
    fig4a_pseudotime = ggplot(data = mc2data)+
    geom_boxplot(aes(x = `celltype`,y=Pseudotime), 
                         width = .3,
                 alpha = .3,outlier.shape= NA)  +
        geom_box(width = .2) +
    xlab('') +
    theme_classic() + 
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 

    ##differentiation score analysis
    fig4a_cytotrace = ggplot(data = mc2data)+
    geom_boxplot(aes(x = `celltype`,y=cytoTRACE), 
                         width = .3,
                 alpha = .3,outlier.shape= NA)  +
        geom_box(width = .2) +
    xlab('') +
    theme_classic() + 
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 


    ##differentiation score betwen KAC eoe and 7mon 

    fig4a_cytotrace_KAC = ggplot(data = filter(mc2data,celltype == 'KAC'))+
    geom_boxplot(aes(x = `timeopoint`,y=cytoTRACE), 
                         width = .3,
                 alpha = .3,outlier.shape= NA)  +
        geom_box(width = .2) +
    xlab('') +
    theme_classic() + 
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
}

##figure4 d bubbleplot showing the representitve markers of identified lienages
if(TRUE) {
library(readxl);library(readr);library(dplyr)
setwd('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/figures/figures_data/epithelial/EpithelialUmap/Epithelial_normalCells/alveolar_airway_andAllfiguresworkdir')
data_top50 = read_tsv('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/figures/figures_data/epithelial/EpithelialUmap/Epithelial_normalCells/alveolar_airway_andAllfiguresworkdir/Fig5-mousedata-bubbleplot-markers-cells-2022-03-01.tsv') %>%
    group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC)
data_sel = read_excel('Fig5-mousedata-bubbleplot-markers-sel.xlsx', sheet = 1);
data_sel = filter(data_sel, !(gene %in% c('Eno1','Cdkn2a','Gnk2','Krt18','Krt7','Cldn18',
                                  'Mboat1','Pdpn','Cdk1','Cenpf','Napsa')) )
rbind(data_top50, data_sel) -> mouseMarkers
write_tsv(mouseMarkers, paste0('TableS9-mouseMarkers-',Sys.Date(),'.tsv'))

obj = sct.data1
gene0 = read_excel('Fig5-mousedata-bubbleplot-markers-sel.xlsx', sheet = 1);
gene0 = filter(gene0, !(gene %in% c('Eno1','Cdkn2a','Gnk2','Krt18','Krt7','Cldn18',
                                  'Mboat1','Pdpn','Cdk1','Cenpf','Napsa')) )
p1 = c('AT1',
       'AT2',
       'Alveolar progenitor',
       'Basal',
       'Ciliated',
       'Club and Secretory',
       'Neuroendocrine',
       'Proliferating',
       'Tuft',
       'Tumor cells')
gene0$cluster = factor(gene0$cluster,
                       levels = p1,)
gene0 = arrange(gene0, cluster,p_val_adj)
gene1 = gene0$gene
obj = ScaleData(obj, features = unique(gene1))
setdiff(gene1, rownames(obj))
p<-DotPlot(obj, features = unique(gene1))
data = p$data

data$id = factor(data$id,
                 levels = c(
                   p1
                 ))
data = na.omit(data)
data$features.plot = factor(as.character(data$features.plot),
                            levels = unique(gene1))

colors = c('#efedf5','#bcbddc','#756bb1')
colors = c('#fcfbfd', '#efedf5', '#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#4a1486')
colors1 = colorRampPalette(colors)(50)
barp <-ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
  #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
  ## scale_fill_viridis(option = 'plasma')+
  ## scale_fill_brewer(direction = -1) +
  ## scale_fill_gradient(low = "grey60", high = "#e48759") +
  scale_fill_gradientn(colors = colors1) +
  labs(x='',y='')+
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+             ## to tune the size of circles
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )
ggsave(paste0('Fig5-mouse-umap-bubblePlot',Sys.Date(),'.pdf'),barp,height=4,width=10,useDingbats = F)
}

##figure4 e barplots showing the fraction changes among the timepoint and treatment
if(TRUE) {
    fractionDat = read_tsv('')
    fig4e = ggplot(data = fractionDat)+
    geom_bar(aes(x = `Sample`,y=Fraction), 
                         width = .3,
                 alpha = .3,outlier.shape= NA)  +
    xlab('') +
    theme_classic() + facet_wrap(~celltype) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
}

##figure4 g barplots showing the fraction Krasg12d mutation carrying cells among cell lineages
if(TRUE) {
     ##bottome part of the panel showing the fraction of G12D mutaiton cell fraction
        data = read_excel('SuppTable.xlsx',sheet = 11)
        fig4g = ggplot(data = data)+
    geom_bar(aes(x = `Group`,y=Fraction), 
                         width = .3)
        geom_label(`number of cells with Kras G12D mutation`) + 
    xlab('') +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
}





