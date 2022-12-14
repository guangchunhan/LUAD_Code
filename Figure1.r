##figure 1
##figure 1a was a plot for illustratoin purpose, it's generated with biorender, ms ppt and illustrator

##figure 1b bubble plot
if(TRUE) { ##code for figure 1b bubble plot
    library(ggplot2)
    library('viridis')
    ## Fig1b-markerGene-expression.rds: bubble plot data object generated from seurat object, seurat object too big to be uploade.
    data  = readRDS('Input_data/Fig1b-markerGene-expression-v1.rds')
    data = filter(data, !(id %in% c('proliferating_alveolar','proliferating_basal','proliferating_club')))
    p1 = c('Malignant','AT1',  'AT2',
           'AIC','SDP','Tuft','Ciliated','NE', 'Ionocyte','Basal','Club')
    data$id = factor(data$id,
                     levels = c(
                         p1
                     ))
    data$features.plot = factor(as.character(data$features.plot),
                                levels = c(
                                    'AGER','CAV1','PDPN','CLIC5', ##AT1, CLIC5 in lung cancer altlas paper
                                    'SFTPB','SFTPC','KRT7','NKX2-1','HBEGF', ## AIC highly expressed
                                    'SFTPD','ETV5','WIF1','HHIP', ##AT2 cells 
                                    'ASCL2','MGST2','PTGS1','GNAT3',##tuft
                                    'CAPS','PIFO','FOXJ1','CCDC78','TP73', ##ciliated
                                    'CALCA','ASCL1','CHGA','ENO2','SYP', ##neuroendocrine
                                    'ASCL3','FOXI1','CEL','ATP6V1B1',##ionocyte
                                    'KRT17','TP63','KRT5','KRT15', ##basal cells
                                    'SCGB3A2','SCGB1A1','SCGB3A1','WFDC2','CYP2F1' ##Club and secretory
                                ))

    colors = c('#efedf5','#bcbddc','#756bb1')
    colors = c('#fcfbfd', '#efedf5', '#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#4a1486')
    colors1 = colorRampPalette(colors)(50)
    Fig1b_bubblePlot <-ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
        geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +  
        scale_fill_gradientn(colors = colors1) +
        labs(x='',y='')+
        scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+             
        theme(panel.grid.major = element_line(colour = "grey90",size=0.2),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )
    print(Fig1b_bubblePlot)
    ## ggsave(paste0('Fig1B-bubblePlot.pdf'),Fig1b_bubblePlot,height=4,width=10,useDingbats = F)

    ##figure 1b cell type hierarachical tree
    ##Fig1b-hc.rds: PCA input data for hierachical clustering
    mat.pca = readRDS('Input_data/Fig1b-hc-v1.rds')
    library(dendextend)
    library(circlize)
    hcAirway = hclust(dist(t(as.matrix((mat.pca[,c("Basal","Ciliated","Club","Ionocyte",
                                                   "NE","Tuft")])))),method = 'ward.D')
    hcAirwayd = as.dendrogram(hcAirway)
    hcAirwayd %>% set("leaves_pch", 19)  -> hcAirwayd
    pdf(paste0('Fig1bAirwayCells-dendplot-',Sys.Date(),'.pdf'),useDingbats = F)
    par(mar = c(20,2,1,1))
    plot(hcAirwayd,type = 'triangle')
    colored_bars(colors = rep('black',6),dend = hcAirwayd, sort_by_labels_order = F)
    dev.off()


    hcAlveoli = hclust(dist(t(as.matrix((mat.pca[,c("AT1","AT2","AIC","SDP")])))),method = 'ward.D')
    hcAlveolid = as.dendrogram(hcAlveoli)
    hcAlveolid %>% set("leaves_pch", 19)  -> hcAlveolid
    pdf(paste0('Fig1bAlveoliCells-dendplot-',Sys.Date(),'.pdf'),useDingbats = F)
    par(mar = c(20,2,1,1))
    plot(hcAlveolid,type = 'triangle')
    colored_bars(colors = rep('black',4),dend = hcAlveolid, sort_by_labels_order = F)
    dev.off()
}

##PCA visualizationin fig1e, UMAPs in figure 1c ,1d, 1f,1g,1k umap visualization plots were generated from seurat object and imported into JMP Pro v15.

##figure 1f,g: violinplots
if(TRUE) {
    ##fig1F: violinplots of cnv score and smoking status
    vlnData = readRDS('fig1f_g_violnPlot-v1.rds') ## data with cnvscore/cytotrace score for malignant cells 
    vlnData$smoking1 = factor(vlnData$smoking1, levels = c("smoker","never smoker"))
    smokingstatus=c('#509109','#ECE807'); names(smokingstatus)=c("never smoker","smoker");
    Fig1f_CNV_VlnBySmoke = ggplot(data = vlnData,aes(x=`smoking1`,y=`aviv_cnvscore_hg38`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  smoking1)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('inferCNVscore') +
        theme_classic() +
        scale_fill_manual(values = smokingstatus) + 
        theme(axis.text.x = element_text(angle = 90)) + theme_classic()
    print(Fig1f_CNV_VlnBySmoke)

    ##show cnv score violin plot by mutation status three categories
    mutationcolor1= c('KRAS' = '#ae3135','EGFR' = '#3d74b6','EGFR_KRAS_WT' = '#b4b4b5') 
    Fig1f_CNV_VlnByMutation = ggplot(data = vlnData,aes(x=`mutation1`,y=`aviv_cnvscore_hg38`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  mutation1)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('inferCNVscore') +
        theme_classic() +
        scale_fill_manual(values = mutationColor1) + 
        theme(axis.text.x = element_text(angle = 90)) + theme_classic()

    ##Fig1g: cytoTRACE score vlnplot by mutation group
    vlnData1 = vlnData;vlnData1$mutation1 = factor(vlnData1$mutation1,
                                                   levels = c('EGFR_KRAS_WT',
                                                              'EGFR','KRAS'))
    Fig1gCyto_VlnByMutation = ggplot(data = vlnData1,aes(x=`mutation1`,y=`CytoTRACEScore`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  mutation1)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('cytoTRACEscore') +
        theme_classic() +
        scale_fill_manual(values = mutationColor1) + 
        theme(axis.text.x = element_text(angle = 90)) + theme_classic()

    library(ggridges)
    ##Fig1G Ridge plot of differentiation score( cytotrace score) by patient 
    sampleColor = c('#737fb0','#c16360','grey93', '#a38861',
                    'grey78','#3d74b6', '#9ebcde', '#bfd8f0',
                    '#c69a58', '#e79941','grey63',
                    'grey48', '#ac3135','#c4cee7')
    names(sampleColor) =  paste0('P',c(1:10,11,13:15))
    vlnData$patient = factor(vlnData$patient, levels = rev(c('P10','P14','P2','P8',
                                                             'P1','P7','P15','P6',
                                                             'P4','P9',
                                                             'P13','P3','P11','P5')))
    Fig1g_RidgePlot = ggplot(data = vlnData,aes(x = CytoTRACEScore, y = patient)) +
        ggridges::geom_density_ridges(aes(fill = patient),alpha = .8) +
        scale_fill_manual(values = sampleColor) + 
        scale_x_reverse() +
        labs(x = '',
             y = '') +
        ggridges::theme_ridges(grid = T) 
}

##figure 1h,i,j: metaprogram analysis
if(TRUE) {
    ##Fig1h: unsupservised way to learn about metaprogram score using the heatmap
    metaprogramData=readRDS('Input_Data/Fig1h_heatmapData.rds')
    MPmat = metaprogramData[[1]]; ##contains the metaprogram scores derived from cytotrace analysi
    sDat = metaprogramData[[2]] ##contains the sample meta information
    sweep(MPmat, 1, apply(MPmat, 1, mean)) -> mat1
    mat1 = sweep(mat1, 1, FUN = '/',STATS = apply(MPmat, 1 ,sd))
    mat1[ which(mat1 > 2,2) ] = 2
    mat1[ which(mat1 < -2,2) ] = -2
    mutationColor= c('KRAS' = '#ae3135','EGFR' = '#3d74b6','MET' = '#78c679','Other' = '#b4b4b5') ##updated MET color as suggested by linghua
    smokingstatus=c('#509109','#ECE807')
    names(smokingstatus)=c("never smoker","smoker");
    sampleColor = c('#737fb0','#c16360','grey93', '#a38861',
                    'grey78','#3d74b6', '#9ebcde', '#bfd8f0',
                    '#c69a58', '#e79941','grey63',
                    'grey48', '#ac3135','#c4cee7')
    names(sampleColor) =  paste0('P',c(1:10,11,13:15))
    aCol = list(mutation = mutationColor,
                smoking1=smokingstatus,
                patient = sampleColor,aviv_cnvscore_hg38 = c('#243160','#efdf6a'),
                CytoTRACEScore = colorRampPalette(c('#642b6c','#6fb657','#ba3330'))(50))
    callback = function(hc, mat){
        sv = svd(t(mat))$v[,1]
        dend = reorder(as.dendrogram(hc), wts = sv)
        as.hclust(dend)
    }
    library(pheatmap)
    pheatmap(mat1,
             scale = 'none',
             col = colorRampPalette((c('#4b509c','white','#d24945')))(100),
             annotation_col=sDat,
             annotation_colors = aCol,
             clustering_method = 'ward.D',
             cluster_cols = T,
             cluster_rows = T,
             cutree_cols = 5,
             show_colnames = F,height = 12,width = 18,
             border_color = NA,filename = paste0('Fig1h-metaprogram-heatmap.png')) 

    ##figure 1i and 1j
    fig1i_j_data = readRDS('Input_Data/Fig1i_j_MPVlnData.rds') ##long format data of metaprogram scores for malignant cells
    library(ggplot2)
    Fig1i_j_VlnPlots = ggplot(data = fig1i_j_data,aes(x=`k5_rename`,y=`MPscore`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  k5_rename)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('score') +
        theme_classic() +
        ##scale_fill_manual(values = colors) + 
        facet_wrap(~MPs,scales = 'free') + 
        theme(axis.text.x = element_text(angle = 90))
    print(Fig1i_j_VlnPlots)
    ggsave('Figure1i_j-selectedMPs-VlnPlot.pdf',plots,width = 7,
           height = 12)
}

##figure 1k: visualization results of malignant cell analysis of P14 incuding umap and mc2 dimension reduction was loaded into JMP Pro v15 and visualized accordingly.




