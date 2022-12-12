##figure 1

##R version:4.1.0

##figure 1b bubble plot
if(TRUE) { ##code for figure 1b bubble plot
    library(Seurat)
    library(readxl);library(readr);library(dplyr)
    obj = readRDS('allCellSeuratObj-.rds')
    library(ggplot2)
    library('viridis')
    DefaultAssay(obj)='RNA'
    obj = NormalizeData(obj, normalization.method = 'LogNormalize', scale.factor = 10000)
    ##genes below are selected from the marker gene lists
    gene = c(
        'AGER','CAV1','PDPN','CLIC5', 
        'SFTPB','SFTPC','KRT7','NKX2-1','HBEGF', 
        'SFTPC','SFTPD','ETV5','WIF1','HHIP', ##AT2 cells, ETV5 is used in nature paper
        'KRT17','TP63','KRT5','KRT15', ##basal cells
        'CAPS','PIFO','FOXJ1','CCDC78','TP73', ##ciliated
        'SCGB3A2','SCGB1A1','SCGB3A1','WFDC2','CYP2F1', ##Club and secretory 
        'CALCA','ASCL1','CHGA','ENO2','SYP', ##neuroendocrine
        'ASCL2','MGST2','PTGS1','GNAT3',##tuft
        'ASCL3','FOXI1','CEL','ATP6V1B1'##ionocyte
    )
    gene1 = gene
    obj = ScaleData(obj, features = unique(gene1))
    setdiff(gene1, rownames(obj))
    p<-DotPlot(obj, features = unique(gene1))
    data = p$data
    p1 = c('alveoli_AT1',  'alveoli_AT2',
           'alveoli_AVP','alveoli_Bronchio',
           'airway_Tufft','airway_Ciliated','airway_Neuroendocrine', 'airway_Ionocyte','airway_Basal',
           'airway_ClubandSecretory','malignant',
           'proliferating_alveolar','proliferating_basal',
           'proliferating_club')
    data$id = factor(data$id,
                     levels = c(
                         p1
                     ))
    data = na.omit(data)
    data$features.plot = factor(as.character(data$features.plot),
                                levels = c(
                                    'AGER','CAV1','PDPN','CLIC5', ##AT1 cells, hopx is not very specific, use CLIC5 which also be used in nature lung cancer altlas paper
                                    'SFTPB','SFTPC','KRT7','NKX2-1','HBEGF', ##Alveolar progenitor (combine of AT cells and NKX2-1, HB-EGF etc)
                                    'SFTPD','ETV5','WIF1','HHIP', ##AT2 cells, ETV5 is used in nature paper
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
    ggsave(paste0('Fig1B-bubblePlot.pdf'),barp,height=4,width=10,useDingbats = F)


    ##figure 1b cell type hierarachical tree clustering
    mat.pseudobulk = readRDS('NormalEpithelial-bulkHarmony-P8P9-correct-2022-01-03.rds')
    library(dendextend)
    library(circlize)
    hcAirway = hclust(dist(t(as.matrix((mat.pseudobulk[,c("airway_Basal","airway_Ciliated","airway_ClubandSecretory","airway_Ionocyte",
                                                          "airway_Neuroendocrine","airway_Tufft")])))),method = 'ward.D')
    hcAirwayd = as.dendrogram(hcAirway)
    hcAirwayd %>% set("leaves_pch", 19)  -> hcAirwayd
    pdf(paste0('NormalEpithelialAirwayCells-dendplot-',Sys.Date(),'.pdf'),useDingbats = F)
    par(mar = c(20,2,1,1))
    plot(hcAirwayd,type = 'triangle')
    colored_bars(colors = rep('black',6),dend = hcAirwayd, sort_by_labels_order = F)
    dev.off()


    hcAlveoli = hclust(dist(t(as.matrix((mat.pseudobulk[,c("alveoli_AT1","alveoli_AT2","alveoli_AVP","alveoli_Bronchio")])))),method = 'ward.D')
    hcAlveolid = as.dendrogram(hcAlveoli)
    hcAlveolid %>% set("leaves_pch", 19)  -> hcAlveolid
    pdf(paste0('NormalEpithelialAlveoliCells-dendplot-',Sys.Date(),'.pdf'),useDingbats = F)
    par(mar = c(20,2,1,1))
    plot(hcAlveolid,type = 'triangle')
    colored_bars(colors = rep('black',4),dend = hcAlveolid, sort_by_labels_order = F)
    dev.off()
}

##figure 1c and d umap visualization plots
if(TRUE) {
    library(Seurat)
    library(readxl);library(readr);library(dplyr)
    library(ggplot2)
    library('viridis')
    ##read in the umap dimention reduction data
    dimUmap = read_tsv('EpiPosumap_malignantCells_n17064recent_withP8P9hg38_cleanded_no_batchCorrection2022-03-07.tsv')

    umap$patient = gsub('_.*','',gsub('T_.*','',gsub('-.*$','',umap$orig.ident)))
    df = data.frame(patient = paste0('P',c(1:10,11,13:15)),
                    mutation = c('EGFR','KRAS','Other','MET','Other','EGFR','EGFR','EGFR','MET','KRAS',
                                 'Other','Other','KRAS','EGFR'))
    umap = left_join(umap,df, by = 'patient') ##smoking status
    umap$patient = factor(umap$patient, levels =paste0('P',c(1:10,11,13:15) ))
    mutationColor= c('KRAS' = '#ae3135','EGFR' = '#3d74b6','MET' = '#78c679','Other' = '#b4b4b5') ##updated MET color as suggested by linghua
    g <- ggplot(umap,
                aes(x = UMAP_1, y=UMAP_2, colour = mutation)) + scale_color_manual(values = mutationColor) +
        geom_point(show.legend = T, size = .3) + theme_classic()  + guides(colour = guide_legend(override.aes = list(size=10)))

    ## blank_theme
    ggsave(paste0('Fig1d-bymutationCat-',Sys.Date(),'.png'),
           g,  height = 4,width = 5,units = 'in',
           dpi = 300)
    ##by clusterID
    umap$`snn_default_param_res.0.4` = paste0('c',umap$`snn_default_param_res.0.4`)
    g <- ggplot(umap,
                aes(x = UMAP_1, y=UMAP_2, colour = `snn_default_param_res.0.4`))+
        ## scale_color_manual(values = smokingstatus) +
        geom_point(show.legend = T, size = .3) + theme_classic()  + guides(colour = guide_legend(override.aes = list(size=10)))

    ## blank_theme
    ggsave(paste0('Fig1d-byClusterID-',Sys.Date(),'.pdf'),
           g)
}

##figure 1e
if(TRUE) {
    ## scatterplot for malignant cells 3d of PCAs
    library(readr);library(readxl)
    maligScatterData = read_excel('malignantCells_scatterplot_3d_P8P9Hg38Ref.xlsx')
    ##load updated nc17064 pca
    obj = readRDS('EpiPos_malignantCells_n17064recent_cleanded_withP8P9hg38_2022-03-07.rds')
    pca = Embeddings(obj,reduction = 'pca'); pca = pca[,1:3]
    cellids = intersect(maligScatterData$ID,rownames(pca)); pos = na.omit(match(maligScatterData$ID,cellids))
    maligScatterData = maligScatterData[pos,]
    match(maligScatterData$ID, rownames(pca)) -> pos1
    maligScatterData$PC_1 = pca[pos1, 'PC_1'];maligScatterData$PC_2 = pca[pos1, 'PC_2'];maligScatterData$PC_3 = pca[pos1, 'PC_3']
    maligScatterData$`mutation category`[maligScatterData$`mutation category` == 'EML4-ALK'] = 'Other' 
    maligScatterData$`mutation category`[maligScatterData$`patient ID` %in% c("P4","P9")] = 'MET'
    write_tsv(maligScatterData, paste0('malignantCells_scatterplot_3d_P8P9Hg38Ref_n17064recent-',Sys.Date(),'.tsv')) ## use jmp to generate pca plot in figure 2
    library(scatterplot3d)
    set.seed(417)
    library(plotly)
    mycol=c('#BE1E2D','#1B75BC','#B4B4B4','#AA875B'); names(mycol) = c('KRAS','EGFR','Other','MET')
    fig <- plot_ly(maligScatterData, x = ~PC_3, y = ~PC_2, z = ~PC_1, color = ~`mutation category`,
                   ##fig <- plot_ly(maligScatterData, x = ~PC_1, y = ~PC_3, color = ~`mutation category`,
                   colors =mycol)
    fig <- fig %>% add_markers(size = 3) 
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                       yaxis = list(title = 'PC2'),zaxis = list(title = 'PC3')
                                       ))

    fig
}

##figure 1f,g: harmony corrected version of the malignant cell UMAP
if(TRUE) {
    ##the umap dimension reduction results were stored in a data frame:
    umap = read_tsv(paste0('Epithelial_malignant-n17064recent_-Harmony-UMAP_0.1_2022-03-07.tsv'))
    umap$patient = gsub('_.*','',gsub('T_.*','',gsub('-.*$','',umap$orig.ident)))
    df = data.frame(patient = paste0('P',c(1:10,11,13:15)),
                    mutation = c('EGFR','KRAS','Other','MET','Other','EGFR','EGFR','EGFR','MET','KRAS',
                                 'Other','Other','KRAS','EGFR'))
    umap = left_join(umap,df, by = 'patient')
    sampleColor = c('#737fb0','#c16360','grey93', '#99d8c9', ##udpate MET color
                    'grey78','#3d74b6', '#9ebcde', '#bfd8f0',
                    '#addd8e', '#e79941','grey63', ##udpate MET color
                    'grey48', '#ac3135','#c4cee7')
    names(sampleColor) =  paste0('P',c(1:10,11,13:15))
    umap$UMAP_2 = 0-umap$UMAP_2
    g <- ggplot(umap,
                aes(x = UMAP_1, y=UMAP_2, colour = patient)) +
        geom_point(show.legend = T, size = .3) +
        scale_color_manual(values = sampleColor) + theme_classic()
    ## blank_theme
    ggsave(paste0('Fig1f-left-Harmony-UMAP-bymutation-',Sys.Date(),'.png'),
           g,  height = 9,width = 10,units = 'in',
           dpi = 300)

    ##cnv score and smoking status
    cnvdat = read_tsv(paste0('Epithelial_malignant-Harmony-UMAP_0.1_2022-01-03_cnvscore.tsv')) %>% select(-CytoTRACEScore) ##
    cnvdat$patient = gsub('_.*','',gsub('T_.*','',gsub('-.*$','',cnvdat$orig.ident))); cnvdat = select(cnvdat, ID, aviv_cnvscore_hg38)
    umap = left_join(umap,cnvdat, by = c('ID' = 'ID') ) 
    dat = read_excel('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/figures/figures_data/clinicalvariableplot/sample_clinical.xlsx')
    smokepackyr = read_excel('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/figures/figures_data/clinicalvariableplot/P1 P16 smoking details_Updated_12_22_2021.xlsx') %>% select(-Gender)
    dat = left_join(dat, smokepackyr,
                    by = c('Patient #' = 'scRNAseq P#'))
    umap = left_join(umap,dat, by = c('patient' = 'Patient #') ) ##dat is smoking

    plots = ggplot(data = umap,aes(x=`smoking1`,y=`aviv_cnvscore_hg38`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  smoking1)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('inferCNVscore') +
        theme_classic() +
        scale_fill_manual(values = smokingstatus) + 
        theme(axis.text.x = element_text(angle = 90)) + theme_classic()
    ggsave('figure2-malignantCells-n17064recent-cnvscore-bySmokingstatus.pdf',plots,width = 4,
           height = 4)

    ##show cnv score violin plot by mutation status
    mutationColor= c('KRAS' = '#ae3135','EGFR' = '#3d74b6','MET' = '#78c679','Other' = '#b4b4b5') ##updated MET color as suggested by linghua
    plots = ggplot(data = umap,aes(x=`mutation`,y=`aviv_cnvscore_hg38`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  mutation)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('inferCNVscore') +
        theme_classic() +
        scale_fill_manual(values = mutationColor) + 
        theme(axis.text.x = element_text(angle = 90)) + theme_classic()
    ggsave('figure2-malignantCells-n17064recent-cnvscore-byMutationStatus.pdf',plots,width = 4,
           height = 4)

    ##show cnv score violin plot by mutation status three categories
    umap1 = umap;umap1$mutation1 = umap$mutation; umap1$mutation1[umap$mutation %in% c('MET','Other')] = 'EGFR_KRAS_WT'
    mutationColor1= c('KRAS' = '#ae3135','EGFR' = '#3d74b6','EGFR_KRAS_W' = '#b4b4b5') ##updated MET color as suggested by linghua
    plots = ggplot(data = umap1,aes(x=`mutation1`,y=`aviv_cnvscore_hg38`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  mutation1)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('inferCNVscore') +
        theme_classic() +
        scale_fill_manual(values = mutationColor1) + 
        theme(axis.text.x = element_text(angle = 90)) + theme_classic()
    ggsave('figure2-malignantCells-n17064recent-cnvscore-byMutationStatus-three-categories.pdf',plots,width = 4,
           height = 4)

    ##cytoTRACE score
    ##joint cytotrace score
    cytotrace = read_tsv('cytotrace-EpiPosTumorSamplesCells-P8P9hg38Ref-2022-01-06.tsv')
    cytotrace$ID = gsub('^_','',gsub('-1$','',cytotrace$ID))
    umap = left_join(umap,cytotrace, by = c('ID' = 'ID'))
    ##write_tsv(umap[,-c(78)], paste0('figure2-malignantCells-17216-allinone-table-',Sys.Date(),'.tsv'))
    write_tsv(umap[,-c(78)], paste0('figure2-malignantCells-n17064recent-allinone-table-',Sys.Date(),'.tsv'))
    ##show cytotrace score by mutation group
    plots = ggplot(data = umap,aes(x=`mutation`,y=`CytoTRACEScore`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  mutation)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('cytoTRACEscore') +
        theme_classic() +
        ## scale_fill_manual(values = smokingstatus) + 
        theme(axis.text.x = element_text(angle = 90)) + theme_classic()
    ggsave('figure2-malignantCells-n17064recent-cytoTRACescore-byMutationstatus.pdf',plots,width = 4,
           height = 4)

    ##joint cytotrace score with three categories
    cytotrace = read_tsv('cytotrace-EpiPosTumorSamplesCells-P8P9hg38Ref-2022-01-06.tsv')
    cytotrace$ID = gsub('^_','',gsub('-1$','',cytotrace$ID))
    umap1 = left_join(umap1,cytotrace, by = c('ID' = 'ID'))
    ##write_tsv(umap[,-c(78)], paste0('figure2-malignantCells-17216-allinone-table-',Sys.Date(),'.tsv'))
    ##write_tsv(umap[,-c(78)], paste0('figure2-malignantCells-n17064recent-allinone-table-',Sys.Date(),'.tsv'))
    ##show cytotrace score by mutation group
    plots = ggplot(data = umap1,aes(x=`mutation1`,y=`CytoTRACEScore`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  mutation1)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('cytoTRACEscore') +
        theme_classic() +
        scale_fill_manual(values = mutationColor1) + 
        theme(axis.text.x = element_text(angle = 90)) + theme_classic()
    ggsave('figure2-malignantCells-n17064recent-cytoTRACescore-byMutationstatus-three-category.pdf',plots,width = 4,
           height = 4)
    library(ggridges)
    ##show cytotrace score by smoking group
    plots = ggplot(data = umap,aes(x=`smoking1`,y=`CytoTRACEScore`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  smoking1)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('cytoTRACEscore') +
        theme_classic() +
        scale_fill_manual(values = smokingstatus) + 
        theme(axis.text.x = element_text(angle = 90)) + theme_classic()
    ggsave('figure2-malignantCells-n17064recent-cytoTRACescore-bySmokingstatus.pdf',plots,width = 4,
           height = 4)

    ##fig1g, ridge plots by patient
    ##cytotrace ridge plot
    gp1.fi = ggplot(data = umap,aes(x = CytoTRACEScore, y = mutation)) +
        geom_density_ridges(aes(fill = mutation),alpha = .8) +
        ##scale_fill_manual(values = c('#abdda4',rep('#f46d43',14))) + 
        labs(x = '',
             y = '') +
        theme_ridges(grid = T) 


    ggsave(paste0('Figure2-n17064recent-cnv-cytotrace-by-mutationstatus-ridgePlot-',Sys.Date(),'.pdf'),gp1.fi,
           useDingbats = F, width = 8,height = 10)

    ##rigde plot by patient 
    sampleColor = c('#737fb0','#c16360','grey93', '#a38861',
                    'grey78','#3d74b6', '#9ebcde', '#bfd8f0',
                    '#c69a58', '#e79941','grey63',
                    'grey48', '#ac3135','#c4cee7')
    names(sampleColor) =  paste0('P',c(1:10,11,13:15))
    umap$patient = factor(umap$patient, levels = rev(c("P2","P10","P14",
                                                       "P1","P6","P7","P8","P15",
                                                       "P4","P9","P3",'P5','P11',"P13")))
    gp1.fi = ggplot(data = umap,aes(x = CytoTRACEScore, y = patient)) +
        geom_density_ridges(aes(fill = patient),alpha = .8) +
        scale_fill_manual(values = sampleColor) + 
        labs(x = '',
             y = '') +
        theme_ridges(grid = T) 


    ggsave(paste0('Figure2-n17064recent-cnv-cytotrace-by-mutationstatus-ridgePlot-eachpatient-',Sys.Date(),'.pdf'),gp1.fi,
           useDingbats = F, width = 4,height = 8)

    ##order by level of cytotrace
    group_by(umap, patient,mutation) %>% summarise(meanCytotrace = mean(CytoTRACEScore)) ->patientLvlCytoTRACE
    patientLvlCytoTRACE$mutation = factor(patientLvlCytoTRACE$mutation, levels = c('KRAS','EGFR','MET','Other'))
    umap$patient = factor(umap$patient, levels = rev(arrange(patientLvlCytoTRACE, mutation,desc(meanCytotrace)) %>% .$patient))
    gp1.fi = ggplot(data = umap,aes(x = CytoTRACEScore, y = patient)) +
        geom_density_ridges(aes(fill = patient),alpha = .8) +
        scale_fill_manual(values = sampleColor) + 
        labs(x = '',
             y = '') +
        theme_ridges(grid = T) 


    ggsave(paste0('Figure2-n17064recent-cnv-cytotrace-by-mutationstatus-ridgePlot-eachpatient-order-bylvl-',Sys.Date(),'.pdf'),gp1.fi,
           useDingbats = F, width = 4,height = 8)

    ##checked on 06/30/2022
    umap1$patient = factor(umap1$patient, levels = rev(c("P2","P10","P14",
                                                         "P1","P6","P7","P8","P15",
                                                         "P4","P9","P3",'P5','P11',"P13")))
    gp1.fi = ggplot(data = umap1,aes(x = CytoTRACEScore, y = patient)) +
        ggridges::geom_density_ridges(aes(fill = patient),alpha = .8) +
        scale_fill_manual(values = sampleColor) + 
        labs(x = '',
             y = '') +
        ggridges::theme_ridges(grid = T) 


    ggsave(paste0('Figure2-n17064recent-cnv-cytotrace-by-mutationstatus-ridgePlot-eachpatient-order-bylvl-',Sys.Date(),'.pdf'),gp1.fi,
           useDingbats = F, width = 4,height = 8)

    ##boxplot of cytotrace between EGFR and KRAS, others
    ##cytotrace ridge plot
    gp1.fi = ggplot(data = umap,aes(y = CytoTRACEScore, x = mutation)) +
        geom_violin(aes(fill = mutation),alpha = .8) +
        geom_boxplot(width = .2, outlier.shape = NA) +
        ##scale_fill_manual(values = c('#abdda4',rep('#f46d43',14))) + 
        labs(x = '',
             y = '') +
        theme_classic() 


    ggsave(paste0('Figure2-n17064recent-cnv-cytotrace-by-mutationstatus-violin-',Sys.Date(),'.pdf'),gp1.fi,
           useDingbats = F, width = 6,height =3)
}

##figure 1h,i,j: metaprogram analysis
if(TRUE) {
##Fig2x: unsupservised way to learn about metaprogram score
    umap = read_tsv('Figure2-n17064recent-cnv-cytotrace-vs-smoking-mutationstatus-2022-03-14.tsv') ##contains Meta program scores
MPmat = select(umap, one_of(paste0('MP_',c(1:21,30,31)))) %>% as.matrix(); rownames(MPmat) = umap$ID
sDat = select(umap, patient,smoking1, mutation,aviv_cnvscore_hg38,CytoTRACEScore)  %>% as.data.frame(); rownames(sDat) = umap$ID
MPmat = t(MPmat)
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
         ##         color = colorramppalette(c("navy", "grey60", "yellow"))(50),
         ##col = colorRampPalette(rev(c('red','orange','yellow','black')))(100),
         col = colorRampPalette((c('#4b509c','white','#d24945')))(100),
         annotation_col=sDat,
        ## clustering_callback = callback,
         annotation_colors = aCol,
         clustering_method = 'ward.D',
         cluster_cols = T,
         ## gaps_col = gapscol,
         cluster_rows = T,
         cutree_cols = 5,
         show_colnames = F,height = 12,width = 18,
         border_color = NA,filename = paste0('Fig1h-metaprogram-heatmap-',Sys.Date(),'.pdf') ) 

    ##figure 1i and 1j
    ##fig2: after separate into 5 groups, compare MP3, 21 and 30
vlndat = alluvidat0; vlndat$k5_rename = rep(NA,nrow(vlndat))
vlndat$k5_rename[vlndat$k5 == 5] = 'c2';
vlndat$k5_rename[vlndat$k5 == 2] = 'c4';
vlndat$k5_rename[vlndat$k5 == 3] = 'c1';
vlndat$k5_rename[vlndat$k5 == 1] = 'c3';
vlndat$k5_rename[vlndat$k5 == 4] = 'c5';
vlndat1 = cbind(vlndat[,c('aviv_cnvscore_hg38', 'CytoTRACEScore','k5_rename'),drop = F], t(MPmat)[rownames(vlndat),]);
vlndat1.tibble = data.frame(ID = rownames(vlndat1),vlndat1)
##write_tsv(vlndat1.tibble, paste0('Fig2-malignantCells-nc17261-MPscores-heamapgrps-',Sys.Date(),'.tsv'))
write_tsv(vlndat1.tibble, paste0('Fig2-malignantCells-nc17604-MPscores-heamapgrps-',Sys.Date(),'.tsv'))

library(ggplot2)
vlndat1.long = tidyr::gather(vlndat1, key = MPs, value = MPscore, - k5_rename)
plots = ggplot(data = vlndat1.long,aes(x=`k5_rename`,y=`MPscore`)) +
  geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  k5_rename)) +
  geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
  xlab('') + ylab('score') +
  theme_classic() +
  ##scale_fill_manual(values = colors) + 
  facet_wrap(~MPs,scales = 'free') + 
  theme(axis.text.x = element_text(angle = 90))
ggsave('figure2-selectedMPs-byUnsupervisedGroup-nc17604.pdf',plots,width = 7,
       height = 12)
}

##figure 1k: malignant cell analysis of P14
if(TRUE) {
P14_data = read_tsv('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/figures/figures_data/epithelial/EpithelialUmap/Epithelial_normalCells/alveolar_airway_andAllfiguresworkdir/tumorSamples_clustering/persampleClustering/P14_perSampleUmapWithAnnotatoin2021-06-03.jmp')
##add code for KRAS mutation

##add code for KRAS signature score

##add code for  differentiation score

##monocle results
library(Seurat);library(monocle)
ClusterObj <- p14seurat
obj = ClusterObj
data <- Seurat::GetAssayData(obj, slot = 'counts') 
pd <- data.frame(data = obj@meta.data)
pd1 = pd
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- fData
metaData = pd1
expr_matrix = data
geneAnnotation = fData
pd <- new("AnnotatedDataFrame", data = metaData)
fd <- new("AnnotatedDataFrame", data = geneAnnotation)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 50))
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds,num_paths = 2)

##add other variables
maligdf = read_tsv('figure2-malignantCells-17216-allinone-table-2022-01-14.tsv')
pos = match(rownames(pData(cds)), maligdf$ID)

pData(cds)$CytoTRACEScore = maligdf$CytoTRACEScore[pos]
pData(cds)$aviv_cnvscore_hg38 = maligdf$aviv_cnvscore_hg38[pos]

cnvclst = read_tsv('Fig2-P14-cnv-refT-tumorCells-heatmap-cellmetadf.tsv') ##didn't do this update as it's time consuming
pos = match(rownames(pData(cds)), cnvclst$id)

pData(cds)$cnvGrp = paste0('cnvGrp',cnvclst$cnvHeatmapGrp[pos])
pdf(paste0('Fig2-P14-tumorCells-ntotalmalig17064-monocle2-plots-',Sys.Date(),'.pdf'),
    useDingbats = F)
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "CytoTRACEScore")
plot_cell_trajectory(cds, color_by = "aviv_cnvscore_hg38")
plot_cell_trajectory(cds, color_by = "cnvGrp")
dev.off()
saveRDS(cds, paste0('Fig2-P14-tumorCells-ntotalmalig17064-monocle2-obj-',Sys.Date(),'.rds'))
cds = readRDS('Fig2-P14-tumorCells-ntotalmalig17064-monocle2-obj-2022-01-14.rds')
mc2df = as.data.frame( t(cds@reducedDimS) ); mc2df$ID = rownames(mc2df)
pseudotimdf = pData(cds)
mc2df = cbind(mc2df, pseudotimdf[mc2df$ID,])
write_tsv(mc2df, paste0('Fig2-P14-tumorCells-ntotalmalig17064-monocle2-coord-',Sys.Date(),'.tsv'))

##output coordinate
cds1 <- reduceDimension(cds, max_components = 2, method = 'ICA')
cds1 <- orderCells(cds1,num_paths = 1)
pdf(paste0('Fig2-P14-tumorCells-ntotalmalig17064-monocle2-ICA-plots-',Sys.Date(),'.pdf'),
    useDingbats = F)
plot_cell_trajectory(cds1, color_by = "Pseudotime")
plot_cell_trajectory(cds1, color_by = "State")
plot_cell_trajectory(cds1, color_by = "CytoTRACEScore")
plot_cell_trajectory(cds1, color_by = "aviv_cnvscore_hg38")
plot_cell_trajectory(cds1, color_by = "cnvGrp")
dev.off()
}




