##figure 2a
if(TRUE) {
    colors = c(
        '#e05d5c',
        '#66c4d1',
        '#85af5e','#e78ac3',
        '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462',
        '#9fd1c7', '#bdbad7'
    )
    names(colors) = c(
        'AIC','AT1', 'AT2',  'SDP', 
        'Basal', 'Ciliated', 'Club', 'Ionocyte', 'NE', 'Tuft')
    perFieldFrac = readRDS('Input_Data/Fig2a_FracBarplotData.rds')
    perFieldFrac$celltype = factor(perFieldFrac$celltype,
                                   levels = c(
                                       'AIC','AT1','AT2','SDP','Basal',
                                       'Ciliated','Club', 'Ionocyte', 'NE', 'Tuft'
                                   ))
    perFieldFrac$Field = factor(perFieldFrac$Field,
                                levels = rev(levels(perFieldFrac$Field)))
    library(ggplot2)
    ##make a barplot
    Fig2a_FracBarplot = ggplot(perFieldFrac,aes(x=Field))
    Fig2a_FracBarplot = Fig2a_FracBarplot + geom_bar(stat='identity',aes(y=`Frac`,fill = celltype))
    Fig2a_FracBarplot = Fig2a_FracBarplot + scale_fill_manual(values = colors)
    Fig2a_FracBarplot = Fig2a_FracBarplot + theme_bw() + theme(panel.border = element_blank(),
                                                               panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank(),
                                                               axis.text.x = element_text(angle=60,hjust=1,size=10),
                                                               axis.title = element_text(size = 10),
                                                               axis.line = element_line(colour = "black"),
                                                               legend.position='right'
                                                               )
    Fig2a_FracBarplot = Fig2a_FracBarplot + ylim(0,1.02)
    ##to remove space between y axis and y=0
    Fig2a_FracBarplot = Fig2a_FracBarplot + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) 
    Fig2a_FracBarplot = Fig2a_FracBarplot + ylab('% (of epithelial cells)') + xlab('')
    print(Fig2a_FracBarplot)
}

##figure 2b
if(TRUE) {
    ##dot plot for fraction changes along sample field
    library(ggpubr)
    sampleColor = c('tumor' = '#df5c5b',
                    'distant' ='#7bbde6' ,
                    'intermediate' = '#7a7979',
                    'adjacent' = '#c0bfc0')
    perSampleFrac = readRDS('Input_Data/Fig2b_fractionData-v1.rds')
    perSampleFrac$Field = factor(perSampleFrac$Field,
                                 levels = rev(c('tumor','adjacent','intermediate','distant')))
    perSampleFrac$celltype = factor(perSampleFrac$celltype,
                                    levels = c('AT2','AIC'))
    Fig2b = ggpaired(perSampleFrac,
                     x= 'Field',y = 'Frac',id = 'patient',
                     color = 'Field',line.color = 'grey',line.size = 0,
                     palette = sampleColor,
                     scales = 'free',
                     xlab = '',facet.by = 'celltype',
                     ylab = 'Proportion (among epithelial cells)') + stat_compare_means(paired = T)
    print(Fig2b)
}

##fig2c,fig2d monocle2 results and fig2e umap on the left were visualized using JMP Pro v15.

## fig2e bubbleplot
if(TRUE) {
    ##bubbleplot
    library(ggplot2)
    data = readRDS('Input_Data/Fig2e_bubblePlotData.rds')
    data$features.plot = factor(data$features.plot,
                                levels = c('CDKN1A', 'CLDN4', 'PLAUR', 'CDKN2A', 'KRT8'))
    data$id1 = rep(NA,nrow(data))
    data$id1[data$id == 0] = 'Other AICs'
    data$id1[data$id == 1] = 'KAC'
    colors = c('#fcfbfd', '#efedf5', '#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#4a1486')
    colors1 = colorRampPalette(colors)(50)
    Fig2e_bubblePlot <-ggplot(data, aes(y = features.plot,x = id1)) +        ## global aes
        geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
        scale_fill_gradientn(colours=colors1,
                             limits = c(-1,1))+       ## color of the corresponding aes
        labs(x='',y='')+
        scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+             ## to tune the size of circles
        theme(panel.grid.major = element_line(colour = "grey90",size=0.2), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )
    print(Fig2e_bubblePlot)
}

## fig2f was generated by importing the cytoTRACE prediction scores into JMP Pro v15 as is described in methods part.

## figure2g: KACs among epithelial cells
if(TRUE) {
    ##fraction of KACs in each sample
    Frac_KAC = readRDS('Input_Data/Fig2g_KACFrac_data.rds')
    library(ggpubr)
    Frac_KAC$Field = factor(Frac_KAC$Field,
                            levels = rev(c('tumor','adjacent','intermediate','distant')))
    sampleColor = c('tumor' = '#df5c5b',
                    'distant' ='#7bbde6' ,
                    'intermediate' = '#7a7979',
                    'adjacent' = '#c0bfc0')
    Fig2g = ggpaired(Frac_KAC,
                     x= 'Field',y = 'Frac',id = 'patient',
                     color = 'Field',line.color = 'grey',line.size = 0,
                     palette = sampleColor,
                     scales = 'free',
                     xlab = '',
                     ylab = 'Proportion (among epithelial cells)') + stat_compare_means(paired = T)
    print(Fig2g)
}

##figure2h: galaxy plot showing cell density of KAC in UMAP separating LUAD and normal lung tissue
if(TRUE) {
    galaxyTheme_black = function(base_size = 12, base_family = "") {
        theme_grey(base_size = base_size, base_family = base_family) %+replace%
            theme(
                                        # Specify axis options
                axis.line = element_blank(),  
                axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
                axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
                axis.ticks = element_line(color = "white", size  =  0.2),  
                axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
                axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
                axis.ticks.length = unit(0.3, "lines"),   
                                        # Specify legend options
                legend.background = element_rect(color = NA, fill = "black"),  
                legend.key = element_rect(color = "white",  fill = "black"),  
                legend.key.size = unit(1.2, "lines"),  
                legend.key.height = NULL,  
                legend.key.width = NULL,      
                legend.text = element_text(size = base_size*0.8, color = "white"),  
                legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
                legend.position = "none",  
                legend.text.align = NULL,  
                legend.title.align = NULL,  
                legend.direction = "vertical",  
                legend.box = NULL, 
                                        # Specify panel options
                panel.background = element_rect(fill = "black", color  =  NA),  
                panel.border = element_rect(fill = NA, color = "white"),  
                ##panel.grid.major = element_line(color = "grey35"),  
                panel.grid.major = element_blank(),  
                ##panel.grid.minor = element_line(color = "grey20"),  
                panel.grid.minor = element_blank(),  
                panel.spacing = unit(0.5, "lines"),   
                                        # Specify facetting options
                strip.background = element_rect(fill = "grey30", color = "grey10"),  
                strip.text.x = element_text(size = base_size*0.8, color = "white"),  
                strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
                                        # Specify plot options
                plot.background = element_rect(color = "black", fill = "black"),  
                plot.title = element_text(size = base_size*1.2, color = "white"),  
                plot.margin = unit(rep(1, 4), "lines")
                
            )
        
    }
    library(ggplot2);library(viridis)
    library(ggpubr)
    data = readRDS('Input_Data/Fig2h_KAC_UmapData.rds')
    data$Field1 = factor(data$Field1,levels = c('tumor','NL'))
    Fig2h_Umap = ggplot(data = data,
                        aes(x = UMAP_1, y = UMAP_2))
    Fig2h_Umap = Fig2h_Umap + stat_density_2d(aes(fill = ..density..), geom = "raster",
                                              contour = F)
    Fig2h_Umap = Fig2h_Umap + geom_point(color = 'white',size = .02)
    Fig2h_Umap = Fig2h_Umap + facet_wrap(~Field1,ncol = 2)
    Fig2h_Umap = Fig2h_Umap + scale_fill_viridis(option="magma")
    Fig2h_Umap = Fig2h_Umap + galaxyTheme_black()
    print(Fig2h_Umap) ##suggest to save in a figure with large width and height, cause the cell dots plotted in a small figure affect the visualization of the actual density
}
##figure2h barplot fraction of LUAD coming cells was generated using JMP Pro v15

##figure2j CNV score spatial plot: need ST data share permission, don't run, wait for data release
if(FALSE) {
    library(patchwork)
    library(cowplot)
    library(infercnv)
    library(data.table)
    library(dplyr)
    library(Seurat)
    library(ggplot2)
    ST_LUAD = readRDS('Input_Data/P14_LUAD_STobj.rds')
    library("ggplot2")
    pdf("Fig2j_LUAD_ST_cnvscore.pdf",width = 7*4)
    SpatialFeaturePlot(ST_LUAD, features = "CNV") + theme(legend.position = "right")
    dev.off()
}

##figure2k: violinplots and scatterplots showing correlation between kras signature score
## and alveolar signature score level in alveolar cells
if(TRUE) {
    colors = c(
        '#89B4E0','#E27525',
        '#66c2a5',
        '#91b566'
    )
    names(colors) = c('Other AICs','KAC','AT1', 'AT2')
    library(ggplot2)

    ##KRT8 signature in tumors vs. kras signature score in tumor samples scatterplots
    Fig2k_scatterData = readRDS('Input_Data/Fig2k_scatterData.rds')
    Fig2k_scatterData$celltype = gsub('alveoli_','',Fig2k_scatterData$ct1)
    Fig2k_scatterData$celltype[Fig2k_scatterData$celltype=='KRT8_AVP'] = 'KAC'
    Fig2k_scatterData$celltype[Fig2k_scatterData$celltype=='AVP'] = 'Other AICs'

    Fig2k_KAC_KRAS_scatterPlot <-ggplot(dplyr::filter(Fig2k_scatterData,celltype == 'KAC'),
                                        aes(x = ownKRAS,y = `krt8_tp100`)) +        ## global aes
        geom_point(color = 'grey60') +
        labs(x='',y='')+
        theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")
    print(Fig2k_KAC_KRAS_scatterPlot)

    Fig2kData = readRDS('Input_Data/Fig2k_VlnPlotData-v1.rds')
    Fig2kData$celltype = factor(Fig2kData$celltype, levels = c('AT1', 'KAC','Other AICs','AT2'))
    Fig2k_MP31_alveolarSig = ggplot(data = Fig2kData,aes(x=`celltype`,y=`score`)) +
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  celltype)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        xlab('') + ylab('') +
        theme_classic() +
        scale_fill_manual(values = colors) + facet_wrap(~MPs,ncol = 4, scales = 'free') + 
        theme(axis.text.x = element_text(angle = 90))
    print(Fig2k_MP31_alveolarSig)

    Fig2k_KAC_MP31Alveolar_scatterPlot <-ggplot(dplyr::filter(Fig2k_scatterData,celltype == 'KAC'),
                                                aes(x = AT2_altlas,y = `krt8_tp100`)) +        ## global aes
        geom_point(color = 'grey60') +
        labs(x='',y='')+
        theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")
    print(Fig2k_KAC_MP31Alveolar_scatterPlot)

}

##figure2l: kac signature analysis in KM LUAD vs EM LUAD
if(TRUE) {
    colors = c('EGFR' = '#325fa9','KRAS' = '#9d202d')
    Fig2l_KACSigInKACsData = readRDS('Input_Data/Fig2l_KACinKACs_data.rds')
    Fig2l_KACSigInKACsVln<-ggplot(Fig2l_KACSigInKACsData, aes(x = mutation,y = `KACSig`)) + ## global aes
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill = mutation)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        labs(x='',y='')+
        scale_fill_manual(values = colors)
    ##showing KAC signature score in KACs
    print(Fig2l_KACSigInKACsVln)

    ##quantify KAC signature score in tumor cells
    data1 = readRDS('Input_Data/Fig2l_KACinMalig_data.rds')
    Fig2l_KACSigInMaligVln<-ggplot(data1, aes(x = mutation,y = `KACSig`)) +        ## global aes
        geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill = mutation)) +
        geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
        labs(x='',y='')+
        scale_fill_manual(values = colors)
    print(Fig2l_KACSigInMaligVln)
}

##figure2m: TCGA and pre-neoplasia validation
if(TRUE) {
    ##TCGA LUAD cohort
    tcga_data = readRDS('Input_Data/Fig2m_TCGAcohort.rds')
    Fig2m_KACsig_inTCGA = ggpaired(tcga_data,
                                   x= 'Sample',y = 'value',id = 'ShortID',
                                   color = 'Sample',line.color = 'grey',line.size = .4,
                                   paletter = c('LUAD' = '#e3615e',
                                                'Normal' = '#7dbee6'
                                                ),
                                   xlab = '',
                                   ylab = '') +
        scale_color_manual(values  = c('LUAD' = '#e3615e',
                                       'Normal' = '#7dbee6'
                                       ))
    print(Fig2m_KACsig_inTCGA)

    ##preneoplasia cohort boxplot of KAC signature
    PreNeoData = readRDS('Input_Data/Fig2m_PreNeoCohort.rds')
    Fig2m_KACsig_inPreNeoplasia = ggpaired(PreNeoData,
                                           x= 'Group',y = 'sigscore',id = 'group1',
                                           color = 'Group',line.color = 'grey',line.size = .4,
                                           paletter = c('LUAD' = '#e3615e',
                                                        'AAH' = '#FEcfc8',
                                                        'NL' = '#7dbee6'
                                                        ),
                                           xlab = '',
                                           ylab = 'signature score') +
        scale_color_manual(values  = c('LUAD' = '#e3615e',
                                       'AAH' = '#FEcfc8',
                                       'NL' = '#7dbee6'
                                       ))
    print(Fig2m_KACsig_inPreNeoplasia)
}


##figure2o: multivariate analysis of KAC signature/Other AIC signature/disease stage/Age in TCGA LUAD cohort
if(TRUE) {
    ##run a multivariate analysis
    library(survival)
    survdata1 = readRDS('Input_Data/Fig2o_tcga_survivalData.rds')
    multivarcox=broom::tidy(coxph(Surv(OS.time_6yrs,OS_6yrs) ~ KAC + OtherAIC + Stage + Age,data = survdata1 ),exponentiate = T,conf.int = .95)
    library(rmeta)
    max.conf.int = 10
    min.conf.int = -10
    os.hazard.data1 = multivarcox %>%
        na.omit() %>%
        dplyr::rename('HR' = estimate) %>%
        mutate(q = p.adjust(p.value,'BH')) %>%
        dplyr::arrange(q)
    os.hazard.data1$comp = 'tbd'
    tabletext = cbind(c('term',as.character(os.hazard.data1$term)),
                      c('Comparison',os.hazard.data1$comp),
                      c('HR',signif(os.hazard.data1$HR,2)),
                      c('P',signif(os.hazard.data1$q,2)),
                      c('cox.p.adjust',signif(os.hazard.data1$q,2)))
    tabletext = tabletext[,-c(2,4)]
    tabletext[,1] = c('term','Stage',
                      'KAC',
                      'Age',
                      'Other AICs')

    library(rmeta)
    library(metafor)
    forestplot(tabletext,
               mean = c(NA,os.hazard.data1$HR),
               lower = c(NA,os.hazard.data1$conf.low),
               upper = c(NA,os.hazard.data1$conf.high),
               zero=1,is.summary=F,
               clip=c(0,5),
               xlog=F,xticks = c(0.5, 1, 2,3,4),
               col=meta.colors(box="royalblue",line="darkblue", summary="royalblue")
               )
}

