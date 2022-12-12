##figure 2a
if(TRUE) {
    write_tsv(perSampleFrac, paste0('Fig4PerPatientcompostionbarplot-P8P9hg38rRef-data-',Sys.Date(),'.tsv'))
##st test
ksRes = perSampleFrac %>%
group_by(ct1) %>%
    do({
        tryCatch(broom::tidy(kruskal.test(x = .$Frac ,g = .$Field)),
                 error = function(e) data.frame('error'))
    }) 
ksRes$fdr = p.adjust(ksRes$p.value,'BH')
write_tsv(ksRes, paste0('Fig4PerPatientcompostionbarplot-P8P9hg38rRef-data-ksTest-',Sys.Date(),'.tsv'))

##pairwise test
tRes = perSampleFrac %>%
  group_by(ct1) %>%
  do({
    tryCatch(broom::tidy(pairwise.wilcox.test(x = .$Frac ,g = .$Field)),
             error = function(e) data.frame(as.character(e)))
  }) 
tRes$fdr = p.adjust(tRes$p.value,'BH')
Param = perSampleFrac %>%
    group_by(ct1,Field) %>%
summarise(n = n(),
          med = median(Frac),
          mean = mean(Frac))
tRes1 = left_join(tRes,filter(Param,Field == 'tumor'),by = c('ct1'='ct1')) %>%
    left_join(filter(Param,Field == 'adjacent'),by = c('ct1'='ct1')) %>%
    left_join(filter(Param,Field == 'intermediate'),by = c('ct1'='ct1')) %>%
    left_join(filter(Param,Field == 'distant'),by = c('ct1'='ct1'))

write_tsv(tRes1, paste0('Fig4PerPatientcompostionbarplot-P8P9hg38rRef-data-wilcxTest-',Sys.Date(),'.tsv'))


##lined dot plot for fig4 
library(ggpubr)
ct1sig = filter(ksRes, fdr <1) %>% .$ct1 ##generated all
perSampleFrac$Field = factor(perSampleFrac$Field,levels = c(
  'tumor','adjacent','intermediate','distant'
))
sampleColor = c('tumor' = '#df5c5b',
                'distant' ='#7bbde6' ,
                'intermediate' = '#7a7979',
                'adjacent' = '#c0bfc0')
gp = ggpaired(filter(perSampleFrac,ct1 %in% ct1sig),
              x= 'Field',y = 'Frac',id = 'patient',
              color = 'Field',line.color = 'grey',line.size = .4,
              palette = sampleColor,
              scales = 'free',
              xlab = '',facet.by = 'ct1',
              ylab = 'Proportion (among epithelial cells)') + stat_compare_means(paired = T)
ggsave(paste0('Fig4SuppCellTypes-pairedboxplot-P8P9hg38rRef-',Sys.Date(),'.pdf'),
       gp,useDingbats = F,height = 7,width = 6)


#######################merged samples in each field######
##update: fraction plot (now in figure 3) figure 3 (fig3) fraction barplot no proliferating cells
ct %>%
  filter(!(ct1 %in% c('alveoli_Sprites','proliferating_alveolar',    
                      'proliferating_basal',   
                      'proliferating_club'))) %>%
  mutate(
         Field = as.factor(Field),
         ct1 = as.factor(ct1)) %>%
  group_by(Field, .drop = F) %>%
  summarise(ncellPerPatient = n()) -> ncellPerField
ct %>%
  filter(!(ct1 %in% c('alveoli_Sprites','proliferating_alveolar',    
                      'proliferating_basal',   
                      'proliferating_club'))) %>%
  mutate(
         Field = as.factor(Field),
         ct1 = as.factor(ct1)) %>%
  group_by(Field,ct1, .drop = F) %>%
  summarise(ncellPerctPerPatient = n()) -> ncellPerctPerField

inner_join(ncellPerField, ncellPerctPerField,
           by = c('Field' = 'Field')) %>%
  mutate(Frac = ncellPerctPerPatient/ncellPerPatient) -> perFieldFrac

perFieldFrac$Field = factor(perFieldFrac$Field,
                                levels = c(
                                  'tumor', 'adjacent','intermediate', 'distant'                                 
                                ))
perFieldFrac$ct1 = factor(perFieldFrac$ct1,
                          levels = c('alveoli_AVP',
                            'alveoli_AT1', 'alveoli_AT2',  'alveoli_Bronchio','airway_Basal', 'airway_Ciliated',## 'proliferating_alveolar',
                            'airway_ClubandSecretory', 'airway_Ionocyte', 'airway_Neuroendocrine', 'airway_Tufft' ##  'proliferating_basal', 'proliferating_club'
                            ))

colors = c(
  '#8da0cb',
  '#66c2a5',
  '#fc8d62',
  
  '#e78ac3','#e37e59', '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462',
  '#9fd1c7', '#bdbad7'
)
names(colors) = c(
  
  'alveoli_AVP','alveoli_AT1', 'alveoli_AT2',  'alveoli_Bronchio', 'proliferating_alveolar', 
  'airway_Basal', 'airway_Ciliated', 'airway_ClubandSecretory', 'airway_Ionocyte', 'airway_Neuroendocrine', 'airway_Tufft',
  'proliferating_basal', 'proliferating_club')
library(ggplot2)
##make a barplot
perFieldFrac = na.omit(perFieldFrac) ##remove P1 intermediate since we don't have it
gp = ggplot(perFieldFrac,aes(x=Field))
gp = gp + geom_bar(stat='identity',aes(y=`Frac`,fill = ct1))
gp = gp + scale_fill_manual(values = colors)
gp = gp + theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.text.x = element_text(angle=60,hjust=1,size=10),
                             axis.title = element_text(size = 10),
                             axis.line = element_line(colour = "black"),
                             legend.position='right'
)
gp = gp + ylim(0,1.02)
##to remove space between y axis and y=0
gp = gp + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) 
gp = gp + ylab('% (of epithelial cells)') + xlab('')
##ggsave(paste0('Fig4PerFieldCompostionBarplot-P8P9hg38-',Sys.Date(),'.pdf'),
##       gp,height = 6,width = 6)
ggsave(paste0('Fig2a-PerFieldCompostionBarplot-noProlif-P8P9hg38-',Sys.Date(),'.pdf'),
      gp,height = 6,width = 6)
}

##figure 2b
if(TRUE) {
    ##lined dot plot for fig4 
library(ggpubr)
ct1sig = filter(ksRes, fdr <1) %>% .$ct1 ##generated all
perSampleFrac$Field = factor(perSampleFrac$Field,levels = c(
  'tumor','adjacent','intermediate','distant'
))
sampleColor = c('tumor' = '#df5c5b',
                'distant' ='#7bbde6' ,
                'intermediate' = '#7a7979',
                'adjacent' = '#c0bfc0')
gp = ggpaired(filter(perSampleFrac,ct1 %in% ct1sig),
              x= 'Field',y = 'Frac',id = 'patient',
              color = 'Field',line.color = 'grey',line.size = .4,
              palette = sampleColor,
              scales = 'free',
              xlab = '',facet.by = 'ct1',
              ylab = 'Proportion (among epithelial cells)') + stat_compare_means(paired = T)
ggsave(paste0('Fig2B-pairedboxplot-P8P9hg38rRef-',Sys.Date(),'.pdf'),
       gp,useDingbats = F,height = 7,width = 6)
}

##figure 2c&d: monocle trajectory analysis
if(TRUE) {
    ##code was run using monocle with default parameters and DDRTree model was used for trejactory construction
    ## the inferred trajectory coordinates was exported and visualized using jmp.
}

## figure2e bubbleplot
if(TRUE) {
##bubbleplot
library(ggplot2)
data = read_tsv('avpobj_exp_bubbleData.tsv') ##check
colors = c('#fcfbfd', '#efedf5', '#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#4a1486')
colors1 = colorRampPalette(colors)(50)
barp <-ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
  scale_fill_gradientn(colours=colors1)+       ## color of the corresponding aes
  labs(x='',y='')+
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+             ## to tune the size of circles
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )

ggsave(paste0('Fig4-KRT8-vbubbleplot-',Sys.Date(),'.pdf'), barp,height = 3, width = 6)
}


## figure2g: KACs among epithelial cells
if(TRUE) {
    ##lined dot plot for fig4 
    perSampleFrac_KAC = read_tsv('KAC_frac.tsv')
library(ggpubr)
ct1sig = filter(ksRes, fdr <1) %>% .$ct1 ##generated all
perSampleFrac$Field = factor(perSampleFrac$Field,levels = c(
  'tumor','adjacent','intermediate','distant'
))
sampleColor = c('tumor' = '#df5c5b',
                'distant' ='#7bbde6' ,
                'intermediate' = '#7a7979',
                'adjacent' = '#c0bfc0')
gp = ggpaired(filter(perSampleFrac,ct1 %in% ct1sig),
              x= 'Field',y = 'Frac',id = 'patient',
              color = 'Field',line.color = 'grey',line.size = .4,
              palette = sampleColor,
              scales = 'free',
              xlab = '',facet.by = 'ct1',
              ylab = 'Proportion (among epithelial cells)') + stat_compare_means(paired = T)
ggsave(paste0('Fig2g-pairedboxplot-P8P9hg38rRef-',Sys.Date(),'.pdf'),
       gp,useDingbats = F,height = 7,width = 6)
}

##figure2h: galaxy plot showing density of AIC umap
if(TRUE) {
    ##fig4 galaxy plot of avp subclusters
umap = read_tsv('AllnormalEpithelial_p8p9hg38_avpobjnc10033UMAP_-2022-01-20.tsv')
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

g.overlay = ggplot(data = umap,
                   aes(x = UMAP_1, y = UMAP_2))
g.overlay = g.overlay + stat_density_2d(aes(fill = ..density..), geom = "raster",
                                        contour = F)
g.overlay = g.overlay + geom_point(color = 'white',size = .02)
g.overlay = g.overlay + facet_wrap(~Field,ncol = 2)
g.overlay = g.overlay + scale_fill_viridis(option="magma")
g.overlay = g.overlay + galaxyTheme_black()
ggsave(paste0('Fig4-avpcell-galaxyPlot-',Sys.Date(),'.png'),
       g.overlay)

umap$Field1 = umap$Field; umap$Field1[umap$Field %in% c('distant','intermediate','adjacent')] = 'NL'
g.overlay = ggplot(data = umap,
                   aes(x = UMAP_1, y = UMAP_2))
g.overlay = g.overlay + stat_density_2d(aes(fill = ..density..), geom = "raster",
                                        contour = F)
g.overlay = g.overlay + geom_point(color = 'white',size = .02)
g.overlay = g.overlay + facet_wrap(~Field1,ncol = 1)
g.overlay = g.overlay + scale_fill_viridis(option="magma")
g.overlay = g.overlay + galaxyTheme_black()
ggsave(paste0('Fig4-avpcell-galaxyPlot-tvsNL-',Sys.Date(),'.png'),
       g.overlay,width = 6,height = 12)
}

##figure2k: violinplots showing kras and alveolar signature score level
if(TRUE) {
    MPscore = read_tsv('AllnormalEpithelial_AlveolarUMAP-correctP8P9_nc187914-metaprogramScore-2022-01-24.tsv')
ct = read_tsv('normalEpithelialcells_celltype_table_update_P8P9hg38.tsv') ##update P8P9
ct = filter(ct, !(ct1 %in% c('airway_Airway_sprite','alveoli_Sprites')))
gc()
krt8cells = filter(read_tsv('AllnormalEpithelial_p8p9hg38_avpobjnc10033UMAP_-2022-01-20.tsv')) %>%
  filter(snn_default_param_res.0.1 == 1)
oldkrt8cells = read_tsv('krt8_nc1754_cellID.txt') %>% mutate(ID = gsub('-1$','',gsub('^_','',ID)))
ct$krt8 = rep('Other',nrow(ct))
ct$krt8[ct$ID %in% krt8cells$ID] = 'KRT8'
ct$ct1[ct$krt8 == 'KRT8'] = 'KRT8_AVP'
MPtable = inner_join(ct,MPscore, by = c('ID'='cellID')) %>%
  filter(ct1 %in% c('alveoli_AT1', 'alveoli_AT2','alveoli_AVP','KRT8_AVP'))

colors = c(
  '#89B4E0','#E27525',
  '#66c2a5',
  '#fc8d62'
)
names(colors) = c(
  
  'alveoli_AVP','KRT8_AVP','alveoli_AT1', 'alveoli_AT2')
library(ggplot2)
MPtable.long = tidyr::gather(select(MPtable,ID, ct1,Field.x, patient.x, one_of(paste0('MP_',1:41))), key = MPs,
                             value = score, -ct1, -ID,-Field.x,-patient.x)
plots = ggplot(data = MPtable.long,aes(x=`ct1`,y=`score`)) +
  geom_violin(scale = 'width', alpha = .8,lwd = 0.1,aes(fill =  ct1)) +
  geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
  xlab('') + ylab('') +
  theme_classic() +
  scale_fill_manual(values = colors) + facet_wrap(~MPs,ncol = 4, scales = 'free') + 
  theme(axis.text.x = element_text(angle = 90))
ggsave('figure2-k-MP31.pdf',plots,width = 14,
       height = 24)

##KRT8 signature in tumors vs. kras signature score in tumor samples
##alternative: scatterplots
plotx<-ggplot(MPtable, aes(x = ownKRAS,y = `krt8_tp100`)) +        ## global aes
  geom_point(aes(color = ct1)) +
  labs(x='',y='')+
  scale_color_manual(values = colors)+facet_wrap(~ct1, ncol = 2, scales = 'free')+
  theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")

ggsave(paste0('figure2-k-krasSig-vs-krt8Sig-updated-scatter-',Sys.Date(),'.pdf'),plotx,
       useDingbats = F)
}


##figure2l: kac signature analysis in KM LUAD vs EM LUAD
if(TRUE) {
##quantify KAC signature score in tumor cells
data = read_tsv('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/figures/figures_data/epithelial/EpithelialUmap/Epithelial_normalCells/alveolar_airway_andAllfiguresworkdir/figure4-KRT8sig-updatekrt8vsalveolar-malignantCells-data-2022-02-11.tsv')
maligtable = read_tsv('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/nature_revision/Reviewer2/mutationScreening/prepare_bam_files/malignantCells_nc17064.tsv')
data1 = inner_join(data,maligtable,by = c( 'ID' = 'ID'))
plotx<-ggplot(data1, aes(x = patient,y = `KRT8top100`)) +        ## global aes
    geom_violin(scale = 'width', alpha = .8,lwd = 0.1) +
    geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
    labs(x='',y='')+
    scale_color_manual(values = colors)+
    theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")

ggsave(paste0('R2Q4-malignant-KACSig-betweenKRASvsEGFR-samplelvl-',Sys.Date(),'.pdf'),plotx,
       useDingbats = F)

plotx<-ggplot(data1, aes(x = mutation,y = `KRT8top100`)) +        ## global aes
    geom_violin(scale = 'width', alpha = .8,lwd = 0.1) +
    geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
    labs(x='',y='')+
    scale_color_manual(values = colors)+
    theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")

ggsave(paste0('R2Q4-malignant-KACSig-betweenKRASvsEGFR-grplvl-',Sys.Date(),'.pdf'),plotx,
       useDingbats = F)

 mutationTab = data.frame(patient = paste0('P',1:16),mutation=NA)
mutationTab$mutation = rep('Other', nrow(mutationTab))
mutationTab$mutation[mutationTab$patient %in% c('P2','P10','P14')] = 'KRAS'
mutationTab$mutation[mutationTab$patient %in% c('P1','P6','P7','P8','P15')] = 'EGFR'
mutationTab$mutation[mutationTab$patient %in% c('P4',"P9")] = 'MET'
MPtable = inner_join(MPtable, mutationTab, by = c('patient.x' = 'patient'))
gpdat = filter(MPtable,mutation %in% c('KRAS','EGFR') & ct1 == 'KRT8_AVP')
gpdat$patient.x = factor(gpdat$patient.x,levels = c('P2', 'P10', 'P14', 'P6', 'P7', 'P15', 'P1', 'P8'))
plotx<-ggplot(gpdat, aes(x = patient.x,y = `krt8_tp100`)) +        ## global aes
    geom_violin(scale = 'width', alpha = .8,lwd = 0.1) +
    geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
    labs(x='',y='')+
    scale_color_manual(values = colors)+
    theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")

ggsave(paste0('R2Q4-KACSig-betweenKRASvsEGFR-',Sys.Date(),'.pdf'),plotx,
       useDingbats = F)
wilcox.test(filter(gpdat,mutation == 'KRAS') %>% .$krt8_tp100,y = filter(gpdat,mutation == 'EGFR') %>% .$krt8_tp100)
plotx<-ggplot(gpdat, aes(x = mutation,y = `krt8_tp100`)) +        ## global aes
    geom_violin(scale = 'width', alpha = .8,lwd = 0.1) +
    geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
    labs(x='',y='')+
    scale_color_manual(values = colors)+
    theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")

ggsave(paste0('R2Q4-KACSig-betweenKRASvsEGFR-grouplevel-',Sys.Date(),'.pdf'),plotx,
       useDingbats = F)
}

##figure2m: TCGA and pre-neoplasia validation
if(TRUE) {
##quantify KAC signature score in tumor cells
data = read_tsv('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/figures/figures_data/epithelial/EpithelialUmap/Epithelial_normalCells/alveolar_airway_andAllfiguresworkdir/figure4-KRT8sig-updatekrt8vsalveolar-malignantCells-data-2022-02-11.tsv')
maligtable = read_tsv('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/nature_revision/Reviewer2/mutationScreening/prepare_bam_files/malignantCells_nc17064.tsv')
data1 = inner_join(data,maligtable,by = c( 'ID' = 'ID'))
plotx<-ggplot(data1, aes(x = patient,y = `KRT8top100`)) +        ## global aes
    geom_violin(scale = 'width', alpha = .8,lwd = 0.1) +
    geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
    labs(x='',y='')+
    scale_color_manual(values = colors)+
    theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")

ggsave(paste0('R2Q4-malignant-KACSig-betweenKRASvsEGFR-samplelvl-',Sys.Date(),'.pdf'),plotx,
       useDingbats = F)

plotx<-ggplot(data1, aes(x = mutation,y = `KRT8top100`)) +        ## global aes
    geom_violin(scale = 'width', alpha = .8,lwd = 0.1) +
    geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
    labs(x='',y='')+
    scale_color_manual(values = colors)+
    theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")

ggsave(paste0('R2Q4-malignant-KACSig-betweenKRASvsEGFR-grplvl-',Sys.Date(),'.pdf'),plotx,
       useDingbats = F)

 mutationTab = data.frame(patient = paste0('P',1:16),mutation=NA)
mutationTab$mutation = rep('Other', nrow(mutationTab))
mutationTab$mutation[mutationTab$patient %in% c('P2','P10','P14')] = 'KRAS'
mutationTab$mutation[mutationTab$patient %in% c('P1','P6','P7','P8','P15')] = 'EGFR'
mutationTab$mutation[mutationTab$patient %in% c('P4',"P9")] = 'MET'
MPtable = inner_join(MPtable, mutationTab, by = c('patient.x' = 'patient'))
gpdat = filter(MPtable,mutation %in% c('KRAS','EGFR') & ct1 == 'KRT8_AVP')
gpdat$patient.x = factor(gpdat$patient.x,levels = c('P2', 'P10', 'P14', 'P6', 'P7', 'P15', 'P1', 'P8'))
plotx<-ggplot(gpdat, aes(x = patient.x,y = `krt8_tp100`)) +        ## global aes
    geom_violin(scale = 'width', alpha = .8,lwd = 0.1) +
    geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
    labs(x='',y='')+
    scale_color_manual(values = colors)+
    theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")

ggsave(paste0('R2Q4-KACSig-betweenKRASvsEGFR-',Sys.Date(),'.pdf'),plotx,
       useDingbats = F)
wilcox.test(filter(gpdat,mutation == 'KRAS') %>% .$krt8_tp100,y = filter(gpdat,mutation == 'EGFR') %>% .$krt8_tp100)
plotx<-ggplot(gpdat, aes(x = mutation,y = `krt8_tp100`)) +        ## global aes
    geom_violin(scale = 'width', alpha = .8,lwd = 0.1) +
    geom_boxplot(color = 'grey60',width = .2,outlier.shape= NA)  +
    labs(x='',y='')+
    scale_color_manual(values = colors)+
    theme_classic() +geom_smooth(method=lm, se=FALSE)+ ggpubr::stat_cor(method = "pearson")

ggsave(paste0('R2Q4-KACSig-betweenKRASvsEGFR-grouplevel-',Sys.Date(),'.pdf'),plotx,
       useDingbats = F)
}

##figure2n: 
if(TRUE) {
    ## done by jmp: add the data table here
}

##figure2o: multivariate analysis 
if(TRUE) {
    ##compare TCGA KRAS MUT vs. WT the KAC signature score and other AIC signature score
done using jmp file:/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/nature_revision/Reviewer2/R2Q2-Fig2I TCGA KRAS KAC signature score.jmp
##survival analysis
survdata = read_tsv('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/figures/figures_data/epithelial/EpithelialUmap/Epithelial_normalCells/alveolar_airway_andAllfiguresworkdir/TCGA-LUAD-survivalData-grpbyKRT8SigSig-2022-02-11.tsv')
survdata1 = survdata
fit = survfit(Surv(OS.time_6yrs,OS_6yrs) ~ KRT8SigSigQ2, data = survdata1)
mypalette = c('blue','brown');names(mypalette) = c('Q2=1','Q2=2')
ttpal = mypalette
ggsurvplot(
    fit, 
    data = survdata1, 
    ## xlim = c(0,100),
    size = 1,                 # change line size
    conf.int = F,          # Add confidence interval
    pval = TRUE,              # Add p-value
    ## palette = ttpal,
    risk.table = T,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    xlab = 'Time in months',  
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw()      # Change ggplot2 theme
) -> g
##pdf(paste0('TCGA-LUAD-KRT8SigSigQ4-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6)
pdf(paste0('TCGA-LUAD-KRT8SigSigQ2-6yrs-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6) ##really nice survival results
print(g,newpage = F)
dev.off()

survdata1.nostageIV = filter(survdata,stage != 'Stage IV' )
fit = survfit(Surv(OS.time_6yrs,OS_6yrs) ~ KRT8SigSigQ2, data = survdata1.nostageIV)
mypalette = c('blue','brown');names(mypalette) = c('Q2=1','Q2=2')
ttpal = mypalette
ggsurvplot(
    fit, 
    data = survdata1, 
    ## xlim = c(0,100),
    size = 1,                 # change line size
    conf.int = F,          # Add confidence interval
    pval = TRUE,              # Add p-value
    ## palette = ttpal,
    risk.table = T,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    xlab = 'Time in months',  
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw()      # Change ggplot2 theme
) -> g
##pdf(paste0('TCGA-LUAD-KRT8SigSigQ4-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6)
pdf(paste0('TCGA-LUAD-nostageIV-KRT8SigSigQ2-6yrs-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6) ##really nice survival results
print(g,newpage = F)
dev.off()

survdata1.stageI = filter(survdata,stage == 'Stage I' )
fit = survfit(Surv(OS.time_6yrs,OS_6yrs) ~ KRT8SigSigQ2, data = survdata1.stageI)
mypalette = c('blue','brown');names(mypalette) = c('Q2=1','Q2=2')
ttpal = mypalette
ggsurvplot(
    fit, 
    data = survdata1, 
    ## xlim = c(0,100),
    size = 1,                 # change line size
    conf.int = F,          # Add confidence interval
    pval = TRUE,              # Add p-value
    ## palette = ttpal,
    risk.table = T,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    xlab = 'Time in months',  
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw()      # Change ggplot2 theme
) -> g
##pdf(paste0('TCGA-LUAD-KRT8SigSigQ4-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6)
pdf(paste0('TCGA-LUAD-stageI-KRT8SigSigQ2-6yrs-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6) ##really nice survival results
print(g,newpage = F)
dev.off()

##group by other AIC signature score
kacTCGAsignature = read_tsv('TCGA-KRT8signature-mcpcounter-nsample563-otherAICSignature-updated-2022-06-27.tsv')
df = read_excel('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/mouse_data_1006/clusterEachLineage4th/Epithelial/KO samples only/TCGA_LUAD/TCGALUAD_mutation_Avp_N498.xlsx')
df1 = df
sigScore = kacTCGAsignature %>% filter(ID %in% df1$ID_ExpressionData) %>%
    mutate(ShortID = substr(ID, 1,12),typeCode = substr(ID,14,15)) 
survdata1 = inner_join(sigScore, survdata1, by = c('ShortID' = 'ShortID'))

survdata1 = survdata1 %>%
    mutate(
        OtherAICQ4 = ntile(OtherAIC_tp100, 4),
        OtherAICQ2 = ntile(OtherAIC_tp100,2))

fit = survfit(Surv(OS.time_6yrs,OS_6yrs) ~ OtherAICQ2, data = survdata1)
mypalette = c('blue','brown');names(mypalette) = c('Q2=1','Q2=2')
ttpal = mypalette
ggsurvplot(
    fit, 
    data = survdata1, 
    ## xlim = c(0,100),
    size = 1,                 # change line size
    conf.int = F,          # Add confidence interval
    pval = TRUE,              # Add p-value
    ## palette = ttpal,
    risk.table = T,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    xlab = 'Time in months',  
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw()      # Change ggplot2 theme
) -> g
##pdf(paste0('TCGA-LUAD-KRT8SigSigQ4-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6)
pdf(paste0('TCGA-LUAD-OtherAICSigSigQ2-6yrs-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6) ##really nice survival results
print(g,newpage = F)
dev.off()

survdata1.nostageIV = filter(survdata1,stage != 'Stage IV' )
fit = survfit(Surv(OS.time_6yrs,OS_6yrs) ~ OtherAICQ2, data = survdata1.nostageIV)
mypalette = c('blue','brown');names(mypalette) = c('Q2=1','Q2=2')
ttpal = mypalette
ggsurvplot(
    fit, 
    data = survdata1, 
    ## xlim = c(0,100),
    size = 1,                 # change line size
    conf.int = F,          # Add confidence interval
    pval = TRUE,              # Add p-value
    ## palette = ttpal,
    risk.table = T,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    xlab = 'Time in months',  
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw()      # Change ggplot2 theme
) -> g
##pdf(paste0('TCGA-LUAD-KRT8SigSigQ4-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6)
pdf(paste0('TCGA-LUAD-nostageIV-OtherAICSigSigQ2-6yrs-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6) ##really nice survival results
print(g,newpage = F)
dev.off()

survdata1.stageI = filter(survdata1,stage == 'Stage I' )
fit = survfit(Surv(OS.time_6yrs,OS_6yrs) ~ OtherAICQ2, data = survdata1.stageI)
mypalette = c('blue','brown');names(mypalette) = c('Q2=1','Q2=2')
ttpal = mypalette
ggsurvplot(
    fit, 
    data = survdata1, 
    ## xlim = c(0,100),
    size = 1,                 # change line size
    conf.int = F,          # Add confidence interval
    pval = TRUE,              # Add p-value
    ## palette = ttpal,
    risk.table = T,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    xlab = 'Time in months',  
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw()      # Change ggplot2 theme
) -> g
##pdf(paste0('TCGA-LUAD-KRT8SigSigQ4-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6)
pdf(paste0('TCGA-LUAD-stageI-OtherAICSigSigQ2-6yrs-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6) ##really nice survival results
print(g,newpage = F)
dev.off()

##group by KAC signatuer high vs othe AIC signature hig
write_tsv(survdata1, paste0('TCGA-LUAD-OtherAICsignatureScore-',Sys.Date(),'.tsv'))
##survival between the combination of KRT8sig and OtherAIC sig
cor.test(survdata1$OtherAIC_tp100,survdata$krt8_tp100) ##correlation score: .29, p < 0.001
group = paste0(survdata1$KRT8SigSigQ2,'-',survdata1$OtherAICQ2)
survdata1$grp1 = group

fit = survfit(Surv(OS.time_6yrs,OS_6yrs) ~ grp1, data = survdata1)
mypalette = c('blue','brown');names(mypalette) = c('Q2=1','Q2=2')
ttpal = mypalette
ggsurvplot(
    fit, 
    data = survdata1, 
    ## xlim = c(0,100),
    size = 1,                 # change line size
    conf.int = F,          # Add confidence interval
    pval = TRUE,              # Add p-value
    ## palette = ttpal,
    risk.table = T,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    xlab = 'Time in months',  
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw()      # Change ggplot2 theme
) -> g
##pdf(paste0('TCGA-LUAD-KRT8SigSigQ4-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6)
pdf(paste0('TCGA-LUAD-byKACOtherAICCombination-6yrs-OSKMplot-',Sys.Date(),'.pdf'),height = 6,width = 6) ##really nice survival results
print(g,newpage = F)
dev.off()

##run a multivariate analysis
library(survival)
survdata1$age1 = rep(NA, nrow(survdata1));
survdata1$age1[survdata1$age > 65] = 'age2' ##mean age == 65.44
survdata1$age1[survdata1$age < 65] = 'age1'
survdata1$stage1 = rep(NA, nrow(survdata1));
survdata1$stage1[survdata1$stage != 'Stage I' ] = 'stage2' 
survdata1$stage1[survdata1$stage  == 'Stage I'] = 'stage1'

survdata1$stage2 = rep(NA, nrow(survdata1));
survdata1$stage2[survdata1$stage != 'Stage IV' ] = 'stage123' 
survdata1$stage2[survdata1$stage  == 'Stage IV'] = 'stage4'

multivarcox=broom::tidy(coxph(Surv(OS.time_6yrs,OS_6yrs) ~ KRT8SigSigQ2 + OtherAICQ2 + stage1 + age1,data = survdata1 ),exponentiate = T,conf.int = .95)
library(rmeta)
max.conf.int = 10
min.conf.int = -10
os.hazard.data1 = multivarcox %>%
    na.omit() %>%
    dplyr::rename('HR' = estimate) %>%
    mutate(q = p.adjust(p.value,'BH')) %>%
    arrange(q)
os.hazard.data1$comp = 'tbd'
tabletext = cbind(c('term',as.character(os.hazard.data1$term)),
                  c('Comparison',os.hazard.data1$comp),
                  c('HR',signif(os.hazard.data1$HR,2)),
                  c('P',signif(os.hazard.data1$q,2)),
                  c('cox.p.adjust',signif(os.hazard.data1$q,2)))
library(rmeta)
library(metafor)
pdf(paste0('R2-Q2-otheAIC-multivariate-OR-methoda-',Sys.Date(),'.pdf'),width = 12)
forestplot(tabletext,
           mean = c(NA,os.hazard.data1$HR),
           lower = c(NA,os.hazard.data1$conf.low),
           upper = c(NA,os.hazard.data1$conf.high),
           zero=1,is.summary=F,
           clip=c(0,5),
           xlog=F,xticks = c(0.5, 1, 2,3,4),
           col=meta.colors(box="royalblue",line="darkblue", summary="royalblue")
           )
dev.off()
}

