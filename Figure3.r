##figure3 b
if(TRUE) {
    library(ggplot2)
    fractionData = readr::read_tsv('Input_Data/Fig3_FractionChange_dataInput.txt')
    fractionData$treatment = factor(fractionData$treatment,
                                    levels = c('Saline','NNK'))

fractionData$timepoint = factor(fractionData$timepoint,
                                    levels = c('EOE','7 Mo'))

    fig3b_left = ggplot(data = dplyr::filter(fractionData,celltype == 'Malignant'),
                        aes(x = `treatment`,y=Fraction))+
    geom_boxplot( width = .3,
                 alpha = .3,outlier.shape= NA)  +
        geom_point() +
    xlab('') +
    theme_classic() + facet_wrap(~timepoint) + theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
    print(fig3b_left)


fig3b_right = ggplot(data = dplyr::filter(fractionData,celltype == 'KAC'),
                        aes(x = `treatment`,y=Fraction))+
    geom_boxplot( width = .3,
                 alpha = .3,outlier.shape= NA)  +
        geom_point() +
    xlab('') +
    theme_classic() + facet_wrap(~timepoint) + theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
    print(fig3b_right)

}

##figure3d was generated using JMP Pro v15 with umap 2d dimention reduction data loaded.

##figure3 e violinplots of cnv scores and Kras mutation fration barplot
##Figure3 e
if(TRUE) {
        cnvdata = readr::read_tsv('Input_data/Fig3e_cnv_dataInput.txt')
        colors = c('AT1' = '#b56cab','AT2' = '#0d7072','KAC' = '#a59c35','Malignant' ='#de1e26' )
    fig3etop = ggplot(data = cnvdata,aes(x = `celltype`,y=cnvscore))+
        geom_jitter(width = .2) +
        geom_violin(scale = 'width',width=.5,aes(fill = celltype),
                    alpha = .8,outlier.shape= NA)  +
        scale_fill_manual(values = colors) +
    xlab('') +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
        print(fig3etop)

        ##bottome part of the panel showing the fraction of G12D mutaiton cell fraction(Data from supplementary Table S11)
        data = readr::read_tsv('Input_Data/Fig3e_KrasFraction_dataInput.txt')
        fig3ebottom = ggplot(data = data,aes(x = `Group`,y=Fraction))+
    geom_bar(width = .3,stat = 'identity',aes(fill = Group)) +
    xlab('') +
        scale_fill_manual(values =colors) +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
        print(fig3ebottom)
}

