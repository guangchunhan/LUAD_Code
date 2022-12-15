##figure3 b
if(TRUE) {
    library(ggplot2)
    fractionData = read_tsv('Input_Data/Fig3_FractionChange_dataInput.txt')
    fractionData$treatment = factor(fractionData$treatment,
                                    levels = c('Saline','NNK'))

fractionData$timepoint = factor(fractionData$timepoint,
                                    levels = c('EOE','7 Mo'))

    fig3b_left = ggplot(data = filter(fractionData,celltype == 'Malignant'),
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


fig3b_right = ggplot(data = filter(fractionData,celltype == 'KAC'),
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
        cnvdata = read_tsv('Input_data/Fig3e_cnv_dataInput.txt')
    fig3etop = ggplot(data = cnvdata,aes(x = `celltype`,y=cnvscore))+
        geom_jitter(width = .2) +
        geom_violin(scale = 'width',width=.5,
                    alpha = .8,outlier.shape= NA)  +
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
        data = read_tsv('Input_Data/Fig3e_KrasFraction_dataInput.txt')
        fig3ebottom = ggplot(data = data,aes(x = `Group`,y=Fraction))+
    geom_bar(width = .3,stat = 'identity') +
    xlab('') +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
        print(fig3ebottom)
}

