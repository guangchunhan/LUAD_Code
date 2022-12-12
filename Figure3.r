##figure3 b generated using jmp with the following data
if(TRUE) {
    fractionData = read_tsv('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/mouse_data_1006/Grant figures2021 Sept/scripts/Fig1D data barplots of fraction.jmp')
    fig3b = ggplot(data = fractionData)+
    geom_boxplot(aes(x = `treatment`,y=Fraction), 
                         width = .3,
                 alpha = .3,outlier.shape= NA)  +
        geom_jitter() +
    xlab('') +
    theme_classic() + facet_grid(`Column 18`~timepoint) + ##Column 18 contains the most recent celltype information
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 
}

##figure3 d was generated using jmp v15 with umap 2d dimention reduction data loaded.

##figure3 e violinplots of cnv scores and Kras mutation fration barplot
if(TRUE) {
        cnvdata = read_tsv('/Users/ghan1/Library/CloudStorage/Box-Box/Humam/LUAD(P1-P16)/figures/figures_data/epithelial/EpithelialUmap/Epithelial_normalCells/alveolar_airway_andAllfiguresworkdir/Figure4 mousedata cnvscore.jmp')
    fig3etop = ggplot(data = cnvdata)+
    geom_violin(aes(x = `Column 18`,y=aviv_score_refT), 
                         width = .3,
                 alpha = .3,outlier.shape= NA)  +
        geom_boxplot(width = .2) +
    xlab('') +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='none',
        axis.text.x = element_text(angle=0)
    ) 


        ##bottome part of the panel showing the fraction of G12D mutaiton cell fraction
        data = read_excel('SuppTable.xlsx',sheet = 11)
        fig3ebottom = ggplot(data = data)+
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

