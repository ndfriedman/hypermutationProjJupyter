#written by Noah Friedman
#the code to plot all the figures for figure 1
library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(stringr)




emptyTheme <- theme(axis.line = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())
#
###
######
###########
#######
###
#

#PLOT FIGURE 1A

#please refer to 
#hypermutationAnalysisProject/plottingScripts/plotAndDefineHypermutationThresholds.R
#and
#/Users/friedman/Desktop/WORK/hypermutationProjectJupyterScripts/scriptsToGenerateFigures/generate_mut_classification_figure#1a#.ipynb

#
###
######
###########
#######
###
#

#PLOT FIGURE 1B

plot_n_cases_figure <- function(df){
  
  p <- ggplot(df, aes(x=reorder(cancerType, -nTotal)))+
    geom_bar(aes(y=nTotal, fill='Total Cases'), stat='identity')+
    #optional include the n high mut burden as well
    geom_bar(aes(y=-1*nHypermutated - nHighMutBurden, fill='High mutation burden'), stat='identity')+
    
    geom_bar(aes(y=-1*nHypermutated, fill='Hypermutated'), stat='identity')+
    
    theme(axis.text.x = element_text(angle=90))+
    theme(axis.ticks.x = element_blank())+
    scale_fill_manual(values=c('#858585', 'black', 'light gray'))+
    scale_y_continuous(breaks= c(-500,0,2000, 4000, 6000), labels=c('500', '0', '2000', '4000', '6000'))+
    emptyTheme+
    guides(fill=guide_legend(title="Tumor Classification"))+
    xlab('Cancer Type')+
    ylab('Total Cases')
  return(p)
}


dfCancerTypes <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/figure1bCancerTypeSummary.tsv', sep='\t', header=TRUE)
p <- plot_n_cases_figure(dfCancerTypes)
p <- plot_grid(p, ggplot()+labs(caption='plotFigure1b_1c.R\ngenerate_cancer_types_and_sigs_figure_#figure1b#'), nrow=2, rel_heights = c(1,.1))
ggsave('~/Desktop/plot.pdf', plot=p,  width = 6, height = 6)


#
###
######
###########
#######
###
#

#PLOT FIGURE 1C

plot_signatures_figure <- function(df){
  p <- ggplot(df, aes(x=reorder(cancerType, -orderingVal), y=nCases))+
    
    geom_bar(aes(fill=
                   factor(signature, levels=c('APOBEC', 'MMR',
                                              'other', 'POLE', 'POLE&MMR', 'SMOKING', 'TMZ', 'UV')))
             , stat='identity')+
    theme(axis.text.x = element_text(angle=90))+
    theme(axis.ticks.x = element_blank())+
    xlab('Cancer Type')+
    scale_fill_manual(values=c("#FF0000","#267574",
                               'gray', "#ADFF2F","#00dd5d", '#ffb347', '#2A52BE', "#FFF600"))+
    guides(fill=guide_legend(title="Signature"))+
    emptyTheme
  return(p)
}

dfSigs <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/figure1cSignatureSummary.tsv', sep='\t', header=TRUE)
p <- plot_signatures_figure(dfSigs)
p <- plot_grid(p, ggplot()+labs(caption='plotFigure1b_1c.R\ngenerate_cancer_types_and_sigs_figure_#figure1b#'), nrow=2, rel_heights = c(1,.1))
ggsave('~/Desktop/plot.pdf', plot=p,  width = 6, height = 6)


#
###
######
###########
#######
###
#

#PLOT FIGURE 1D

plot_data <- function(df){
  plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nHotspots))+
  #plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nOncMuts))+
    geom_boxplot(fatten = NULL, outlier.shape=NA)+
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    #scale_colour_manual(values =  c('black', "#267574", 'gray', '#ADFF2F', "#9acd32", '#2A52BE'), name="Dominant\nSignature")+
    #ylab('N Driver Mutations')+
    ylab('N hotspot mutations')+
    xlab('Cancer Type')+
    emptyTheme+
    geom_jitter(aes(colour=factor(cohort,
                                  levels = c('normal_Colorectal', 'hyper_Colorectal', 'normal_Endometrial', 'hyper_Endometrial',
                                             'normal_Glioma', 'hyper_Glioma', 'normal_Other', 'hyper_Other'))),
                shape=16, position=position_jitter(0.1), alpha=0.75)+
    
    scale_color_manual(values =c('orange', '#b36200', 'lavender', '#301934', '#add8e6', 'blue', 'gray', '#333333'), name='Cohort')
  
  return(plt)
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/nOncMutByCohort.tsv', sep='\t', header=TRUE)
plt <- plot_data(df)
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 6, height = 6, units = c("in"), limitsize = FALSE)



#
####
#########
##############
#########
####
#

#plot figure 1e

emptyTheme <- theme(axis.line = element_blank(),
                    #axis.text = element_blank(),
                    #axis.ticks = element_blank(),
                    #axis.title = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/hotspotObsExpected.tsv',sep = '\t', header=TRUE)

plt <- ggplot(df)+
  #ggplot(df[df$nmutIM341 > 20,])+
   geom_smooth(aes(x=nmutIM341, y=nTotalHotspots, colour='Observed'),se=FALSE, span = 50)+
                         geom_smooth(aes(x=nmutIM341, y=nHotspotsExpected, colour='Expected'), se=FALSE)+
                          xlab('N Nonsynonymous Mutations\nin IMPACT 341 Genes')+
                          ylab('N hotspots per case')+
                          #ggtitle('Observed and expected hotspot burden\nin hypermutated tumors')+
                          emptyTheme+
                          labs(colour = "Number of Hotspots:")+
                          scale_color_manual(values=c('gray', 'black'))


ggsave('~/Desktop/plot.pdf', plot=plt,  width = 4, height = 5, units = c("in"))





#
###
######
#############
##################
##########################
#################
############
#######
###
#

#WORK on Monday Dec 2

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/essentialGenesMuts.tsv', ,sep = '\t', header=TRUE)

ggplot(df, aes(x= nmut))+

  #geom_smooth(aes(y=nEssentialDoubleNormed, colour='essentialDoubleRate'), span=25)+
  #geom_smooth(aes(y=nNeutralDoubleNormed, colour='neutralDoubleRate'), span=25)+
  #geom_smooth(aes(y=nTSGDoubleNormed, colour='tsgDoubleRate'), span=25)
  
  geom_point(aes(y=nEssentialDoubleNormed, colour='essentialDoubleRate'))+
  geom_point(aes(y=nNeutralDoubleNormed, colour='neutralDoubleRate'))+
  geom_point(aes(y=nTSGDoubleNormed, colour='tsgDoubleRate'))











                                                 