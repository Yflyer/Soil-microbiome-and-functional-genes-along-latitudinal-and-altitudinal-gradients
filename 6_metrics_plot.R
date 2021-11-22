source('0_function.R')

#################### save dir #######################
folder = paste0('6_metrics_effect')
dir.create(folder)
####### 
dt = read.csv(file = '6_metrics_effect/stat_data.csv',row.names = 1)
dt$Study = factor(dt$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
dt$Metric = factor(dt$Metric,levels = c('Composition','Abundance','Diversity'))

# get test data ,Log.Dist=mean(Log.Dist)
# test plot
#test = filter(dt,Method=='Bray-Curits'&Study=='Altitude (HKV)') 
test = dt %>% group_by(Method,Comparison,Study,Pair,Metric) %>% summarise(Dist=mean(Dist))
#data=test


library(ggsignif)
library(ggpubr)
library(ggplot2)
### plot box
PTs = theme(axis.text.x=element_blank(),
            axis.title=element_blank(),
            legend.position="bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 8),
            panel.background=element_blank(),
            panel.border = element_rect(colour = "grey", fill=NA, size=1))
my_comparisons <- list(  c("Abundance", "Diversity"), c("Composition", "Abundance"),c("Composition", "Diversity") )
data=dt
## plot
st_mst <- ggplot(aes(x=Metric,y=1-Dist),data = test)+geom_violin(aes(x=Metric,y=1-Dist),data = test,width = 0.5,scale = "width",trim = F,size=1)+geom_jitter(aes(x=Metric,y=1-Dist),data = test,shape=16,alpha=0.7,position=position_jitter(0.2),size=0.5)
# Default method = "kruskal.test" for multiple groups
st_mst+PTs+facet_grid(Method~Study,scales = "free")+stat_compare_means(comparisons = my_comparisons,label = 'p.signif',vjust=0.45) #+ylim(0,60) ,label.y.npc = c(0.2,0.5)
st_mst+PTs+facet_grid(Method~Study,scales = "free")+ expand_limits(y=1.05)
ggsave(filename = paste0(folder,'/Box_plot and site-test.pdf'),dpi=600,width=17,height=20,units='cm')

#  +stat_compare_means(label.y =0.7)

