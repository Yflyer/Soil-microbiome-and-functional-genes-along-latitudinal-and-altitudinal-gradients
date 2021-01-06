library(geosphere)
source('0_function.R')
set.seed(999)

#################### save dir #######################
folder = paste0('3_tnst')
dir.create(folder)
# This code is to investigte the gene decay relation ship
#################### Latitude (NA) ############################
load('3_tnst/tnst_mr_Latitude (North America)_GeoChip')
load('project_Latitude (North America)_GeoChip')

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
group_tnst_mr = tnst_mr$index.pair.grp
group_by(group_tnst_mr,group) %>% summarise(mean = mean(MST.ij.bray))
mean(group_tnst_mr$MST.ij.bray)

project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

# sort the site along the gradient
Order = ENV[['Lat']] %>% order(.) # Lat grad run this
G.Grad = ENV[['Site']][Order] %>% unique(.)
group_tnst_mr$Study=prefix
group_tnst_mr$group=factor(group_tnst_mr$group,levels = G.Grad)

assign(paste0('group_tnst_mr',prefix),group_tnst_mr)

#################### Altitude (HKV) ########################
load('3_tnst/tnst_mr_Altitude (HKV)_GeoChip')
load('project_Altitude (HKV)_GeoChip')

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
group_tnst_mr = tnst_mr$index.pair.grp
group_by(group_tnst_mr,group) %>% summarise(mean = mean(MST.ij.bray))
mean(group_tnst_mr$MST.ij.bray)
#################### Altitude (SNNR) ########################
load('3_tnst/tnst_mr_Altitude (SNNR)_GeoChip')
load('project_Altitude (SNNR)_GeoChip')

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
group_tnst_mr = tnst_mr$index.pair.grp
group_by(group_tnst_mr,group) %>% summarise(mean = mean(MST.ij.bray))
mean(group_tnst_mr$MST.ij.bray)

project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

# sort the site along the gradient
Order = ENV[['Elevation']] %>% order(.) # Ele grad run this
G.Grad = ENV[['Site']][Order] %>% unique(.)
group_tnst_mr$Study=prefix
group_tnst_mr$group=factor(group_tnst_mr$group,levels = G.Grad)

assign(paste0('group_tnst_mr',prefix),group_tnst_mr)

group_tnst_mr = Reduce(function(x, y) rbind(x,y),list(`group_tnst_mrAltitude (HKV)`,`group_tnst_mrLatitude (North America)`,`group_tnst_mrAltitude (SNNR)`))
group_tnst_mr$Study = factor(group_tnst_mr$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
##### MST Group plot
########################################################
library(ggthemes)
library(viridis)
PTs <- theme(#text = element_text(size=12),
  #axis.text.x = element_text(size=16, hjust = 0.5,vjust = 0.5),
  panel.background=element_blank(),
  legend.position = "none",
  panel.border = element_rect(colour = "grey", fill=NA, size=1))
### mst plot
st_mr_mst <- ggplot(aes(x=group,y=MST.ij.bray),data = group_tnst_mr)+
  geom_boxplot(aes(color=Study),lwd=0.7)+
  geom_jitter(aes(color=Study),shape=16,alpha=0.3,position=position_jitter(0.2),size=1)
plot = st_mr_mst+PTs+ylim(0,1)+facet_grid(~Study,scales = 'free')
theme_set(theme_bw()+theme(legend.position = "none"))
plot+scale_color_economist()
#plot+scale_color_viridis(discrete=T,alpha=0.9)
# +stat_boxplot(geom='errorbar',width=0.15)
ggsave(filename = paste0(folder,'/','Group_',note,'_MST.pdf'),dpi=900,width=19.5,height=7.5,units='cm')
ggsave(filename = paste0(folder,'/','Group_',note,'_MST.jpg'),dpi=900,width=19.5,height=7.5,units='cm')


