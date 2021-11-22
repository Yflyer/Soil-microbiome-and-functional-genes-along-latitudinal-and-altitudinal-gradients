library(geosphere)
source('0_function.R')
set.seed(999)

#################### save dir #######################
folder = paste0('3_tnst')
dir.create(folder)
# This code is to investigte the gene decay relation ship
#################### Latitude (NA) ############################
load(file = '3_tnst/Latitude (North America)_GeoChip/tnst_Latitude (North America)_GeoChip')
load('project_Latitude (North America)_GeoChip')

tnst = get(paste0('tnst_',prefix,'_',note))
group_tnst = tnst$index.pair.grp
group_by(group_tnst,group) %>% summarise(mean = mean(MST.ij.ruzicka))
mean(group_tnst$MST.ij.ruzicka)

project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

# sort the site along the gradient
Order = ENV[['Lat']] %>% order(.) # Lat grad run this
G.Grad = ENV[['Site']][Order] %>% unique(.) %>% as.character(.)
group_tnst$Study=prefix
levels(group_tnst$group) <- c(levels(group_tnst$group), "CWT")
group_tnst$group[group_tnst$group=="CBF"]="CWT"

assign(paste0('group_tnst',prefix),group_tnst)
assign(paste0('G.Grad_',prefix),G.Grad)
#################### Altitude (HKV) ########################
load(file = '3_tnst/Altitude (HKV)_GeoChip/tnst_Altitude (HKV)_GeoChip')
load('project_Altitude (HKV)_GeoChip')

tnst = get(paste0('tnst_',prefix,'_',note))
group_tnst = tnst$index.pair.grp
group_by(group_tnst,group) %>% summarise(mean = mean(MST.ij.ruzicka))
mean(group_tnst$MST.ij.ruzicka)

project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

# sort the site along the gradient
Order = ENV[['Elevation']] %>% order(.) # Ele grad run this
G.Grad = ENV[['Site']][Order] %>% unique(.) %>% as.character(.)
group_tnst$Study=prefix


assign(paste0('group_tnst',prefix),group_tnst)
assign(paste0('G.Grad_',prefix),G.Grad)
#################### Altitude (SNNR) ########################
load(file = '3_tnst/Altitude (SNNR)_GeoChip/tnst_Altitude (SNNR)_GeoChip')
load('project_Altitude (SNNR)_GeoChip')

tnst = get(paste0('tnst_',prefix,'_',note))
group_tnst = tnst$index.pair.grp
group_by(group_tnst,group) %>% summarise(mean = mean(MST.ij.ruzicka))
mean(group_tnst$MST.ij.ruzicka)

project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

# sort the site along the gradient
Order = ENV[['Elevation']] %>% order(.) # Ele grad run this
G.Grad = ENV[['Site']][Order] %>% unique(.) %>% as.character(.)
group_tnst$Study=prefix


assign(paste0('group_tnst',prefix),group_tnst)
assign(paste0('G.Grad_',prefix),G.Grad)

group_tnst = Reduce(function(x, y) rbind(x,y),list(`group_tnstAltitude (HKV)`,`group_tnstLatitude (North America)`,`group_tnstAltitude (SNNR)`))
group_tnst$Study = factor(group_tnst$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
group_tnst$group = factor(group_tnst$group,levels = c(`G.Grad_Latitude (North America)`,`G.Grad_Altitude (SNNR)`,`G.Grad_Altitude (HKV)`))
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
st_mst <- ggplot(aes(x=group,y=MST.ij.ruzicka),data = group_tnst)+
  geom_boxplot(aes(color=Study),lwd=0.7)+
  geom_jitter(aes(color=Study),shape=16,alpha=0.3,position=position_jitter(0.2),size=1)
plot = st_mst+PTs+ylim(0,1)+facet_grid(~Study,scales = 'free')
theme_set(theme_bw()+theme(legend.position = "none"))
plot+scale_color_economist()
#plot+scale_color_viridis(discrete=T,alpha=0.9)
# +stat_boxplot(geom='errorbar',width=0.15)
ggsave(filename = paste0(folder,'/','Group_',note,'_MST.pdf'),dpi=900,width=16,height=7,units='cm')
ggsave(filename = paste0(folder,'/','Group_',note,'_MST.jpg'),dpi=900,width=16,height=7,units='cm')


