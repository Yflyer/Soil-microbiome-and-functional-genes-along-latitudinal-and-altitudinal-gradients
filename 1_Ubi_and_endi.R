source('0_function.R')
library(scales)

#################### save dir #######################
folder = paste0('1_abund-occur')
dir.create(folder)
#########################################################################
# This code is to investigte the gene distribution relation ship
load('project_Latitude (North America)_GeoChip')
load('project_Altitude (SNNR)_GeoChip')
load('project_Altitude (HKV)_GeoChip')
######################################
load('project_Latitude (North America)_OTU')
load('project_Altitude (SNNR)_OTU')
load('project_Altitude (HKV)_OTU')
######################################

############################## repeat following codes on three projects ######
project = get(paste0('project_',prefix,'_',note)) #%>% rarefy_even_depth(.,sample.size = 10000)

OTU <- otu_table(project)
ENV <- sample_data(project)
# define levels
Site = ENV$Site %>% as.factor(.)
# = ENV$forest_name %>% as.factor(.)

# NA existed
nrow(OTU)
OTU = OTU %>% .[rowSums(.)>0,] %>% .[sample(c(1:nrow(OTU)),size = 14000),] 
nrow(OTU)
#########################################################################
########################### region search ###############################
site_count_table=group_count(table = OTU, group = Site)
# get the mean relative abundance
Occupancy = rowSums(site_count_table)
Abundance=apply(OTU, 1, function(x) sum(x[x>0])/sum(x>0))
Site_occupancy=c()
occur_frequency=1
for (i in 1:nlevels(Site)) {
  Site_occupancy[which(apply(site_count_table, 1, function(x) sum(x>=occur_frequency)==i))]=paste(i)
}
ab_oc_data=data.frame(Occupancy,Abundance,Site_occupancy)
ab_oc_data$Site_occupancy=factor(ab_oc_data$Site_occupancy,levels = c(1:nlevels(Site))) %>% droplevels(.)

# plot data process
ab_oc_data$Log.abundance =log(Abundance)
ab_oc_data$Log.occupancy =log(Occupancy)
ab_oc_data$Study=prefix
ab_oc_data$Data=note
ab_oc_data$Type="Occasional"
ab_oc_data$Type[ab_oc_data$Site_occupancy==nlevels(Site)]="Ubiquitous"
ab_oc_data$Type[ab_oc_data$Site_occupancy==1]="Endemic"
ab_oc_data$Type=as.factor(ab_oc_data$Type)
# check
ab_oc_data[is.na(ab_oc_data)]
assign(paste0('ab_oc_data_',prefix),ab_oc_data)
##################### repeat above codes on three projects #####################

########################## merge three ab—oc data
ab_oc_data = Reduce(rbind,list(`ab_oc_data_Latitude (North America)`,`ab_oc_data_Altitude (SNNR)`,`ab_oc_data_Altitude (HKV)`))
ab_oc_data$Study = factor(ab_oc_data$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
########################## counts of the type
Type_count = group_by(ab_oc_data,Type,Study) %>% summarise(Counts=n()) %>% as.data.frame()
Total = summary(ab_oc_data$Study)
Type_count = data.frame(Type_count,Total)
Type_count$Type_ratio = scales::label_percent(accuracy = 0.001)(Type_count$Counts/Type_count$Total)
write.csv(file = paste0(folder,'/Type_count_',note,'.csv'),Type_count)

#############################################
# AO-plot
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
PTs <- theme(legend.title = element_blank(),
             legend.position="bottom",
             #legend.key.size = unit(2, 'cm'),
             #legend.key.height= unit(2, 'cm'),
             #legend.key.width= unit(2, 'cm'),
             legend.key=element_blank(),
             #axis.text = element_text(size=16, hjust = 0.5,vjust = 0.5),
             panel.background=element_blank(),
             panel.border = element_rect(colour = "grey", fill=NA, size=1))
syc = scale_y_continuous(labels = scales::number_format(accuracy = 0.01,decimal.mark = '.'))
guides = guides(color = guide_legend(override.aes = list(size = 5)))

plot_name='1_AO_log_'
ab_oc <- ggplot(ab_oc_data,aes(x=Log.occupancy,y=Log.abundance))+
  geom_point(alpha = 0.3,size=0.8,shape=16,aes(color=Type))
ab_oc+PTs+facet_grid(~Study,scales = 'free_x')+syc+guides+scale_color_wsj()
ggsave(filename = paste0(folder,'/',plot_name,'_',note,'.jpg'),dpi=900,width=17,height=8.5,units='cm')
ggsave(filename = paste0(folder,'/',plot_name,'_',note,'.pdf'),dpi=900,width=17,height=8.5,units='cm')

plot_name='1_AO_'
ab_oc <- ggplot(ab_oc_data,aes(x=Occupancy,y=Log.abundance))+
  geom_point(alpha = 0.3,size=0.8,shape=16,aes(color=Type))
ab_oc+PTs+facet_grid(~Study,scales = 'free_x')+syc+guides+scale_color_wsj()
ggsave(filename = paste0(folder,'/',plot_name,'_',note,'.jpg'),dpi=900,width=17,height=8.5,units='cm')
ggsave(filename = paste0(folder,'/',plot_name,'_',note,'.pdf'),dpi=900,width=17,height=8.5,units='cm')



