source('0_function.R')

#################### save dir #######################
folder = paste0('1_abund-occur')
dir.create(folder)
#########################################################################
# This code is to investigte the gene distribution relation ship

ao_pack <- function(project){
  ############################## repeat following codes on three projects ######
  OTU <- otu_table(project)
  ENV <- sample_data(project)
  # define levels
  Site = ENV$Site %>% as.factor(.)
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
  ab_oc_data$Log.abundance =log1p(Abundance)
  ab_oc_data$Log.occupancy =log1p(Occupancy)
  ab_oc_data$Study=prefix
  ab_oc_data$Data=note
  ab_oc_data$Type="Occasional"
  ab_oc_data$Type[ab_oc_data$Site_occupancy==nlevels(Site)]="Ubiquitous"
  ab_oc_data$Type[ab_oc_data$Site_occupancy==1]="Endemic"
  ab_oc_data$Type=as.factor(ab_oc_data$Type)
  # check
  ab_oc_data[is.na(ab_oc_data)]
  ab_oc_data
  ##################### repeat above codes on three projects #####################
}
############### GeoChip #############
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
ab_oc_data = ao_pack(project) 
#ab_oc_data = ab_oc_data[sample(1:nrow(ab_oc_data),10000),]
assign(paste0('ab_oc_data_',prefix),ab_oc_data)
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
ab_oc_data = ao_pack(project) 
#ab_oc_data = ab_oc_data[sample(1:nrow(ab_oc_data),10000),]
assign(paste0('ab_oc_data_',prefix),ab_oc_data)
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
ab_oc_data = ao_pack(project) 
ab_oc_data = ab_oc_data[complete.cases(ab_oc_data),]
#ab_oc_data = ab_oc_data[sample(1:nrow(ab_oc_data),10000),]
assign(paste0('ab_oc_data_',prefix),ab_oc_data)


########################## merge three a-o data
ab_oc_data = Reduce(rbind,list(`ab_oc_data_Latitude (North America)`,`ab_oc_data_Altitude (SNNR)`,`ab_oc_data_Altitude (HKV)`))
ab_oc_data$Study = factor(ab_oc_data$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
ab_oc_data$Alpha = 0.6
ab_oc_data$Alpha[ab_oc_data$Type1!='Ubiquitous'] = 0.2
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
             panel.grid=element_blank(),
             panel.background=element_blank(),
             panel.border = element_rect(colour = "grey", fill=NA, size=1))
syc = scale_y_continuous(labels = scales::number_format(accuracy = 0.01,decimal.mark = '.'))
guides = guides(color = guide_legend(override.aes = list(size = 5)))
##### jitter
pp <- ggplot(aes(x=Log.occupancy,y=Log.abundance),data = ab_oc_data)+
  geom_jitter(aes(color=Type,alpha = Alpha),shape=16,position=position_jitter(0.1),size=0.5)
plot = pp+PTs+facet_grid(~Study,scales = 'free_x')
{plot+syc+guides+
    scale_color_manual(values = c('red','gainsboro','darkturquoise'))+
    scale_alpha_identity()
}
ggsave(filename = paste0(folder,'/','Abund_occur_jitter_',note,'.pdf'),dpi=600,width=16,height=8,units='cm')
ggsave(filename = paste0(folder,'/','Abund_occur_jitter_',note,'.jpg'),dpi=600,width=16,height=8,units='cm')
###
##### point
pp <- ggplot(aes(x=Log.occupancy,y=Log.abundance),data = ab_oc_data)+
  geom_point(aes(color=Type,alpha = Alpha),shape=16,size=0.5)
plot = pp+PTs+facet_grid(~Study,scales = 'free_x')
{plot+syc+guides+
    scale_color_manual(values = c('red','gainsboro','darkturquoise'))+
    scale_alpha_identity()
}
ggsave(filename = paste0(folder,'/','Abund_occur_point_',note,'.pdf'),dpi=600,width=16,height=8,units='cm')
ggsave(filename = paste0(folder,'/','Abund_occur_point_',note,'.jpg'),dpi=600,width=16,height=8,units='cm')


############### OTU ##############
load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))
ab_oc_data = ao_pack(project) 
#ab_oc_data = ab_oc_data[sample(1:nrow(ab_oc_data),10000),]
assign(paste0('ab_oc_data_',prefix),ab_oc_data)
load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))
ab_oc_data = ao_pack(project) 
#ab_oc_data = ab_oc_data[sample(1:nrow(ab_oc_data),10000),]
assign(paste0('ab_oc_data_',prefix),ab_oc_data)
load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))
ab_oc_data = ao_pack(project) 
#ab_oc_data = ab_oc_data[sample(1:nrow(ab_oc_data),10000),]
assign(paste0('ab_oc_data_',prefix),ab_oc_data)

########################## merge three a-o data
ab_oc_data = Reduce(rbind,list(`ab_oc_data_Latitude (North America)`,`ab_oc_data_Altitude (SNNR)`,`ab_oc_data_Altitude (HKV)`))
ab_oc_data$Study = factor(ab_oc_data$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
ab_oc_data$Alpha = 0.6
ab_oc_data$Alpha[ab_oc_data$Type!='Endemic'] = 0.2
########################## counts of the type
Type_count = group_by(ab_oc_data,Type,Study) %>% summarise(Counts=n()) %>% as.data.frame()
Total = summary(ab_oc_data$Study)
Type_count = data.frame(Type_count,Total)
Type_count$Type_ratio = scales::label_percent(accuracy = 0.001)(Type_count$Counts/Type_count$Total)
write.csv(file = paste0(folder,'/Type_count_',note,'.csv'),Type_count)

##### jitter
pp <- ggplot(aes(x=Log.occupancy,y=Log.abundance),data = ab_oc_data)+
  geom_jitter(aes(color=Type,alpha = Alpha),shape=16,position=position_jitter(0.1),size=0.5)
plot = pp+PTs+facet_grid(~Study,scales = 'free_x')
{plot+syc+guides+
    scale_color_manual(values = c('red','gainsboro','darkturquoise'))+
    scale_alpha_identity()
}
ggsave(filename = paste0(folder,'/','Abund_occur_jitter_',note,'.pdf'),dpi=600,width=16,height=8,units='cm')
ggsave(filename = paste0(folder,'/','Abund_occur_jitter_',note,'.jpg'),dpi=600,width=16,height=8,units='cm')
###
##### point
pp <- ggplot(aes(x=Log.occupancy,y=Log.abundance),data = ab_oc_data)+
  geom_point(aes(color=Type,alpha = Alpha),shape=16,size=0.5)
plot = pp+PTs+facet_grid(~Study,scales = 'free_x')
{plot+syc+guides+
    scale_color_manual(values = c('red','gainsboro','darkturquoise'))+
    scale_alpha_identity()
}
ggsave(filename = paste0(folder,'/','Abund_occur_point_',note,'.pdf'),dpi=600,width=16,height=8,units='cm')
ggsave(filename = paste0(folder,'/','Abund_occur_point_',note,'.jpg'),dpi=600,width=16,height=8,units='cm')
