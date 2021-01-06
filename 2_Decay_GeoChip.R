library(geosphere)
source('0_function.R')
set.seed(999)

#################### save dir #######################
folder = paste0('2_Decay')
dir.create(folder)
# This code is to investigte the gene decay relation ship
################################################

#################### Altitude (SNNR) ########################
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

jaccard = vegdist(t(OTU), method="jaccard", binary=FALSE) %>% col3m(.,dist_name ='jaccard')
bray = vegdist(t(OTU), method="bray", binary=FALSE) %>% col3m(.,dist_name ='bray')
sorensen = vegdist(t(OTU), method="bray", binary=TRUE) %>% col3m(.,dist_name ='sorensen')
distance = geo_dist(ENV) %>% col3m(.,dist_name = 'distance')
#################### elevation distance
distance.ele = ENV[,'Elevation'] %>% dist(.)/1000
distance.ele = col3m(distance.ele,dist_name = 'distance') 
distance$distance=distance$distance+distance.ele$distance

decay_data = Reduce(function(x, y) inner_join(x, y,by=c('row','col')),list(distance,jaccard,bray,sorensen))

#################### decay model examniation
decay_data$Similarity.bray=1-decay_data$bray
decay_data$ln.distance = log1p(decay_data$distance)
#################################################################
decay_data$Study = prefix
assign(paste0('decay_data_',prefix),decay_data)
################################ SNNR #################################

#################### Altitude (HKV) ############################### ##################################
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) #%>% as.data.frame(.) 

jaccard = vegdist(t(OTU), method="jaccard", binary=FALSE) %>% col3m(.,dist_name ='jaccard')
bray = vegdist(t(OTU), method="bray", binary=FALSE) %>% col3m(.,dist_name ='bray')
sorensen = vegdist(t(OTU), method="bray", binary=TRUE) %>% col3m(.,dist_name ='sorensen')
distance = geo_dist(ENV) %>% col3m(.,dist_name = 'distance')
#################### elevation distance
distance.ele = ENV[,'Elevation'] %>% dist(.)/1000
distance.ele = col3m(distance.ele,dist_name = 'distance') 
distance$distance=distance$distance+distance.ele$distance

decay_data = Reduce(function(x, y) inner_join(x, y,by=c('row','col')),list(distance,jaccard,bray,sorensen))
#################### decay model examniation
decay_data$Similarity.bray=1-decay_data$bray
decay_data$ln.distance = log1p(decay_data$distance)
decay_data$Study = prefix
assign(paste0('decay_data_',prefix),decay_data)

#################### Latitude (NA) ##############################
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) #%>% as.data.frame(.) 

jaccard = vegdist(t(OTU), method="jaccard", binary=FALSE) %>% col3m(.,dist_name ='jaccard')
bray = vegdist(t(OTU), method="bray", binary=FALSE) %>% col3m(.,dist_name ='bray')
sorensen = vegdist(t(OTU), method="bray", binary=TRUE) %>% col3m(.,dist_name ='sorensen')
distance = geo_dist(ENV) %>% col3m(.,dist_name = 'distance')
#################### site-inner distance (for latitudinal dataset)
inner_distance = data.frame()
for (i in 1:nlevels(ENV$Site)) {
  dist = ENV[ENV$Site==levels(ENV$Site)[i],c('X','Y')] %>% dist(.)/1000 
  dist = col3m(dist,dist_name = 'inner_distance')
  inner_distance=rbind(inner_distance,dist)
}
inner_distance = left_join(distance,inner_distance,by=c('row','col'))
inner_distance[is.na(inner_distance)]=0
distance$distance=inner_distance$distance+inner_distance$inner_distance
####################################################
decay_data = Reduce(function(x, y) inner_join(x, y,by=c('row','col')),list(distance,jaccard,bray,sorensen))
decay_data$Similarity.bray=1-decay_data$bray
decay_data$ln.distance = log1p(decay_data$distance)
decay_data$Study = 'Latitude (North America)'
assign(paste0('decay_data_',prefix),decay_data)

model.result = betapart::decay.model(decay_data$Similarity.bray,decay_data$ln.distance,perm=999)
model.result$p.value
############################### NA #########################

decay_data = Reduce(function(x, y) rbind(x,y),list(`decay_data_Altitude (HKV)`,`decay_data_Latitude (North America)`,`decay_data_Altitude (SNNR)`))
decay_data$Study = factor(decay_data$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
# plot
library(ggthemes)
PTs <- theme(#text = element_text(size=12),
             #axis.text.x = element_text(size=16, hjust = 0.5,vjust = 0.5),
             panel.background=element_blank(),
             legend.position = "none",
             panel.border = element_rect(colour = "grey", fill=NA, size=1))

# log.distance
syc = scale_y_continuous(labels = scales::number_format(accuracy = 0.01,decimal.mark = '.'),limits=c(0.4,1))
sxc = scale_x_continuous(labels = scales::number_format(accuracy = 0.1,decimal.mark = '.'))

SDecay <- ggplot(decay_data,aes(x=ln.distance,y=Similarity.bray))+geom_point(alpha = 0.5,stroke = 0,size=1.5,shape=16,aes(color=Study))+geom_smooth(method='lm',aes(color=Study),se = F)
plot = SDecay+facet_grid(~Study,scales = 'free')#+lm_Site+Ps_Sitescm_Site+sfm_Site
theme_set(theme_bw()+theme(legend.position = "none"))
plot+scale_color_economist()+syc+sxc
ggsave(filename = paste0(folder,'/','3_Decay_',note,'_log_bray.pdf'),dpi=900,width=19.5,height=7.5,units='cm')
ggsave(filename = paste0(folder,'/','3_Decay_',note,'_log_bray.jpg'),dpi=900,width=19.5,height=7.5,units='cm')



