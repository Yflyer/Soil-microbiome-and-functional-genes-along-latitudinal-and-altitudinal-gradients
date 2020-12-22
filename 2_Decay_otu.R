library(geosphere)
source('0_function.R')
set.seed(999)

load('project_HA_OTU')
load('project_NA_OTU')
load('project_SNJ_OTU')

#################### save dir #######################
folder = paste0('2_Decay_',note)
dir.create(folder)

##############################################################
#################### Altitude (SNNR)
prefix='SNNR'
project = project_SNJ_OTU
OTU <- otu_table(project)
ENV <- sample_data(project) #%>% as.data.frame(.) 

#################### 
jaccard = vegdist(t(OTU), method="jaccard", binary=FALSE) %>% col3m(.,dist_name ='jaccard')
bray = vegdist(t(OTU), method="bray", binary=FALSE) %>% col3m(.,dist_name ='bray')
sorensen = vegdist(t(OTU), method="bray", binary=TRUE) %>% col3m(.,dist_name ='sorensen')
unifrac = UniFrac(project, weighted=T, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'unifrac')
uw.unifrac = UniFrac(project, weighted=F, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'uw.unifrac')
distance = geo_dist(ENV) %>% col3m(.,dist_name = 'distance')
#################### elevation distance
distance.ele = ENV[,'Elevation'] %>% dist(.)/1000
distance.ele = col3m(distance.ele,dist_name = 'distance') 
distance$distance=distance$distance+distance.ele$distance
#################### 
decay_data = Reduce(function(x, y) inner_join(x, y,by=c('row','col')),list(distance,jaccard,bray,sorensen,unifrac,uw.unifrac))

#################### decay model examniation
decay_data$Similarity.bray=1-decay_data$bray
decay_data$ln.distance = log1p(decay_data$distance)
model.result = betapart::decay.model(decay_data$Similarity.bray,decay_data$ln.distance,perm=999)
model.result$p.value
# snj: slope -0.1938805 p 0.5793936 p.p 0.6823177 p.b 0.3316683
#################################################################
decay_data$Study = 'Altitude (SNNR)'
assign(paste0('decay_data_',prefix),decay_data)

#################################################################
#################### Altitude (HKV)
prefix='HKV'
project = project_HA_OTU
OTU <- otu_table(project)
ENV <- sample_data(project) #%>% as.data.frame(.) 
#################### 
jaccard = vegdist(t(OTU), method="jaccard", binary=FALSE) %>% col3m(.,dist_name ='jaccard')
bray = vegdist(t(OTU), method="bray", binary=FALSE) %>% col3m(.,dist_name ='bray')
sorensen = vegdist(t(OTU), method="bray", binary=TRUE) %>% col3m(.,dist_name ='sorensen')
unifrac = UniFrac(project, weighted=T, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'unifrac')
uw.unifrac = UniFrac(project, weighted=F, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'uw.unifrac')
distance = geo_dist(ENV) %>% col3m(.,dist_name = 'distance')
#################### elevation distance
distance.ele = ENV[,'Elevation'] %>% dist(.)/1000
distance.ele = col3m(distance.ele,dist_name = 'distance') 
distance$distance=distance$distance+distance.ele$distance
#################### 
decay_data = Reduce(function(x, y) inner_join(x, y,by=c('row','col')),list(distance,jaccard,bray,sorensen,unifrac,uw.unifrac))
#################### decay model examniation
decay_data$Similarity.bray=1-decay_data$bray
decay_data$ln.distance = log1p(decay_data$distance)
model.result = betapart::decay.model(decay_data$Similarity.bray,decay_data$ln.distance,perm=999)
model.result$p.value
# ha: slope -0.1283425 p =0.008,P.PERM = 0.313 p.b = 0.174
####################################################
decay_data$Study = 'Altitude (HKV)'
assign(paste0('decay_data_',prefix),decay_data)

#################################################################
#################### Latitude (NA)
prefix='NA'
project = project_NA_OTU
OTU <- otu_table(project)
ENV <- sample_data(project) #%>% as.data.frame(.) 
#################### 
jaccard = vegdist(t(OTU), method="jaccard", binary=FALSE) %>% col3m(.,dist_name ='jaccard')
bray = vegdist(t(OTU), method="bray", binary=FALSE) %>% col3m(.,dist_name ='bray')
sorensen = vegdist(t(OTU), method="bray", binary=TRUE) %>% col3m(.,dist_name ='sorensen')
unifrac = UniFrac(project, weighted=T, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'unifrac')
uw.unifrac = UniFrac(project, weighted=F, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'uw.unifrac')
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
#################### 
decay_data = Reduce(function(x, y) inner_join(x, y,by=c('row','col')),list(distance,jaccard,bray,sorensen,unifrac,uw.unifrac))
#################### decay model examniation
decay_data$Similarity.bray=1-decay_data$bray
decay_data$ln.distance = log1p(decay_data$distance)
model.result = betapart::decay.model(decay_data$Similarity.bray,decay_data$ln.distance,perm=999)
model.result$p.value
# na: slope -0.0847966 p <0.001 p.p 0.012 p.b 0.003 p 0.0002408902
####################################################
decay_data$Study = 'Latitude (North America)'
assign(paste0('decay_data_',prefix),decay_data)



# plot
PTs <- theme(#text = element_text(size=12),
             #axis.text.x = element_text(size=16, hjust = 0.5,vjust = 0.5),
             panel.background=element_blank(),
             legend.position = "none",
             panel.border = element_rect(colour = "grey", fill=NA, size=1))

decay_data = Reduce(function(x, y) rbind(x,y),list(decay_data_HKV,decay_data_NA,decay_data_SNNR))

decay_data$Study = factor(decay_data$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))

# log.distance
SDecay <- ggplot(decay_data,aes(x=ln.distance,y=Similarity.bray))+
  geom_point(alpha = 0.5,stroke = 0.3,size=1.5,shape=21,aes(fill=Study))+geom_smooth(method='lm',aes(color=Study),se = F)
SDecay+PTs+ylim(c(0,1))+facet_grid(~Study,scales = 'free')#+lm_Site+Ps_Sitescm_Site+sfm_Site
ggsave(filename = paste0(folder,'/','3_Decay_',note,'_log_bray.pdf'),dpi=900,width=19.5,height=9,units='cm')
ggsave(filename = paste0(folder,'/','3_Decay_',note,'_log_bray.jpg'),dpi=900,width=19.5,height=9,units='cm')

# distance
SDecay <- ggplot(decay_data,aes(x=distance,y=Similiarity_jaccard))+
  geom_point(alpha = 0.5,stroke = 0.1,size=1,shape=21,aes(fill=Study))+geom_smooth(method='lm',aes(color=Study))
SDecay+PTs+ylim(c(0,1))+facet_grid(~Study,scales = 'free')#+lm_Site+Ps_Sitescm_Site+sfm_Site
ggsave(filename = paste0('3_Decay_',prefix,'_unlog_jaccard.pdf'),dpi=900,width=19.5,height=9,units='cm')


