library(geosphere)
source('0_function.R')
set.seed(999)

#################### save dir #######################
folder = paste0('3_tnst')
dir.create(folder)
# This code is to investigte the gene decay relation ship
#################### Latitude (NA) ############################
load('3_tnst/tnst_mr_Latitude (North America)_OTU')
load('project_Latitude (North America)_OTU')

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
pair_tnst_mr = tnst_mr$index.pair
############ for test #########################
#pair_tnst_mr = reshape_col3m(pair_tnst_mr,value = 'MST.ij.bray',symmetric = T) %>% dist(.) %>% col3m(.)
#colnames(pair_tnst_mr)=c('name2','name1','MST.ij.bray')

project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

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
########### merge distance and ST ratio ####################
check_match = sort(pair_tnst_mr$name1)==sort(distance$col) 
all(check_match)

check_match = sort(paste(pair_tnst_mr$name1,pair_tnst_mr$name2))==sort(paste(distance$col,distance$row)) 
all(check_match)

pair_tnst_mr$tag=paste(pair_tnst_mr$name1,pair_tnst_mr$name2)
distance$tag=paste(distance$col,distance$row)
pair_tnst_mr = inner_join(pair_tnst_mr, distance,by=c('tag'))
pair_tnst_mr$Study=prefix
assign(paste0('pair_tnst_mr',prefix),pair_tnst_mr)

#################### Altitude (HKV) ########################
load('3_tnst/tnst_mr_Altitude (HKV)_OTU')
load('project_Altitude (HKV)_OTU')
tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
pair_tnst_mr = tnst_mr$index.pair

project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

distance = geo_dist(ENV) %>% col3m(.,dist_name = 'distance')
#################### elevation distance
distance.ele = ENV[,'Elevation'] %>% dist(.)/1000
distance.ele = col3m(distance.ele,dist_name = 'distance') 
distance$distance=distance$distance+distance.ele$distance

########### merge distance and ST ratio ####################
check_match = sort(pair_tnst_mr$name1)==sort(distance$col) 
all(check_match)

check_match = sort(paste(pair_tnst_mr$name1,pair_tnst_mr$name2))==sort(paste(distance$col,distance$row)) 
all(check_match)

pair_tnst_mr$tag=paste(pair_tnst_mr$name1,pair_tnst_mr$name2)
distance$tag=paste(distance$col,distance$row)
pair_tnst_mr = inner_join(pair_tnst_mr, distance,by=c('tag'))
pair_tnst_mr$Study=prefix
assign(paste0('pair_tnst_mr',prefix),pair_tnst_mr)
#################### Altitude (SNNR) ########################
load('3_tnst/tnst_mr_Altitude (SNNR)_OTU')
load('project_Altitude (SNNR)_OTU')

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
pair_tnst_mr = tnst_mr$index.pair
############ for test #########################
#pair_tnst_mr = reshape_col3m(pair_tnst_mr,value = 'MST.ij.bray',symmetric = T) %>% dist(.) %>% col3m(.)
#colnames(pair_tnst_mr)=c('name2','name1','MST.ij.bray')

project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

distance = geo_dist(ENV) %>% col3m(.,dist_name = 'distance')
#################### elevation distance
distance.ele = ENV[,'Elevation'] %>% dist(.)/1000
distance.ele = col3m(distance.ele,dist_name = 'distance') 
distance$distance=distance$distance+distance.ele$distance

########### merge distance and ST ratio ####################
check_match = sort(pair_tnst_mr$name1)==sort(distance$col) 
all(check_match)

check_match = sort(paste(pair_tnst_mr$name1,pair_tnst_mr$name2))==sort(paste(distance$col,distance$row)) 
all(check_match)

pair_tnst_mr$tag=paste(pair_tnst_mr$name1,pair_tnst_mr$name2)
distance$tag=paste(distance$col,distance$row)
pair_tnst_mr = inner_join(pair_tnst_mr, distance,by=c('tag'))
pair_tnst_mr$Study=prefix
assign(paste0('pair_tnst_mr',prefix),pair_tnst_mr)

pair_tnst_mr = Reduce(function(x, y) rbind(x,y),list(`pair_tnst_mrAltitude (HKV)`,`pair_tnst_mrLatitude (North America)`,`pair_tnst_mrAltitude (SNNR)`))
pair_tnst_mr$Study = factor(pair_tnst_mr$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
pair_tnst_mr$ln.distance = log1p(pair_tnst_mr$distance)
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
# log.distance
syc = scale_y_continuous(labels = scales::number_format(accuracy = 0.01,decimal.mark = '.'))
sxc = scale_x_continuous(labels = scales::number_format(accuracy = 0.1,decimal.mark = '.'))

SDecay <- ggplot(pair_tnst_mr,aes(x=ln.distance,y=MST.ij.bray))+geom_point(alpha = 0.5,stroke = 0,size=1.5,shape=16,aes(color=Study))+geom_smooth(method='lm',aes(color=Study),se = F)
plot = SDecay+facet_grid(~Study,scales = 'free')#+lm_Site+Ps_Sitescm_Site+sfm_Site
theme_set(theme_bw()+theme(legend.position = "none"))
plot+scale_color_economist()+syc+sxc
#plot+scale_color_viridis(discrete=T,alpha=0.9) +ylim(c(0,1))
# +stat_boxplot(geom='errorbar',width=0.15)
ggsave(filename = paste0(folder,'/','Decay_',note,'_MST.pdf'),dpi=900,width=19.5,height=7.5,units='cm')
ggsave(filename = paste0(folder,'/','Decay_',note,'_MST.jpg'),dpi=900,width=19.5,height=7.5,units='cm')


