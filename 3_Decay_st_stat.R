library(geosphere)
library(ieggr)
source('0_function.R')
set.seed(999)

#################### save dir #######################
folder = paste0('3_tnst')
dir.create(folder)
# This code is to investigte the gene decay relation ship
#################### Latitude (NA) ############################
load('3_tnst/tnst_mr_Latitude (North America)_GeoChip')
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
pair_tnst_mr = tnst_mr$index.pair
m.dist = reshape_col3m(pair_tnst_mr,value = 'MST.ij.bray',symmetric = F) %>% as.dist(.)
#################################################
m.grad = (distHaversine(ENV['A100E',c('Lon','Lat')],ENV[,c('Lon','Lat')]) / 1000) %>% as.matrix(.)
rownames(m.grad)=rownames(ENV)
##################################################
Model.decay=tdecay(time=m.grad, dis = m.dist, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)
##################################################
# statistic r: 0.1142 Significance: 0.001
# R2.adj 0.0156019 p.perm 0.0049


##############################################
load('3_tnst/tnst_mr_Latitude (North America)_OTU')
load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
pair_tnst_mr = tnst_mr$index.pair
m.dist = reshape_col3m(pair_tnst_mr,value = 'MST.ij.bray',symmetric = F) %>% as.dist(.)
#################################################
m.grad = (distHaversine(ENV['A100E',c('Lon','Lat')],ENV[,c('Lon','Lat')]) / 1000) %>% as.matrix(.)
rownames(m.grad)=rownames(ENV)
##################################################
Model.decay=tdecay(time=m.grad, dis = m.dist, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)
# statistic r: 0.5546 Significance: 0.001
# R2.adj 0.295176 p.perm 0.0009

########################## SNNR #######################
load('3_tnst/tnst_mr_Altitude (SNNR)_GeoChip')
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
pair_tnst_mr = tnst_mr$index.pair
m.dist = reshape_col3m(pair_tnst_mr,value = 'MST.ij.bray',symmetric = F) %>% as.dist(.)
#################################################
m.grad = as.matrix(ENV[,'Elevation'])/1000
##################################################
Model.decay=tdecay(time=m.grad, dis = m.dist, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)
# R2.adj 0.000293 p.perm 0.301698
# statistic r: -0.04699 Significance: 0.8593

##############################################
load('3_tnst/tnst_mr_Altitude (SNNR)_OTU')
load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
pair_tnst_mr = tnst_mr$index.pair
m.dist = reshape_col3m(pair_tnst_mr,value = 'MST.ij.bray',symmetric = F) %>% as.dist(.)
#################################################
m.grad = as.matrix(ENV[,'Elevation'])/1000
##################################################
Model.decay=tdecay(time=m.grad, dis = m.dist, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)
# statistic r: 0.5652 Significance: 0.001
# R2.adj 0.2711979 p.perm 0.000999

########################## HKV #######################
load('3_tnst/tnst_mr_Altitude (HKV)_GeoChip')
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
pair_tnst_mr = tnst_mr$index.pair
m.dist = reshape_col3m(pair_tnst_mr,value = 'MST.ij.bray',symmetric = F) %>% as.dist(.)
#################################################
m.grad = as.matrix(ENV[,'Elevation'])/1000
##################################################
Model.decay=tdecay(time=m.grad, dis = m.dist, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)
# R2.adj 0.001304516 p.perm 0.5774226
# statistic r: 0.1576 Significance: 0.06
##############################################
load('3_tnst/tnst_mr_Altitude (HKV)_OTU')
load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) 

tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
pair_tnst_mr = tnst_mr$index.pair
m.dist = reshape_col3m(pair_tnst_mr,value = 'MST.ij.bray',symmetric = F) %>% as.dist(.)
#################################################
m.grad = as.matrix(ENV[,'Elevation'])/1000
##################################################
Model.decay=tdecay(time=m.grad, dis = m.dist, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)
# R2.adj 0.173167 p.perm 0.0009
# statistic r: 0.4021 Significance: 0.001



