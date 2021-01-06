library(ieggr)
source('0_function.R')

#################### save dir #######################
folder = paste0('2_Decay')
dir.create(folder)
#########################################################################
# This code is to investigte the gene distribution relation ship

######################################
######################################

########################## NA #######################
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) #%>% as.data.frame(.) 

bray = vegdist(t(OTU), method="bray", binary=FALSE)
#################################################
m.grad = (distHaversine(ENV['A100E',c('Lon','Lat')],ENV[,c('Lon','Lat')]) / 1000) %>% as.matrix(.)
rownames(m.grad)=rownames(ENV)
##################################################
Model.decay=tdecay(time=m.grad, comm = t(OTU), dist.method = "bray", abundance.weighted = T, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary

write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)

# R2.adj 0.02221462 p.perm 0.007992008
# statistic r: 0.1232 Significance: 0.001

##############################################
load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
bray = vegdist(t(OTU), method="bray", binary=FALSE)
##################################################
Model.decay=tdecay(time=m.grad, comm = t(OTU), dist.method = "bray", abundance.weighted = T, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)

# statistic r: 0.6951 Significance: 0.001
# R2.adj 0.3761369 p.perm 0.0009

########################## SNNR #######################
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) #%>% as.data.frame(.) 
bray = vegdist(t(OTU), method="bray", binary=FALSE)
#################################################
m.grad = as.matrix(ENV[,'Elevation'])/1000
##################################################
Model.decay=tdecay(time=m.grad, comm = t(OTU), dist.method = "bray", abundance.weighted = T, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)

# R2.adj -0.0004267441 p.perm 0.8251748
# statistic r: -0.0404 Significance: 0.695
##############################################
load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
bray = vegdist(t(OTU), method="bray", binary=FALSE)

Model.decay=tdecay(time=m.grad, comm = t(OTU), dist.method = "bray", abundance.weighted = T, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)

# statistic r: 0.5652 Significance: 0.001
# R2.adj 0.2711979 p.perm 0.000999
 
########################## HKV #######################
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) #%>% as.data.frame(.) 
bray = vegdist(t(OTU), method="bray", binary=FALSE)
#################################################
m.grad= as.matrix(ENV[,'Elevation'])/1000
##################################################
Model.decay=tdecay(time=m.grad, comm = t(OTU), dist.method = "bray", abundance.weighted = T, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)

# R2.adj 0.001304516 p.perm 0.5774226
# statistic r: 0.1576 Significance: 0.06
##############################################
load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))
OTU <- otu_table(project)
ENV <- sample_data(project) #%>% as.data.frame(.) 
bray = vegdist(t(OTU), method="bray", binary=FALSE)
#################################################
m.grad= as.matrix(ENV[,'Elevation'])/1000

##################################################
Model.decay=tdecay(time=m.grad, comm = t(OTU), dist.method = "bray", abundance.weighted = T, treat = NULL, perm.test = T, boot.strap = T, rand = 1000, Dmax = 1)
Model.decay$summary
write.csv(file = paste0(folder,'/Decay_stat_',prefix,note,'.csv'),Model.decay$summary)

# R2.adj 0.173167 p.perm 0.0009
# statistic r: 0.4021 Significance: 0.001



