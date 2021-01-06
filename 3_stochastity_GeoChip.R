library(NST)
library(dplyr)
library(phyloseq)
source('0_function.R')

#########################################################################
# This code is to investigte the statisics of biodata 
###############################
folder = '3_tnst'
dir.create(folder)

############################################
load('project_Latitude (North America)_GeoChip')
project = project = get(paste0('project_',prefix,'_',note))

OTU <- otu_table(project) %>% data.frame(.)  
ENV <- sample_data(project)

for (i in 1:ncol(OTU)) {
  x=as.numeric(OTU[i,]) * 10000
  OTU[i,]=round(x)
}

Site = ENV$Site %>% as.data.frame(.)
rownames(Site)=sample_names(ENV)
colnames(Site)='Site'

print('NST start')
tnst_mr=tNST(comm=t(OTU), group=Site, dist.method="bray",
             abundance.weighted=TRUE, rand=100,
             nworker=80, null.model="PF", between.group=TRUE,
             SES=F, RC=F)

summary(tnst_mr$index.pair$MST.ij.bray)
assign(paste0('tnst_mr_',prefix,'_',note),tnst_mr) 
save(list = c(paste0('tnst_mr_',prefix,'_',note),'prefix','note'),file = paste0(folder,'/tnst_mr_',prefix,'_',note)) 
print('NST Finished')

############################################
load('project_Altitude (SNNR)_GeoChip')
project = project = get(paste0('project_',prefix,'_',note))

OTU <- otu_table(project) %>% data.frame(.) 
ENV <- sample_data(project)

for (i in 1:ncol(OTU)) {
  x=as.numeric(OTU[i,]) * 10000
  OTU[i,]=round(x)
}

Site = ENV$Site %>% as.data.frame(.)
rownames(Site)=sample_names(ENV)
colnames(Site)='Site'

print('NST start')
tnst_mr=tNST(comm=t(OTU), group=Site, dist.method="bray",
             abundance.weighted=TRUE, rand=100,
             nworker=80, null.model="PF", between.group=TRUE,
             SES=F, RC=F)

summary(tnst_mr$index.pair$MST.ij.bray)
assign(paste0('tnst_mr_',prefix,'_',note),tnst_mr) 
save(list = c(paste0('tnst_mr_',prefix,'_',note),'prefix','note'),file = paste0(folder,'/tnst_mr_',prefix,'_',note)) 
print('NST Finished')

############################################
load('project_Altitude (HKV)_GeoChip')
project = project = get(paste0('project_',prefix,'_',note))

OTU <- otu_table(project) %>% data.frame(.) 
ENV <- sample_data(project)

for (i in 1:ncol(OTU)) {
  x=as.numeric(OTU[i,]) * 10000
  OTU[i,]=round(x)
}

Site = ENV$Site %>% as.data.frame(.)
rownames(Site)=sample_names(ENV)
colnames(Site)='Site'

print('NST start')
tnst_mr=tNST(comm=t(OTU), group=Site, dist.method="bray",
             abundance.weighted=TRUE, rand=100,
             nworker=80, null.model="PF", between.group=TRUE,
             SES=F, RC=F)

summary(tnst_mr$index.pair$MST.ij.bray)
assign(paste0('tnst_mr_',prefix,'_',note),tnst_mr) 
save(list = c(paste0('tnst_mr_',prefix,'_',note),'prefix','note'),file = paste0(folder,'/tnst_mr_',prefix,'_',note))
print('NST Finished')

