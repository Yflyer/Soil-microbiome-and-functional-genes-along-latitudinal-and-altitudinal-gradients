library(NST)
library(iCAMP)
library(dplyr)
library(phyloseq)
source('0_function.R')

#########################################################################
# This code is to investigte the statisics of biodata 
############################################
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
###############################
folder = paste0(prefix,'_',note)
dir.create(folder)
########
data.type = 'gene'
threads = 20
######
OTU <- otu_table(project) %>% data.frame(.)
#TAX <- tax_table(project) %>% data.frame(.) 
ENV <- sample_data(project)

if (data.type=='gene') {
  OTU[OTU==0]=NA
  OTU=  apply(OTU,2,function(x) normalization(x))
  OTU = OTU* 1000 #0 %>% round(.)
}
OTU[is.na(OTU)]=0

Site = ENV$Site %>% as.data.frame(.)
rownames(Site)=sample_names(ENV)
colnames(Site)='Site'

snmout=iCAMP::snm.comm(comm = t(OTU), treat = NULL, 
                       rand=1000, alpha=0.05)
write.csv(snmout$stats,file = paste0(folder,'/',prefix,".NeutralModel.Stats.csv"))
write.csv(snmout$ratio.summary,file = paste0(folder,'/',prefix,".NeutralModel.TypeRatio.csv"))

snmout=iCAMP::snm.boot(comm = t(OTU), 
                       rand=1000, alpha=0.05)
write.csv(snmout$stats,file = paste0(folder,'/',prefix,".NeutralModel.boots.Stats.csv"))
write.csv(snmout$summary,file = paste0(folder,'/',prefix,".NeutralModel.summary.csv"))