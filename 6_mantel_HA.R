library(vegan)
source('0_function.R')

#################################
#################### save dir #######################
folder = paste0('6_multi.fun')
dir.create(folder)
######################################
load('project_Altitude (HKV)_GeoChip')
load('project_NA_GeoChip')
load('project_Altitude (SNNR)_GeoChip')
######################################
load('project_HA_OTU')
load('project_NA_OTU')
load('project_SNJ_OTU')
###########################################
project = get(paste0('project_',prefix,'_',note)) 

OTU <- otu_table(project)
ENV <- sample_data(project)
######################################################
beta_tax=vegdist(t(OTU),method="bray",binary=FALSE)

#all = c('Lat','Lon','Rainfall','NH4.N','NO3.N','N','Total.Organic.Carbon','CN_ratio')
#############################################################
# partial mantel
colnames(ENV)
############################################################
all = c('Elevation','Rainfall','Water.potential','NH4.N','NO3.N','N','Total.Organic.Carbon','P','CN_ratio')
geo.list=c('Elevation')
mulfun.list=c('NH4.N','NO3.N','N','Total.Organic.Carbon','P','CN_ratio')
clm.list = c('Water.potential','Rainfall')
partial.geo.list=all[!all %in% geo.list]
partial.mulfun.list=all[!all %in% mulfun.list]
partial.clm.list =all[!all %in% clm.list]
############################################################
env_all = ENV[,c(all)] %>% scale(.) %>% vegdist(.,"euclid")

geographic_distance=ENV[,geo.list] %>% scale(.) %>% vegdist(.,"euclid")
Soil_multifunctionality = ENV[,mulfun.list] %>% scale(.) %>% vegdist(.,"euclid")
climate_condition = ENV[,clm.list] %>% scale(.) %>% vegdist(.,"euclid")

partial.geo = ENV[,partial.geo.list] %>% scale(.) %>% vegdist(.,"euclid")
partial.mulfun = ENV[,partial.mulfun.list] %>% scale(.) %>% vegdist(.,"euclid")
partial.clm = ENV[,partial.clm.list] %>% scale(.) %>% vegdist(.,"euclid")
################# inner correlation ###################
mantel(geographic_distance,climate_condition)
mantel(Soil_multifunctionality,climate_condition)
mantel(Soil_multifunctionality,geographic_distance)
############################################################
mantel.partial(beta_tax,geographic_distance,partial.geo)
mantel.partial(beta_tax,climate_condition,partial.clm)
mantel.partial(beta_tax,Soil_multifunctionality,partial.mulfun)
###############################################################
### mantel table on single var
mantel_all = mantel(beta_tax,env_all)
mantel_r = mantel_all$statistic
mantel_p = mantel_all$signif
for (i in 1:length(all)) {
  mantel_result = ENV[,all[i]] %>% scale(.) %>% vegdist(.,"euclid") %>% mantel(beta_tax,.)
  #partial_result = ENV[,colnames(ENV)[i]] %>% scale(.) %>% vegdist(.,"euclid") %>% mantel(beta_tax,.)
  mantel_r[i+1] = mantel_result$statistic
  mantel_p[i+1] = mantel_result$signif
}
variables = c('All',all)
mantel_dt = data.frame(variables,mantel_r,mantel_p)

write.csv(mantel_dt,file = paste0(folder,'/4_',prefix,note,'_mantel_mr.csv'),row.names = FALSE)