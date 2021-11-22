source('0_function.R')

#################### save dir #######################
folder = paste0('8_non-para_stat')
dir.create(folder)
######################################
pair_adonis = function(OTU,ENV,Bray,Site){
  ## total
  result_dt=data.frame()
  result = adonis(Bray~Site, data=ENV, permutations = 999)
  result_dt[1,c('Site_comparision','R2','P')]=c('Overall',result$aov.tab[1,c('R2','Pr(>F)')])
  ## and then generate pair
  pair = levels(Site)
  pair_combn = combn(pair,2)
  Site_pair=apply(pair_combn,2,function(x)paste0(x[1],'-',x[2]))
  ## calculate pair
  for (i in 1:ncol(pair_combn)) {
    sample_index = colnames(OTU)[Site %in% pair_combn[,i]]
    Dist = dist_subset(Bray,sample_index)
    result = adonis(Dist~Site, data=ENV[sample_index,], permutations = 999)
    result_dt[i+1,] =c(Site_pair[i],result$aov.tab[1,c('R2','Pr(>F)')])
  }
  #result_dt$P = p.adjust(result_dt$P,method = 'holm')
  result_dt
}
pair_anoism = function(OTU,Bray,Site){
  ## total
  result_dt=data.frame()
  result = anosim(Bray,Site, permutations = 999)
  result_dt[1,c('Site_comparision','ANOSIM-R','P')]=c('Overall',result$statistic,result$signif)
  ## and then generate pair
  pair = levels(Site)
  pair_combn = combn(pair,2)
  Site_pair=apply(pair_combn,2,function(x)paste0(x[1],'-',x[2]))
  ## calculate pair
  for (i in 1:ncol(pair_combn)) {
    sample_index = colnames(OTU)[Site %in% pair_combn[,i]]
    site_index = Site[Site %in% pair_combn[,i]]
    Dist = dist_subset(Bray,sample_index)
    result = anosim(Dist,site_index, permutations = 999)
    result_dt[i+1,] =c(Site_pair[i],result$statistic,result$signif)
  }
  #result_dt$P = p.adjust(result_dt$P,method = 'holm')
  result_dt$`ANOSIM-R`=as.numeric(result_dt$`ANOSIM-R`)
  result_dt
}
pair_mrpp = function(OTU,Bray,Site){
  ## total
  result_dt=data.frame()
  result = mrpp(Bray,Site, permutations = 999)
  result_dt[1,c('Site_comparision','MRPP-A','P')]=c('Overall',result$A,result$Pvalue)
  ## and then generate pair
  pair = levels(Site)
  pair_combn = combn(pair,2)
  Site_pair=apply(pair_combn,2,function(x)paste0(x[1],'-',x[2]))
  ## calculate pair
  for (i in 1:ncol(pair_combn)) {
    sample_index = colnames(OTU)[Site %in% pair_combn[,i]]
    site_index = Site[Site %in% pair_combn[,i]]
    Dist = dist_subset(Bray,sample_index)
    result = mrpp(Dist,site_index, permutations = 999)
    result_dt[i+1,] =c(Site_pair[i],result$A,result$Pvalue)
  }
  #result_dt$P = p.adjust(result_dt$P,method = 'holm')
  result_dt$`MRPP-A`=as.numeric(result_dt$`MRPP-A`)
  result_dt
}
labelpstar <- function(x){
  if (x <= 0.001){x = "***"}
  else if (x <= 0.01){x = "**"}
  else if (x <= 0.05){x = "*"} 
  else {x = ""}
}
######################################
stat_pack = function(project){
  OTU <- otu_table(project)
  ENV <- sample_data(project) %>% data.frame(.)
  Site = ENV$Site 
  Bray = vegdist(t(OTU), method="bray", binary=FALSE)
  
  dt_adonis = pair_adonis(OTU,ENV,Bray,Site)
  dt_anoism = pair_anoism(OTU,Bray,Site)
  dt_mrpp = pair_mrpp(OTU,Bray,Site)
  dt_stat = Reduce(function(x, y) inner_join(x, y,by='Site_comparision'),list(dt_adonis,dt_anoism,dt_mrpp))
  dt_stat$R2=paste0(format(round(dt_stat$R2,2), nsmall = 2),sapply(dt_stat$P.x,function(x)labelpstar(x)))
  dt_stat$`ANOSIM-R`=paste0(format(round(dt_stat$`ANOSIM-R`,2), nsmall = 2),sapply(dt_stat$P.y,function(x)labelpstar(x)))
  dt_stat$`MRPP-A`=paste0(format(round(dt_stat$`MRPP-A`,2), nsmall = 2),sapply(dt_stat$P,function(x)labelpstar(x)))
  dt_stat
}
############ get the non-para stat ###################
######################################
level = 'phylum'
#######################################################
load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))
project = tax_glom(project,taxrank = level)
dt_stat = stat_pack(project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))
project = tax_glom(project,taxrank = level)
dt_stat = stat_pack(project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))
project = tax_glom(project,taxrank = level)
dt_stat = stat_pack(project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
######################################
level = 'class'
#######################################################
load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))
project = tax_glom(project,taxrank = level)
dt_stat = stat_pack(project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))
project = tax_glom(project,taxrank = level)
dt_stat = stat_pack(project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))
project = tax_glom(project,taxrank = level)
dt_stat = stat_pack(project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
level = 'reads'
#######################################################
load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))
dt_stat = stat_pack(project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))
dt_stat = stat_pack(project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))
dt_stat = stat_pack(project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
######################################
######################## GeoChip data ################
stat_pack <- function(OTU,project){
  ENV <- sample_data(project) %>% data.frame(.)
  Site = ENV$Site 
  Bray = vegdist(t(OTU), method="bray", binary=FALSE)
  dt_adonis = pair_adonis(OTU,ENV,Bray,Site)
  dt_anoism = pair_anoism(OTU,Bray,Site)
  dt_mrpp = pair_mrpp(OTU,Bray,Site)
  dt_stat = Reduce(function(x, y) inner_join(x, y,by='Site_comparision'),list(dt_adonis,dt_anoism,dt_mrpp))
  dt_stat$R2=paste0(format(round(dt_stat$R2,2), nsmall = 2),sapply(dt_stat$P.x,function(x)labelpstar(x)))
  dt_stat$`ANOSIM-R`=paste0(format(round(dt_stat$`ANOSIM-R`,2), nsmall = 2),sapply(dt_stat$P.y,function(x)labelpstar(x)))
  dt_stat$`MRPP-A`=paste0(format(round(dt_stat$`MRPP-A`,2), nsmall = 2),sapply(dt_stat$P,function(x)labelpstar(x)))
  dt_stat
}
######################################################
level = 'Gene_category'
#######################################################
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
taxa_list = tax_table(project)[,level]
dt = group_sum(otu_table(project),taxa_list,margin = 2)
dt_stat = stat_pack(dt,project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
taxa_list = tax_table(project)[,level]
dt = group_sum(otu_table(project),taxa_list,margin = 2)
dt_stat = stat_pack(dt,project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
taxa_list = tax_table(project)[,level]
dt = group_sum(otu_table(project),taxa_list,margin = 2)
dt_stat = stat_pack(dt,project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
######################################
level = 'Gene'
#######################################################
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
taxa_list = tax_table(project)[,level]
dt = group_sum(otu_table(project),taxa_list,margin = 2)
dt_stat = stat_pack(dt,project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
taxa_list = tax_table(project)[,level]
dt = group_sum(otu_table(project),taxa_list,margin = 2)
dt_stat = stat_pack(dt,project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
taxa_list = tax_table(project)[,level]
dt = group_sum(otu_table(project),taxa_list,margin = 2)
dt_stat = stat_pack(dt,project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
level = 'probes'
#######################################################
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
dt_stat = stat_pack(otu_table(project),project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
dt_stat = stat_pack(otu_table(project),project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)
#######################################################
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
dt_stat = stat_pack(otu_table(project),project)
write.csv(dt_stat,file = paste0(folder,'/Non-para_stat',prefix,'_',note,'_',level,'.csv'),row.names = FALSE)