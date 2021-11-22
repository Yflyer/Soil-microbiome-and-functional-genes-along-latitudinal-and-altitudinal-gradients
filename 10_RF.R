source('0_function.R')
library(randomForest)
library(rfUtilities)
library(rfPermute)
#####################################################
# need R 4.0 ENVIRONMENT
#################### save dir #######################
folder = paste0('10_RF')
dir.create(folder)
######################################
#######################################################
mergedt_pack<-function(project,ENV_var){
  OTU <- otu_table(project) 
  ENV <- sample_data(project) %>% data.frame(.) 
  
  Bray = vegdist(t(OTU), method="bray", binary=FALSE) %>% col3m(.,dist_name ='Bray')
  #################################################################
  Distance = geo_dist(ENV) %>% col3m(.,dist_name = 'Distance')
  Elevation = ENV[,'Elevation'] %>% dist(.)/1000
  Elevation = col3m(Elevation,dist_name = 'Elevation') 
  Distance$Distance=Distance$Distance+Elevation$Elevation
  Merge_dt = Reduce(function(x, y) inner_join(x, y,by=c('row','col')),list(Distance,Bray))
  
  ####################
  ENV_factor = ENV[, ENV_var] #%>% scale(.)
  for (i in colnames(ENV_factor)) {
    pair_factor = ENV_factor[,i] 
    names(pair_factor) = rownames(ENV)
    pair_factor = pair_factor %>% dist(.) %>% col3m(.,dist_name = i)
    Merge_dt = inner_join(Merge_dt, pair_factor,by=c('row','col'))
  }
  Merge_dt
}

rftest_pack <- function(Merge_dt,num.cores){
  Merge_dt <- select_if(Merge_dt,is.numeric)
  dt_rfP <- rfPermute(Bray~., data = Merge_dt, importance = TRUE, ntree = 500, nrep = 999, num.cores = num.cores)
  #######################################################
  specify_decimal(max(dt_rfP$rsq),3) %>% message(.)
  ##########################################################
  importance.scale <- data.frame(importance(dt_rfP, scale = TRUE), check.names = FALSE)
  pvalue.scale <-dt_rfP$pval[ , , 2]
  colnames(pvalue.scale)=c('p.value','node.p.value')
  rf_result = cbind(importance.scale,pvalue.scale)
  rf_result
}

############ get the random forest result ###################
######################################
ENV_var = c("pH","Rainfall","Temperature","NH4.N","NO3.N","TN","TC","CN_ratio","Plant_richness")
######################################
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
merge_dt =mergedt_pack(project,ENV_var)
rf_result = rftest_pack(merge_dt,6)
assign(paste0('rf_result',prefix,'_',note),rf_result)
write.csv(rf_result,file = paste0(folder,'/RF_',prefix,'_',note,'.csv'))

######################################
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
merge_dt =mergedt_pack(project,ENV_var)
rf_result = rftest_pack(merge_dt,6)
assign(paste0('rf_result',prefix,'_',note),rf_result)
write.csv(rf_result,file = paste0(folder,'/RF_',prefix,'_',note,'.csv'))

######################################
ENV_var = c("pH","Rainfall","Temperature","NH4.N","NO3.N","TN","TC","CN_ratio")
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
merge_dt =mergedt_pack(project,ENV_var)
rf_result = rftest_pack(merge_dt,6)
assign(paste0('rf_result',prefix,'_',note),rf_result)
write.csv(rf_result,file = paste0(folder,'/RF_',prefix,'_',note,'.csv'))

######################################
ENV_var = c("pH","Rainfall","Temperature","NH4.N","NO3.N","TN","TC","CN_ratio","Plant_richness")
######################################
load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))
merge_dt =mergedt_pack(project,ENV_var)
rf_result = rftest_pack(merge_dt,6)
assign(paste0('rf_result',prefix,'_',note),rf_result)
write.csv(rf_result,file = paste0(folder,'/RF_',prefix,'_',note,'.csv'))

######################################
load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))
merge_dt =mergedt_pack(project,ENV_var)
rf_result = rftest_pack(merge_dt,6)
assign(paste0('rf_result',prefix,'_',note),rf_result)
write.csv(rf_result,file = paste0(folder,'/RF_',prefix,'_',note,'.csv'))

######################################
ENV_var = c("pH","Rainfall","Temperature","NH4.N","NO3.N","TN","TC","CN_ratio")
load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))
merge_dt =mergedt_pack(project,ENV_var)
rf_result = rftest_pack(merge_dt,6)
assign(paste0('rf_result',prefix,'_',note),rf_result)
write.csv(rf_result,file = paste0(folder,'/RF_',prefix,'_',note,'.csv'))
