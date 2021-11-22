library(lavaan)  #for doing the CFA
#library(semPlot)  #for plotting your CFA
source('0_function.R')
#################### save dir #######################
folder = paste0('11_SEM')
dir.create(folder)
######################################
mergedt_pack=function(project,ENV_var,multifun_index){
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
  ENV_factor = ENV[,ENV_var ] #%>% scale(.)
  for (i in colnames(ENV_factor)) {
    pair_factor = ENV_factor[,i] 
    names(pair_factor) = rownames(ENV)
    pair_factor = pair_factor %>% dist(.) %>% col3m(.,dist_name = i)
    Merge_dt = inner_join(Merge_dt, pair_factor,by=c('row','col'))
  }
  ####################  %>% scale(.)
  Multifuncitonality = ENV[,multifun_index] %>% .[,complete.cases(t(.))] %>% dist(.) %>% col3m(.,dist_name = 'Multifuncitonality')
  Merge_dt = inner_join(Merge_dt, Multifuncitonality,by=c('row','col')) %>% select_if(.,is.numeric) %>% scale(.) %>% as.data.frame(.)
  Merge_dt
}
######################################
ENV_var = c("pH","Rainfall","Temperature") 
multifun_index = c("NH4.N","NO3.N","TN","TC","CN_ratio","Plant_richness") 
#####################
#### original design model:
#### the finally used model will be confirmed according to fitness measures but basically consistent with the original one
sem_Model <- '
  Temperature ~ Distance
  pH ~ Distance
  Rainfall ~ Distance
  Bray ~ Distance + pH + Temperature + rainfall
  Multifuncitonality ~ Bray + pH + Temperature + Rainfall
'
#####################
#######################################################
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
Merge_dt = mergedt_pack(project,ENV_var,multifun_index)
###################### finally used model
sem_Model <- '
  Bray ~  pH + Rainfall + Distance + Temperature##  
  Rainfall~ Temperature 
  Temperature ~ Distance
  #pH ~ Rainfall 
  Multifuncitonality ~ Bray+ Rainfall + Distance+ Temperature+pH#  #  
'
######################
sem.fit <- sem(model = sem_Model,data = Merge_dt)
resid(sem.fit, "cor")
modificationindices(sem.fit, minimum.value = 20)
sink(file = paste0(folder,'/sem_',prefix,note))
fitmeasures(sem.fit)[c('rmsea','cfi','tli')]
summary(sem.fit)
sink()
#partable(sem.fit)
#standardizedsolution(sem.fit)
write.csv(file = paste0(folder,'/sem_dt_',prefix,note,'.csv'),standardizedsolution(sem.fit))

load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))
Merge_dt = mergedt_pack(project,ENV_var,multifun_index)
###################### finally used model
sem_Model <- '
  Bray ~  pH + Rainfall + Distance + Temperature##  
  #pH ~ Distance + Rainfall + Temperature
  Temperature ~ Distance + Rainfall
  Multifuncitonality ~ Bray + pH + Rainfall + Temperature # + Distance
'
######################
sem.fit <- sem(model = sem_Model,data = Merge_dt)
resid(sem.fit, "cor")
modificationindices(sem.fit, minimum.value = 20)
sink(file = paste0(folder,'/sem_',prefix,note))
fitmeasures(sem.fit)[c('rmsea','cfi','tli')]
summary(sem.fit)
sink()
write.csv(file = paste0(folder,'/sem_dt_',prefix,note,'.csv'),standardizedsolution(sem.fit))

#######################################################
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
Merge_dt = mergedt_pack(project,ENV_var,multifun_index)
###################### finally used model
sem_Model <- '
  Bray ~  Rainfall + Distance+ Temperature##  
  Temperature ~ Distance
  #pH ~ Distance + Rainfall
  Rainfall ~ Temperature
  Multifuncitonality ~ Bray + Rainfall  + Temperature + Distance + pH#  + Temperature
'
######################
sem.fit <- sem(model = sem_Model,data = Merge_dt)
resid(sem.fit, "cor")
modificationindices(sem.fit, minimum.value = 20)
sink(file = paste0(folder,'/sem_',prefix,note))
fitmeasures(sem.fit)[c('rmsea','cfi','tli')]
summary(sem.fit)
sink()
write.csv(file = paste0(folder,'/sem_dt_',prefix,note,'.csv'),standardizedsolution(sem.fit))

load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))
Merge_dt = mergedt_pack(project,ENV_var,multifun_index)
###################### finally used model
sem_Model <- '
  Bray ~  pH + Rainfall + Temperature + Distance##  
  Rainfall ~ Distance + Temperature
#  pH ~ Distance  + Temperature
  Multifuncitonality ~ Bray + pH + Rainfall  + Distance # 
'
######################
sem.fit <- sem(model = sem_Model,data = Merge_dt)
resid(sem.fit, "cor")
modificationindices(sem.fit, minimum.value = 20)
sink(file = paste0(folder,'/sem_',prefix,note))
fitmeasures(sem.fit)[c('rmsea','cfi','tli')]
summary(sem.fit)
sink()
write.csv(file = paste0(folder,'/sem_dt_',prefix,note,'.csv'),standardizedsolution(sem.fit))

#######################################################
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
Merge_dt = mergedt_pack(project,ENV_var,multifun_index)
###################### finally used model
sem_Model <- '
  Bray ~  pH + Rainfall #+ Distance#  + Temperature
  #Temperature ~  Distance + Rainfall
  pH ~ Rainfall + Distance
  Rainfall ~ Distance + Temperature
  Multifuncitonality ~ Bray  + Rainfall  + Distance+ Temperature#+ pH+ Distance #  
'
######################
sem.fit <- sem(model = sem_Model,data = Merge_dt)
resid(sem.fit, "cor")
modificationindices(sem.fit, minimum.value = 20)
sink(file = paste0(folder,'/sem_',prefix,note))
fitmeasures(sem.fit)[c('rmsea','cfi','tli')]
summary(sem.fit)
sink()
write.csv(file = paste0(folder,'/sem_dt_',prefix,note,'.csv'),standardizedsolution(sem.fit))

load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))
Merge_dt = mergedt_pack(project,ENV_var,multifun_index)
###################### finally used model
sem_Model <- '
  Bray ~  Rainfall  + pH + Distance + Temperature
  Temperature ~ Distance + Rainfall
  pH ~ Distance + Rainfall
  #Rainfall ~ Distance + Temperature
  Multifuncitonality ~ Bray + Rainfall  + Temperature + Distance
'
######################
sem.fit <- sem(model = sem_Model,data = Merge_dt)
resid(sem.fit, "cor")
modificationindices(sem.fit, minimum.value = 20)
sink(file = paste0(folder,'/sem_',prefix,note))
fitmeasures(sem.fit)[c('rmsea','cfi','tli')]
summary(sem.fit)
sink()
write.csv(file = paste0(folder,'/sem_dt_',prefix,note,'.csv'),standardizedsolution(sem.fit))

