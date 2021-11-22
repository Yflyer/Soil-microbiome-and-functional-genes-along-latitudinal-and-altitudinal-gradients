source('0_function.R')
set.seed(999)

#################### save dir #######################
folder = paste0('3_tnst')
dir.create(folder)
# This code is to investigate the gene decay relation ship
#############################################################
RC_DEP = function(x){
  if (sum(x==(-1))>0.5*length(x)) {
    DEP='Homogenous selection'
    P.DEP=sum(x==(-1))/length(x)
  } else if (sum(x==(1))>0.5*length(x)) {
    DEP='Dispersal limitation'
    P.DEP=sum(x==(1))/length(x)
  } else {
    DEP='Drift'
    P.DEP=1-(sum(x==(-1))+sum(x==(1)))/length(x)
  }
  c(DEP,P.DEP)
}
merge_models_result=function(prefix,snm_dt,tnst){
  pair_tnst = tnst$index.pair
  group_tnst = tnst$index.pair.grp
  between_tnst = tnst$index.pair.between
  snm_result = snm_dt['mean',c('Neutral.uw', 'Below.uw',  'Above.uw','Neutral.wt', 'Below.wt',  'Above.wt')]
  # mean inner ST
  inner.mean.ST=mean(group_tnst$ST.ij.ruzicka)
  # mean inner NST
  inner.mean.NST=mean(group_tnst$NST.ij.ruzicka)
  # stdev inner NST
  inner.sd.NST=sd(group_tnst$NST.ij.ruzicka)
  # mean between ST
  between.mean.ST=mean(between_tnst$ST.ij.ruzicka)
  # mean between NST
  between.mean.NST=mean(between_tnst$NST.ij.ruzicka)
  # stdev between NST
  between.sd.NST=sd(group_tnst$NST.ij.ruzicka)
  ### add RC
  index = paste0(pair_tnst$name1,pair_tnst$name2) %in% paste0(group_tnst$name1,group_tnst$name2)
  ### assign rc to group data
  group_tnst$RC.ruzicka=pair_tnst$RC.ruzicka[index]
  inner.dep = RC_DEP(group_tnst$RC.ruzicka)
  ### outer group
  between.dep = RC_DEP(pair_tnst$RC.ruzicka[!index])
  ### merge result
  models_result= c(prefix,inner.mean.ST,inner.mean.NST,inner.sd.NST,inner.dep,between.mean.ST,between.mean.NST,between.sd.NST,between.dep,snm_result) %>% unlist()
}
####
dt_models=matrix(nrow=6,ncol=20)
#################### Latitude (NA) ############################
load(file = '3_tnst/Latitude (North America)_GeoChip/tnst_Latitude (North America)_GeoChip')
tnst = get(paste0('tnst_',prefix,'_',note))

# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)
snm_stat = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.Stats.csv'),row.names = 1)
snm_stat = snm_stat %>% filter(treatment.id=='All')%>%select(m,Rsqr,AIC) 

models_result=merge_models_result(prefix,snm_dt,tnst)
dt_models[1,]=unlist(c(models_result,snm_stat))


#################### SNNR ############################
load(file = '3_tnst/Altitude (SNNR)_GeoChip/tnst_Altitude (SNNR)_GeoChip')
tnst = get(paste0('tnst_',prefix,'_',note))

# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)
snm_stat = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.Stats.csv'),row.names = 1)
snm_stat = snm_stat %>% filter(treatment.id=='All')%>%select(m,Rsqr,AIC) 

models_result=merge_models_result(prefix,snm_dt,tnst)
dt_models[2,]=unlist(c(models_result,snm_stat))

#################### HKV ############################
load(file = '3_tnst/Altitude (HKV)_GeoChip/tnst_Altitude (HKV)_GeoChip')
tnst = get(paste0('tnst_',prefix,'_',note))

# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)
snm_stat = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.Stats.csv'),row.names = 1)
snm_stat = snm_stat %>% filter(treatment.id=='All')%>%select(m,Rsqr,AIC) 

models_result=merge_models_result(prefix,snm_dt,tnst)
dt_models[3,]=unlist(c(models_result,snm_stat))

#################### Latitude (NA) ############################
load(file = '3_tnst/Latitude (North America)_OTU/tnst_Latitude (North America)_OTU')
tnst = get(paste0('tnst_',prefix,'_',note))

# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)
snm_stat = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.Stats.csv'),row.names = 1)
snm_stat = snm_stat %>% filter(treatment.id=='All')%>%select(m,Rsqr,AIC) 

models_result=merge_models_result(prefix,snm_dt,tnst)
dt_models[4,]=unlist(c(models_result,snm_stat))


#################### SNNR ############################
load(file = '3_tnst/Altitude (SNNR)_OTU/tnst_Altitude (SNNR)_OTU')
tnst = get(paste0('tnst_',prefix,'_',note))

# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)
snm_stat = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.Stats.csv'),row.names = 1)
snm_stat = snm_stat %>% filter(treatment.id=='All')%>%select(m,Rsqr,AIC) 

models_result=merge_models_result(prefix,snm_dt,tnst)
dt_models[5,]=unlist(c(models_result,snm_stat))

#################### HKV ############################
load(file = '3_tnst/Altitude (HKV)_OTU/tnst_Altitude (HKV)_OTU')
tnst = get(paste0('tnst_',prefix,'_',note))

# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)
snm_stat = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.Stats.csv'),row.names = 1)
snm_stat = snm_stat %>% filter(treatment.id=='All')%>%select(m,Rsqr,AIC) 

models_result=merge_models_result(prefix,snm_dt,tnst)
dt_models[6,]=unlist(c(models_result,snm_stat))

colnames(dt_models)=c('Gradient','Inner-group.mean.ST','Inner-group.mean.NST','Inner-group.sd.NST','Inner-group.DEP','Inner-group.P.DEP','Between-group.mean.ST','Between-group.mean.NST','Between-group.sd.NST','Between-group.DEP','Between-group.P.DEP','Neutral.uw', 'Below.uw',  'Above.uw','Neutral.wt','Below.wt','Above.wt','m','Rsqr','AIC')
dt_models
write.csv(file = paste0(folder,'/Multi.models.summary.csv'),dt_models,)
