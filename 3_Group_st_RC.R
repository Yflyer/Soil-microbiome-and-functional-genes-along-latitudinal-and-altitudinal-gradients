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
merge_models_result=function(prefix,snm_dt,group_tnst_mr,pair_tnst_mr){
  snm_result = snm_dt['mean',c('Neutral.uw', 'Below.uw',  'Above.uw','Neutral.wt', 'Below.wt',  'Above.wt')]
  # mean ST
  mean.ST=mean(group_tnst_mr$ST.ij.bray)
  # mean NST
  mean.NST=mean(group_tnst_mr$NST.ij.bray)
  ### add RC
  index = paste0(pair_tnst_mr$name1,pair_tnst_mr$name2) %in% paste0(group_tnst_mr$name1,group_tnst_mr$name2)
  ### assign rc to group data
  group_tnst_mr$RC.bray=pair_tnst_mr$RC.bray[index]
  inner.dep = RC_DEP(group_tnst_mr$RC.bray)
  ### outer group
  outer.dep = RC_DEP(pair_tnst_mr$RC.bray[!index])
  ### merge result
  models_result= c(prefix,mean.ST,mean.NST,inner.dep,outer.dep,snm_result) %>% unlist()
}
####
dt_models=matrix(nrow=6,ncol=13)
#################### Latitude (NA) ############################
load(file = '3_tnst/tnst_mr_Latitude (North America)_GeoChip')
tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
group_tnst_mr = tnst_mr$index.pair.grp
pair_tnst_mr = tnst_mr$index.pair
# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)

models_result=merge_models_result(prefix,snm_dt,group_tnst_mr,pair_tnst_mr)
dt_models[1,]=models_result


#################### SNNR ############################
load(file = '3_tnst/tnst_mr_Altitude (SNNR)_GeoChip')
tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
group_tnst_mr = tnst_mr$index.pair.grp
pair_tnst_mr = tnst_mr$index.pair
# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)

models_result=merge_models_result(prefix,snm_dt,group_tnst_mr,pair_tnst_mr)
dt_models[2,]=models_result

#################### HKV ############################
load(file = '3_tnst/tnst_mr_Altitude (HKV)_GeoChip')
tnst_mr = get(paste0('tnst_mr_',prefix,'_',note))
group_tnst_mr = tnst_mr$index.pair.grp
pair_tnst_mr = tnst_mr$index.pair
# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)

models_result=merge_models_result(prefix,snm_dt,group_tnst_mr,pair_tnst_mr)
dt_models[3,]=models_result

#################### Latitude (NA) ############################
load(file = '3_tnst/tnst_Latitude (North America)_OTU')
tnst_mr = get(paste0('tnst_',prefix,'_',note))
group_tnst_mr = tnst_mr$index.pair.grp
pair_tnst_mr = tnst_mr$index.pair
# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)

models_result=merge_models_result(prefix,snm_dt,group_tnst_mr,pair_tnst_mr)
dt_models[4,]=models_result


#################### SNNR ############################
load(file = '3_tnst/tnst_Altitude (SNNR)_OTU')
tnst_mr = get(paste0('tnst_',prefix,'_',note))
group_tnst_mr = tnst_mr$index.pair.grp
pair_tnst_mr = tnst_mr$index.pair
# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)

models_result=merge_models_result(prefix,snm_dt,group_tnst_mr,pair_tnst_mr)
dt_models[5,]=models_result

#################### HKV ############################
load(file = '3_tnst/tnst_Altitude (HKV)_OTU')
tnst_mr = get(paste0('tnst_',prefix,'_',note))
group_tnst_mr = tnst_mr$index.pair.grp
pair_tnst_mr = tnst_mr$index.pair
# snm result
snm_dt = read.csv(file = paste0('3_tnst/',prefix,'_',note,'/',prefix,'.NeutralModel.TypeRatio.csv'),row.names = 2)

models_result=merge_models_result(prefix,snm_dt,group_tnst_mr,pair_tnst_mr)
dt_models[6,]=models_result

colnames(dt_models)=c('Gradient','Mean.ST','Mean.NST','Inner-group.DEP','Inner-group.P.DEP','Outer-group.DEP','Outer-group.P.DEP','Neutral.uw', 'Below.uw',  'Above.uw','Neutral.wt','Below.wt','Above.wt')
dt_models
write.csv(file = paste0(folder,'/Multi.models.summary.csv'),dt_models,)
