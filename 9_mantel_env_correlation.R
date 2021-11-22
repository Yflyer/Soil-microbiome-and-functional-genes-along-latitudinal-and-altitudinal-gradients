source('0_function.R')
#################### save dir #######################
folder = paste0('9_mantel')
dir.create(folder)
#########################################################################
############### GeoChip #############
library(parallel)
library(foreach)
library(doParallel)
library(doFuture)
registerDoFuture()
plan(multisession,workers=8)
#########################
mantel_pack <- function(otu_dis_list,gene_dis_list){
  result.dt = foreach (i=1:length(otu_dis_list),.combine = 'rbind') %dopar% {
    dis_stat=mantel(gene_dis_list[[i]],otu_dis_list[[i]],method="pearson",permutations=99999,na.rm = T)
    c(dis_stat$statistic,dis_stat$signif)
  }
  result.dt=data.frame(result.dt)
  rownames(result.dt) = names(otu_dis_list)
  colnames(result.dt)=c('Mantel_r','Mantel_P')
  result.dt$Adjust.P = p.adjust(result.dt$Mantel_P, method = "BH")
  result.dt
}

### we need to first fix gene level only change taxa level
level = 'Gene'
load('project_Altitude (HKV)_GeoChip')
load('project_Altitude (SNNR)_GeoChip')
load('project_Latitude (North America)_GeoChip')

level = 'class'
load('project_Altitude (HKV)_OTU')
load('project_Altitude (SNNR)_OTU')
load('project_Latitude (North America)_OTU')

project_gene = get(paste0('project_',prefix,'_',note))
gene_TAX = tax_table(project_gene)
gene_otu = otu_table(project_gene)

gene_factor = gene_TAX[,level] %>% unique(.)
### get gene and otu dist by taxon level
registerDoFuture()
plan(multisession,workers=12)
gene_dis_list = foreach (i=1:length(gene_factor)) %dopar% {
  lineage_index=gene_TAX[,level]==as.character(gene_factor[i])
  sub_otu = gene_otu[lineage_index,]
  if (nrow(sub_otu)>5 & all(colSums(sub_otu)>0)) {
    complete.cases(sub_otu)
    gene_dis = sub_otu %>% t(.) %>% vegdist()
    return(gene_dis)
  } else  {
    return(NA)
  }
}
names(gene_dis_list)=as.character(gene_factor)
sum(is.na(gene_dis_list))

####################
ENV_var = c("pH","Rainfall","Temperature","NH4.N","NO3.N","TN","TC","CN_ratio")  
ENV <- sample_data(project_gene) %>% data.frame(.)

env_dis_list = foreach (i=1:length(ENV_var)) %dopar% {
  env_dist = ENV[,ENV_var[i]] 
  names(env_dist) = rownames(ENV)
  env_dist = env_dist %>% scale(.) %>% dist(.)
}
names(env_dis_list)=ENV_var

result_table=data.frame()
for (i in 1:length(env_dis_list)) {
  result.dt= foreach (j=1:length(gene_dis_list),.combine = 'rbind') %dopar% {
    if (!is.na(env_dis_list[i]) & !is.na(gene_dis_list[j])) {
      dis_stat=mantel(gene_dis_list[[j]],env_dis_list[[i]],method="spearman",permutations=9999,na.rm = T)
      c(dis_stat$statistic,dis_stat$signif)
    } else {
      return(NA)
    }}
   result.dt=data.frame(result.dt)
   colnames(result.dt)=c('Mantel_r','Mantel_P')
   result.dt$ID = names(gene_dis_list)
   result.dt$Env=names(env_dis_list)[i]
   result_table=rbind(result_table,result.dt)
}

result_table$Significance='Correlated'
result_table$Significance[result_table$Mantel_P>0.05]='Uncorrelated'
result_table$Level=level
result_table$Method='comp'
result_table$Study=prefix
result_table=result_table[complete.cases(result_table),]

write.csv(result_table,file = paste0(folder,'/cor-comp-',level,'_',prefix,'.csv'),row.names = T)
prefix
