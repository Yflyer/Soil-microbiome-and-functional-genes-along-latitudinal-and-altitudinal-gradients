source('0_function.R')

#################### save dir #######################
folder = paste0('5_Nitrogen')
dir.create(folder)
#########################################################################
#########################################################################
library(stringr)
gene_filter_pack = function(project){
  project = subset_taxa(project,Gene_category %in% c('Nitrogen','nitrogen'))
  OTU = otu_table(project) #%>% .[tax_level!='Unknown',]
  TAX <- tax_table(project)
  ENV <- sample_data(project)  
  ### add unknown
  tax_level = apply(TAX[,'Lineage'],1,function(x) str_split(x,';')[[1]][1])
  tax_level[is.na(tax_level)]='Unknown'
  tax_level[tax_level=='']='Unknown'
  {
    gene_level = TAX[,'Gene']
    derep_level = paste(gene_level,tax_level) 
    names(derep_level)=rownames(gene_level)
    derep_level = derep_level[tax_level!='Unknown'] # remove unknown
  }
  ###
  OTU[OTU==0]=NA
  OTU = OTU[names(derep_level),] %>% scale(.)
  
  ##################
  Env_fun = ENV[,c('TN','NO3.N','NH4.N')] %>% as.matrix(.) # ,'CN_ratio'
  ##################
  unique_derep_level = unique(derep_level)
  
  Final_dt = data.frame()
  for (j in 1:ncol(Env_fun)) {
    p_value = c()
    slope = c()
    gene_ID = c()
    for (i in 1:nrow(OTU)) {
      # get abundance sum of OTU
      gene.abund = OTU[i,]
      if (sum(gene.abund>0,na.rm = T)>6) {
        # match
        match.index = !is.na(gene.abund)
        gene.abund = gene.abund[match.index]
        Fun_value = Env_fun[,j] %>% .[match.index] %>% scale(.)
        # cor test
        Cor = cor.test(gene.abund,Fun_value,method = 'spearman',exact = F)
        # extract parameters
        slope = c(slope,Cor$estimate)
        p_value =c(p_value,Cor$p.value)
        gene_ID = c(gene_ID,derep_level[names(derep_level)==rownames(OTU)[i]])
      }
    }
    result.lm = data.frame(gene_ID,slope,p_value) # ,gene_lineage
    # p adjust
    adjust.index = p.adjust(result.lm$p_value, method = "BH") < 0.05 # bonferroni
    if (sum(adjust.index)>0) {
      result.lm = result.lm[adjust.index,]
      result.lm[['ENV_fun']]=colnames(Env_fun)[j]
      Final_dt=rbind(result.lm,Final_dt)
    }
  }
  Final_dt$gene_name = sapply(Final_dt[,'gene_ID'],function(x) str_split(x,' ')[[1]][1])
  Final_dt$gene_taxa = sapply(Final_dt[,'gene_ID'],function(x) str_split(x,' ')[[1]][2])
  Final_dt
}
#############################start
############### GeoChip #############
dt_stat=group_stat=data.frame()
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note)) 
dt_stat = gene_filter_pack(project)
dt_stat$Study=prefix
write.csv(dt_stat,file = paste0(folder,'/Singleprobe_detail_',prefix,'_',note,'.csv'),row.names = FALSE)
group_stat = dt_stat %>% group_by(ENV_fun,gene_name) %>% summarise(n=n())
write.csv(group_stat,file = paste0(folder,'/Singleprobe_group-stat_',prefix,'_',note,'.csv'),row.names = FALSE)


dt_stat=group_stat=data.frame()
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note)) 
dt_stat = gene_filter_pack(project)
dt_stat$Study=prefix
write.csv(dt_stat,file = paste0(folder,'/Singleprobe_detail_',prefix,'_',note,'.csv'),row.names = FALSE)
group_stat = dt_stat %>% group_by(ENV_fun,gene_name) %>% summarise(n=n())
write.csv(group_stat,file = paste0(folder,'/Singleprobe_group-stat_',prefix,'_',note,'.csv'),row.names = FALSE)

dt_stat=group_stat=data.frame()
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note)) 
dt_stat = gene_filter_pack(project)
dt_stat$Study=prefix
write.csv(dt_stat,file = paste0(folder,'/Singleprobe_detail_',prefix,'_',note,'.csv'),row.names = FALSE)
group_stat = dt_stat %>% group_by(ENV_fun,gene_name) %>% summarise(n=n())
write.csv(group_stat,file = paste0(folder,'/Singleprobe_group-stat_',prefix,'_',note,'.csv'),row.names = FALSE)
