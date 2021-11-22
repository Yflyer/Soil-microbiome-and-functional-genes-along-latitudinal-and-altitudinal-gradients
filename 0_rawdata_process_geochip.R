source('0_function.R')
### Meta information  
meta <- read.csv('0_data/merge_mapping.csv',header = TRUE,row.names = 1,sep=',',na.strings = "NA")

########################################
prefix = 'Altitude (HKV)'
dt <- read.csv('0_data/hawaii_GeoChip_ln_mr.txt',header = TRUE,row.names = 1,sep='\t')
taxa <- select_if(dt,is.factor)
dt <- select_if(dt,is.numeric)
dt[is.na(dt)]=0

################### singleton filter ################
singleton_frequency=1
singleton_index = c(rowSums(dt>0)<=singleton_frequency)
print(paste0("there are ",sum(singleton_index)," singeltons less than frequency ",singleton_frequency))
dt = dt[!singleton_index,]
##### phylo object
{
OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
TAX <- taxa %>% GeoChip_TAX_reformat(.) %>% as.matrix(.) %>% tax_table(.)
ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
project = phyloseq(OTU,TAX,ENV)
}

#########
### other indexes and add alpha
{
  sample_data(project)$Shannon=diversity(t(otu_table(project)),index='shannon')
  sample_data(project)$Simpson=diversity(t(otu_table(project)),index='simpson')
  sample_data(project)$Richness=specnumber(t(otu_table(project)))
}
##########
note='GeoChip'
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 

########################################
prefix = 'Altitude (SNNR)'
dt <- read.csv('0_data/Shennongjia_GeoChip_ln_mr.txt',header = TRUE,row.names = 1,sep='\t')
taxa <- select_if(dt,is.factor)
dt <- select_if(dt,is.numeric) 
dt[is.na(dt)]=0

################### singleton filter ################
singleton_frequency=1
singleton_index = c(rowSums(dt>0)<=singleton_frequency)
print(paste0("there are ",sum(singleton_index)," singeltons less than frequency ",singleton_frequency))
dt = dt[!singleton_index,]
###############################################
##### phylo object
{
  OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
  TAX <- taxa %>% GeoChip_TAX_reformat(.) %>% as.matrix(.) %>% tax_table(.)
  ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
  project = phyloseq(OTU,TAX,ENV)
}
#########
### other indexes and add alpha
{
  sample_data(project)$Shannon=diversity(t(otu_table(project)),index='shannon')
  sample_data(project)$Simpson=diversity(t(otu_table(project)),index='simpson')
  sample_data(project)$Richness=specnumber(t(otu_table(project)))
}
##########
note='GeoChip'
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 

########################################
prefix = 'Latitude (North America)'
dt <- read.csv('0_data/NA_GeoChip_ln_mr.txt',header = TRUE,row.names = 1,sep='\t')
taxa <- select_if(dt,is.factor)
dt <- select_if(dt,is.numeric)
dt[is.na(dt)]=0

################### singleton filter ################
singleton_frequency=1
singleton_index = c(rowSums(dt>0)<=singleton_frequency)
print(paste0("there are ",sum(singleton_index)," singeltons less than frequency ",singleton_frequency))
dt = dt[!singleton_index,]
##### phylo object
{
  OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
  TAX <- taxa %>% GeoChip_TAX_reformat(.) %>% as.matrix(.) %>% tax_table(.)
  ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
  project = phyloseq(OTU,TAX,ENV)
}
#########
### other indexes and add alpha
{
  sample_data(project)$Shannon=diversity(t(otu_table(project)),index='shannon')
  sample_data(project)$Simpson=diversity(t(otu_table(project)),index='simpson')
  sample_data(project)$Richness=specnumber(t(otu_table(project)))
}
##########
note='GeoChip'
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 
