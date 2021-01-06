library(dplyr)
library(vegan)
library(phyloseq)
library(picante)
library(ieggr)
source('0_function.R')
### Meta information  
meta <- read.csv('0_data/merge_mapping.csv',header = TRUE,row.names = 1,sep=',',na.strings = "NA")

######################################## HKV
prefix = 'Altitude (HKV)'
dt <- read.csv('0_data/hawaii_16S_otu_table.tsv',header = TRUE,row.names = 1,sep='\t')
taxa <- read.csv('0_data/hawaii_16S_otu_taxonomy.tsv',header = TRUE,sep='\t',row.names = 1,na.strings = "NA")
tree <- read_tree('0_data/hawaii_16S_otu_tree.nwk') %>% root(.,1)

# remove non-bacteria taxa
taxa$kingdom %>% summary
taxa = taxa[taxa$kingdom == 'Bacteria',] %>% droplevels(.)
taxa$phylum %>% summary(.)
# check chloroplast
taxa$class [taxa$phylum == 'Cyanobacteria'] %>% droplevels(.) %>% summary(.)

sum(rownames(dt) %in% rownames(taxa))

################### singleton filter ################
singleton_frequency=1
singleton_index = c(rowSums(dt)<=singleton_frequency)
print(paste0("there are ",sum(singleton_index)," singeltons less than frequency ",singleton_frequency))
dt = dt[!singleton_index,]

###############################################
##### phylo object
{
OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
TAX <- taxa %>% as.matrix(.) %>% tax_table(.) 
ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
project = phyloseq(OTU,TAX,ENV,tree)
}

##################################
########## data process (filter)
OTU <- otu_table(project)
ENV = sample_data(project)
### remove site-singleton
OTU = clean_cut(OTU,map=ENV,group = 'Site',cut=1) %>% otu_table(.,taxa_are_rows = T)
OTU = OTU[,order(colnames(OTU))]
colnames(OTU)=rownames(ENV)

########## data process (resample)
set.seed(930827)
project = phyloseq(OTU,TAX,ENV,tree)
colSums(OTU) %>% sort(.) 
project = rarefy_even_depth(project,sample.size = 14000) 


#########
### other indexes and add alpha
{
  sample_data(project)$Shannon=diversity(t(otu_table(project)),index='shannon')
  sample_data(project)$Simpson=diversity(t(otu_table(project)),index='simpson')
  sample_data(project)$Richness=specnumber(t(otu_table(project)))
  PD = pd(t(otu_table(project)), phy_tree(project), include.root = F)
  sample_data(project)$PD = PD$PD 
}
##########
note='OTU'
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 

######################################## SNNR
prefix = 'Altitude (SNNR)'
dt <- read.csv('0_data/Shennongjia_16S_otu_table.csv',header = TRUE,row.names = 1,sep=',')
taxa <- read.csv('0_data/Shennongjia_16S_otu_taxanomy.csv',header = TRUE,sep=',',row.names = 1,na.strings = "NA")
tree <- read_tree('0_data/Shennongjia_16S_otu_tree.nwk') %>% root(.,1)

# remove non-bacteria taxa
taxa$phylum %>% summary
# Euryarchaeota Crenarchaeota Cyanobacteria/Chloroplast
non_bacteria=c('Euryarchaeota', 'Crenarchaeota', 'Cyanobacteria/Chloroplast')
taxa = taxa[taxa$phylum != non_bacteria,] %>% droplevels(.)

sum(rownames(dt) %in% rownames(taxa))

################### singleton filter ################
singleton_frequency=1
singleton_index = c(rowSums(dt)<=singleton_frequency)
print(paste0("there are ",sum(singleton_index)," singeltons less than frequency ",singleton_frequency))
dt = dt[!singleton_index,]

###############################################
##### phylo object
{
  OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
  TAX <- taxa %>% as.matrix(.) %>% tax_table(.) 
  ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
  project = phyloseq(OTU,TAX,ENV,tree)
}

##################################
########## data process (filter)
OTU <- otu_table(project)
ENV = sample_data(project)
### remove site-singleton
OTU = clean_cut(OTU,map=ENV,group = 'Site',cut=1) %>% otu_table(.,taxa_are_rows = T)
OTU = OTU[,order(colnames(OTU))]
colnames(OTU)=rownames(ENV)

########## data process (resample)
set.seed(930827)
project = phyloseq(OTU,TAX,ENV,tree)
colSums(OTU) %>% sort(.) 
project = rarefy_even_depth(project,sample.size = 14000) 

#########
### other indexes and add alpha
{
  sample_data(project)$Shannon=diversity(t(otu_table(project)),index='shannon')
  sample_data(project)$Simpson=diversity(t(otu_table(project)),index='simpson')
  sample_data(project)$Richness=specnumber(t(otu_table(project)))
  PD = pd(t(otu_table(project)), phy_tree(project), include.root = F)
  sample_data(project)$PD = PD$PD 
}
##########
note='OTU'
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note))

######################################## NA
prefix = 'Latitude (North America)'
dt <- read.csv('0_data/NA_16S_otu_table.txt',header = TRUE,row.names = 1,sep='\t')
taxa <- read.csv('0_data/NA_16S_otu_taxonomy.tsv',header = TRUE,sep='\t',row.names = 1,na.strings = "NA")
tree <- read_tree('0_data/NA_16S_otu_tree.nwk') %>% root(.,1)

# remove non-bacteria taxa
taxa$kingdom %>% summary
taxa = taxa[taxa$kingdom == 'Bacteria',] %>% droplevels(.)
taxa$phylum %>% summary(.)
# check chloroplast
taxa$class [taxa$phylum == 'Cyanobacteria'] %>% droplevels(.) %>% summary(.)

sum(rownames(dt) %in% rownames(taxa))

################### singleton filter ################
singleton_frequency=1
singleton_index = c(rowSums(dt)<=singleton_frequency)
print(paste0("there are ",sum(singleton_index)," singeltons less than frequency ",singleton_frequency))
dt = dt[!singleton_index,]

###############################################
##### phylo object
{
  OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
  TAX <- taxa %>% as.matrix(.) %>% tax_table(.) 
  ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
  project = phyloseq(OTU,TAX,ENV,tree)
}

##################################
########## data process (filter)
OTU <- otu_table(project)
ENV = sample_data(project)
### remove site-singleton
OTU = clean_cut(OTU,map=ENV,group = 'Site',cut=1) %>% otu_table(.,taxa_are_rows = T)
OTU = OTU[,order(colnames(OTU))]
colnames(OTU)=rownames(ENV)

########## data process (resample)
set.seed(930827)
project = phyloseq(OTU,TAX,ENV,tree)
colSums(OTU) %>% sort(.) 
project = rarefy_even_depth(project,sample.size = 14000) 


#########
### other indexes and add alpha
{
  sample_data(project)$Shannon=diversity(t(otu_table(project)),index='shannon')
  sample_data(project)$Simpson=diversity(t(otu_table(project)),index='simpson')
  sample_data(project)$Richness=specnumber(t(otu_table(project)))
  PD = pd(t(otu_table(project)), phy_tree(project), include.root = F)
  sample_data(project)$PD = PD$PD 
}
##########
note='OTU'
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 

