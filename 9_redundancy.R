source('0_function.R')
###################################
#################### save dir #######################
folder = paste0('8_redundancy_stat')
dir.create(folder)
# This code is to investigate the redundancy metrice

p_value = c()
slope = c()
###########################################
load('project_Latitude (North America)_GeoChip')
p_gene_subset = get(paste0('project_',prefix,'_',note))
load('project_Latitude (North America)_OTU')
p_otu_subset = get(paste0('project_',prefix,'_',note))
###################################
p_gene_subset = subset_samples(p_gene_subset,sample_names(p_gene_subset) %in% sample_names(p_otu_subset)) 
p_otu_subset = subset_samples(p_otu_subset,sample_names(p_otu_subset) %in% sample_names(p_gene_subset)) 
################################
note='GeoChip'
OTU <- otu_table(p_gene_subset)
ENV <- sample_data(p_gene_subset)
###################################
Alpha = ENV$Richness
Shannon = ENV$Shannon
Dist = dist(t(OTU))
Bray = vegdist(t(OTU))
PC = pcoa(Bray, correction="none",rn=NULL) %>% .$vectors
pcoa1 = PC[,1]   
####################################
assign(paste0('Alpha_',note),Alpha)
assign(paste0('Shannon_',note),Shannon)
assign(paste0('Dist_',note),Dist)
assign(paste0('Bray_',note),Bray)
assign(paste0('pcoa1_',note),pcoa1)

#########################################
note='OTU'
OTU <- otu_table(p_otu_subset)
ENV <- sample_data(p_otu_subset)
######################################
Alpha = ENV$Richness
Shannon = ENV$Shannon
Dist = dist(t(OTU))
Bray = vegdist(t(OTU))
PC = pcoa(Bray, correction="none",rn=NULL) %>% .$vectors
pcoa1 = PC[,1]   
####################################
assign(paste0('Alpha_',note),Alpha)
assign(paste0('Shannon_',note),Shannon)
assign(paste0('Dist_',note),Dist)
assign(paste0('Bray_',note),Bray)
assign(paste0('pcoa1_',note),pcoa1)

Cor = cor.test(Alpha_GeoChip,Alpha_OTU,method = 'spearman',exact = F)
# extract parameters
slope = c(slope,Cor$estimate)
p_value =c(p_value,Cor$p.value)

Cor = cor.test(Shannon_GeoChip,Shannon_OTU,method = 'spearman',exact = F)
# extract parameters
slope = c(slope,Cor$estimate)
p_value =c(p_value,Cor$p.value)

Cor = cor.test(pcoa1_GeoChip,pcoa1_OTU,method = 'spearman',exact = F)
# extract parameters
slope = c(slope,Cor$estimate)
p_value =c(p_value,Cor$p.value)

Cor = mantel(Bray_GeoChip,Bray_OTU)
# extract parameter
slope = c(slope,Cor$statistic)
p_value =c(p_value,Cor$signif)

Cor = mantel(Dist_GeoChip,Dist_OTU)
# extract parameter
slope = c(slope,Cor$statistic)
p_value =c(p_value,Cor$signif)

###########################################
load('project_Altitude (SNNR)_GeoChip')
p_gene_subset = get(paste0('project_',prefix,'_',note))
load('project_Altitude (SNNR)_OTU')
p_otu_subset = get(paste0('project_',prefix,'_',note))
###################################
p_gene_subset = subset_samples(p_gene_subset,sample_names(p_gene_subset) %in% sample_names(p_otu_subset)) 
p_otu_subset = subset_samples(p_otu_subset,sample_names(p_otu_subset) %in% sample_names(p_gene_subset)) 
################################
note='GeoChip'
OTU <- otu_table(p_gene_subset)
ENV <- sample_data(p_gene_subset)
###################################
Alpha = ENV$Richness
Shannon = ENV$Shannon
Dist = dist(t(OTU))
Bray = vegdist(t(OTU))
PC = pcoa(Bray, correction="none",rn=NULL) %>% .$vectors
pcoa1 = PC[,1]   
####################################
assign(paste0('Alpha_',note),Alpha)
assign(paste0('Shannon_',note),Shannon)
assign(paste0('Dist_',note),Dist)
assign(paste0('Bray_',note),Bray)
assign(paste0('pcoa1_',note),pcoa1)

#########################################
note='OTU'
OTU <- otu_table(p_otu_subset)
ENV <- sample_data(p_otu_subset)
######################################
Alpha = ENV$Richness
Shannon = ENV$Shannon
Dist = dist(t(OTU))
Bray = vegdist(t(OTU))
PC = pcoa(Bray, correction="none",rn=NULL) %>% .$vectors
pcoa1 = PC[,1]   
####################################
assign(paste0('Alpha_',note),Alpha)
assign(paste0('Shannon_',note),Shannon)
assign(paste0('Dist_',note),Dist)
assign(paste0('Bray_',note),Bray)
assign(paste0('pcoa1_',note),pcoa1)

Cor = cor.test(Alpha_GeoChip,Alpha_OTU,method = 'spearman',exact = F)
# extract parameters
slope = c(slope,Cor$estimate)
p_value =c(p_value,Cor$p.value)

Cor = cor.test(Shannon_GeoChip,Shannon_OTU,method = 'spearman',exact = F)
# extract parameters
slope = c(slope,Cor$estimate)
p_value =c(p_value,Cor$p.value)

Cor = cor.test(pcoa1_GeoChip,pcoa1_OTU,method = 'spearman',exact = F)
# extract parameters
slope = c(slope,Cor$estimate)
p_value =c(p_value,Cor$p.value)

Cor = mantel(Bray_GeoChip,Bray_OTU)
# extract parameter
slope = c(slope,Cor$statistic)
p_value =c(p_value,Cor$signif)

Cor = mantel(Dist_GeoChip,Dist_OTU)
# extract parameter
slope = c(slope,Cor$statistic)
p_value =c(p_value,Cor$signif)


###########################################
load('project_Altitude (HKV)_GeoChip')
p_gene_subset = get(paste0('project_',prefix,'_',note))
load('project_Altitude (HKV)_OTU')
p_otu_subset = get(paste0('project_',prefix,'_',note))
###################################
p_gene_subset = subset_samples(p_gene_subset,sample_names(p_gene_subset) %in% sample_names(p_otu_subset)) 
p_otu_subset = subset_samples(p_otu_subset,sample_names(p_otu_subset) %in% sample_names(p_gene_subset)) 
################################
note='GeoChip'
OTU <- otu_table(p_gene_subset)
ENV <- sample_data(p_gene_subset)
###################################
Alpha = ENV$Richness
Shannon = ENV$Shannon
Dist = dist(t(OTU))
Bray = vegdist(t(OTU))
PC = pcoa(Bray, correction="none",rn=NULL) %>% .$vectors
pcoa1 = PC[,1]   
####################################
assign(paste0('Alpha_',note),Alpha)
assign(paste0('Shannon_',note),Shannon)
assign(paste0('Dist_',note),Dist)
assign(paste0('Bray_',note),Bray)
assign(paste0('pcoa1_',note),pcoa1)

#########################################
note='OTU'
OTU <- otu_table(p_otu_subset)
ENV <- sample_data(p_otu_subset)
######################################
Alpha = ENV$Richness
Shannon = ENV$Shannon
Dist = dist(t(OTU))
Bray = vegdist(t(OTU))
PC = pcoa(Bray, correction="none",rn=NULL) %>% .$vectors
pcoa1 = PC[,1]   
####################################
assign(paste0('Alpha_',note),Alpha)
assign(paste0('Shannon_',note),Shannon)
assign(paste0('Dist_',note),Dist)
assign(paste0('Bray_',note),Bray)
assign(paste0('pcoa1_',note),pcoa1)

Cor = cor.test(Alpha_GeoChip,Alpha_OTU,method = 'spearman',exact = F)
# extract parameters
slope = c(slope,Cor$estimate)
p_value =c(p_value,Cor$p.value)

Cor = cor.test(Shannon_GeoChip,Shannon_OTU,method = 'spearman',exact = F)
# extract parameters
slope = c(slope,Cor$estimate)
p_value =c(p_value,Cor$p.value)

Cor = cor.test(pcoa1_GeoChip,pcoa1_OTU,method = 'spearman',exact = F)
# extract parameters
slope = c(slope,Cor$estimate)
p_value =c(p_value,Cor$p.value)

Cor = mantel(Bray_GeoChip,Bray_OTU)
# extract parameter
slope = c(slope,Cor$statistic)
p_value =c(p_value,Cor$signif)

Cor = mantel(Dist_GeoChip,Dist_OTU)
# extract parameter
slope = c(slope,Cor$statistic)
p_value =c(p_value,Cor$signif)
