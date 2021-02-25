library(vegan)
load('0_project.rdata')

#########################################################################
# This code is to investigte the statisics of biodata 
#############################start
project = project_SNJ
prefix='snj'

project = project_HA
prefix='ha'

project = project_NA
prefix='na'

rank_names(project)
#########################################
level = 'phylum'

### different level
project = tax_glom(project,taxrank = 'D1')
###########################

### different level
project = tax_glom(project,taxrank = 'D5')
level = 'genus'

#####
project = tax_glom(project,taxrank = 'D6')
level = 'gene'

###################
OTU <- otu_table(project)
TAX <- tax_table(project)
ENV <- sample_data(project)

##########################pair the result
# generate pair
treat = ENV$replicate_site
treat_name = 'replicate_site'

pair = levels(treat)
pair_combn = combn(pair,2)

pair_adonis = data.frame(matrix(ncol=length(pair),nrow=length(pair),dimnames = list(pair,pair)))
for (i in 1:ncol(pair_combn)) {
  # merge the pair subset
  pair_project_1 = subset_samples(project,replicate_site==pair_combn[1,i])
  pair_project_2 = subset_samples(project,replicate_site==pair_combn[2,i])
  pair_project = merge_phyloseq(pair_project_1,pair_project_2)
  # load the data from the merge pair data
  ENV <- sample_data(pair_project)
  OTU <- otu_table(pair_project)
  group = ENV[[treat_name]] %>% as.data.frame(.)
  colnames(group)='group'
  # get the beta and their adonis result
  jaccard = vegdist(t(OTU), method="jaccard", binary=FALSE)
  result = adonis(jaccard~group, data=group, permutations = 999)
  result_table = result$aov.tab
  pair_adonis[pair_combn[1,i],pair_combn[2,i]]=result_table[1,5]
}

write.csv(pair_adonis,paste0('5_stat_adonis_pair_probe_',prefix,'_',level,'.csv'))

pair_anoism = data.frame(matrix(ncol=length(pair),nrow=length(pair),dimnames = list(pair,pair)))
for (i in 1:ncol(pair_combn)) {
  # merge the pair subset
  pair_project_1 = subset_samples(project,replicate_site==pair_combn[1,i])
  pair_project_2 = subset_samples(project,replicate_site==pair_combn[2,i])
  pair_project = merge_phyloseq(pair_project_1,pair_project_2)
  # load the data from the merge pair data
  ENV <- sample_data(pair_project)
  OTU <- otu_table(pair_project)
  group = ENV[[treat_name]]
  # get the beta and their adonis result
  jaccard = vegdist(t(OTU), method="jaccard", binary=FALSE)
  result = anosim(jaccard,group, permutations = 999)
  pair_anoism[pair_combn[1,i],pair_combn[2,i]]=result$signif
}

write.csv(pair_anoism,paste0('5_stat_anoism_pair_probe_',prefix,'_',level,'.csv'))

pair_mrpp = data.frame(matrix(ncol=length(pair),nrow=length(pair),dimnames = list(pair,pair)))
for (i in 1:ncol(pair_combn)) {
  # merge the pair subset
  pair_project_1 = subset_samples(project,replicate_site==pair_combn[1,i])
  pair_project_2 = subset_samples(project,replicate_site==pair_combn[2,i])
  pair_project = merge_phyloseq(pair_project_1,pair_project_2)
  # load the data from the merge pair data
  ENV <- sample_data(pair_project)
  OTU <- otu_table(pair_project)
  group = ENV[[treat_name]]
  # get the beta and their adonis result
  jaccard = vegdist(t(OTU), method="jaccard", binary=FALSE)
  result = mrpp(jaccard,group, permutations = 999)
  pair_mrpp[pair_combn[1,i],pair_combn[2,i]]=result$Pvalue
}
write.csv(pair_mrpp,paste0('5_stat_mrpp_pair_probe_',prefix,'_',level,'.csv'))
head(OTU)
