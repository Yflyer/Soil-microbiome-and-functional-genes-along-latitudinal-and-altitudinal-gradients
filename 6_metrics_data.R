source('0_function.R')

#################### save dir #######################
folder = paste0('6_metrics_effect')
dir.create(folder)
####### load data ##########
load('project_Altitude (HKV)_GeoChip')
load('project_Altitude (SNNR)_GeoChip')
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note)) 

OTU=otu_table(project)
TAX=tax_table(project)
ENV=sample_data(project)
gene_list = TAX[,'Gene'] %>% as.character()

### filter
unique_gene = table(gene_list) %>% as.data.frame() %>% filter(Freq>6) %>% .[['gene_list']] %>% as.character()

### abundance
ab_dt = group_sum(OTU,gene_list,margin = 2) %>% .[unique_gene,]

### diversity
div_dt = list()
for (i in 1:length(unique_gene)) {
  div_index = OTU[gene_list== unique_gene[i],] %>% t(.) %>% diversity(.,index='shannon')
  div_dt[[i]] = div_index
}
div_dt=data.frame(div_dt) %>% t(.) %>% as.data.frame()
row.names(div_dt) = unique_gene
### composition
comp_dt = OTU [gene_list%in% unique_gene,] %>% data.frame()


### use different method to measure abundance, diversity and composition
# abundance
ab_Jaccard = vegdist(t(ab_dt),method = 'jaccard') %>% col3m(dist_name = 'Dist')
ab_Jaccard$Method='Jaccard'
ab_Sorsen = vegdist(t(ab_dt),method = 'bray', binary=T) %>% col3m(dist_name = 'Dist')
ab_Sorsen$Method='Sorensen'
ab_Bray = vegdist(t(ab_dt),method = 'bray') %>% col3m(dist_name = 'Dist')
ab_Bray$Method='Bray-Curits'
ab_Mori = vegdist(t(ab_dt),method = 'horn') %>% col3m(dist_name = 'Dist')
ab_Mori$Method='Morisita-Horn'

ab_result= Reduce(rbind,list(ab_Jaccard,ab_Sorsen,ab_Bray,ab_Mori))
ab_result$Metric='Abundance'
# diversity
div_Jaccard = vegdist(t(div_dt),method = 'jaccard') %>% col3m(dist_name = 'Dist')
div_Jaccard$Method='Jaccard'
div_Sorsen = vegdist(t(div_dt),method = 'bray', binary=T) %>% col3m(dist_name = 'Dist')
div_Sorsen$Method='Sorensen'
div_Bray = vegdist(t(div_dt),method = 'bray') %>% col3m(dist_name = 'Dist')
div_Bray$Method='Bray-Curits'
div_Mori = vegdist(t(div_dt),method = 'horn') %>% col3m(dist_name = 'Dist')
div_Mori$Method='Morisita-Horn'

div_result= Reduce(rbind,list(div_Jaccard,div_Sorsen,div_Bray,div_Mori))
div_result$Metric='Diversity'
# composition
comp_Jaccard = vegdist(t(comp_dt),method = 'jaccard') %>% col3m(dist_name = 'Dist')
comp_Jaccard$Method='Jaccard'
comp_Sorsen = vegdist(t(comp_dt),method = 'bray', binary=T) %>% col3m(dist_name = 'Dist')
comp_Sorsen$Method='Sorensen'
comp_Bray = vegdist(t(comp_dt),method = 'bray') %>% col3m(dist_name = 'Dist')
comp_Bray$Method='Bray-Curits'
comp_Mori = vegdist(t(comp_dt),method = 'horn') %>% col3m(dist_name = 'Dist')
comp_Mori$Method='Morisita-Horn'

comp_result= Reduce(rbind,list(comp_Jaccard,comp_Sorsen,comp_Bray,comp_Mori))
comp_result$Metric='Composition'

Final_result= Reduce(rbind,list(ab_result,div_result,comp_result))
Final_result$row.Site = sapply(Final_result$row,function(x) {ENV$Site[rownames(ENV)==x]})
Final_result$col.Site = sapply(Final_result$col,function(x) {ENV$Site[rownames(ENV)==x]})
Final_result$Comparison = 'Between-group'
Final_result$Comparison [Final_result$row.Site==Final_result$col.Site] = 'Within-group'
Final_result$Study = prefix
write.csv(Final_result,paste0(folder,'/Three-metrics_comparision_data_',prefix,'.csv'),row.names = T)
### other notes
message('Study ',prefix,' included ',length(unique(gene_list)),' genes (total of ',nrow(OTU),' probes; after qualified filter, ',length(unique_gene),' genes remained (total of ',nrow(comp_dt),' probes.')

