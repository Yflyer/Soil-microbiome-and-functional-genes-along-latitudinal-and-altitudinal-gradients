library(dplyr)
library(vegan)
library(phyloseq)
library(ggplot2)
source('0_function.R')
load('0_project.rdata')


project = project_SNJ
prefix='snj'

project = project_HA
prefix='ha'

project = project_NA
prefix='na'

#########################################
level = 'phylum'

rank_names(project)
tax_table(project)
### different level
project = subset_taxa(project,D0=='Bacteria')
project = tax_glom(project,taxrank = 'D1')


sample_sums(project)
project = rarefy_even_depth(project,sample.size = min(sample_sums(project)))
project = rarefy_even_depth(project,sample.size = 20000)
###################
OTU <- otu_table(project)
TAX <- tax_table(project)
ENV <- sample_data(project)

abun_10_taxa = rowMeans(OTU) %>% sort(.,decreasing = TRUE) %>% .[1:10]
OTU = OTU[rownames(OTU)%in% names(abun_10_taxa),]
Taxa = TAX[rownames(TAX) %in% names(abun_10_taxa),'D1'] %>% factor(.)
site = ENV$replicate_site

rownames(TAX)

plot_data = data.frame()
for (i in 1:ncol(OTU)) {
  OTU_abundance = OTU[,i] %>% as.numeric()
  Sample = rep(colnames(OTU)[i],nrow(OTU)) 
  #OTU_abundance = apply(OTU,1,function(x) sum(x[site %in% levels(site)[i]])/summary(site)[i] )
  Site = rep(site[i],nrow(OTU)) 
  OOO = data.frame(OTU_abundance,Sample,Site,Taxa)
  plot_data = rbind(plot_data,OOO)
}


######################
PTs = theme(axis.text.x=element_blank(),
            axis.title=element_blank(),
            legend.position="bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            panel.spacing = unit(0.2, "lines"))
comp = ggplot(data = plot_data,aes(x = Sample, y = OTU_abundance,group = factor(Taxa)))
barset = geom_bar(aes(fill = factor(Taxa)),stat='identity',width=0.6,color='black')
area =  geom_area(aes(fill = factor(Taxa)),stat = "identity",alpha = 0.5)
line = geom_line(aes(colour = factor(Taxa)), position = "stack")
scm = scale_color_manual(values = c('orangered','#43a2ca','#018571','#3C5488','#f4a582','slategray','#80cdc1','#ca0020','#b2abd2','#5e3c99'))
sfm = scale_fill_manual(values = c('orangered','#43a2ca','#018571','#3C5488','#f4a582','slategray','#80cdc1','#ca0020','#b2abd2','#5e3c99'))
comp+PTs +line+barset+area +sfm +facet_grid(~Site,scales = "free", space = "free")
ggsave(filename = paste0('7_comp_',prefix,'.pdf'),dpi=900,width=14,height=19,units='cm')

