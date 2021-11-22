source('0_function.R')

GeoChip_TAX_reformat <- function(TAX){
  TAX=as.matrix(TAX)
  colnames(TAX)=str_to_title(colnames(TAX))
  TAX=TAX[,c("Gene","Gene_category","Subcategory1","Subcategory2","Organism","Lineage")]
  TAX[,"Gene"]=str_to_lower(TAX[,"Gene"])
  TAX[,"Gene_category"]=str_to_lower(TAX[,"Gene_category"])
  ### Geochip Domain
  Domain =  apply(TAX[,'Lineage'],1,function(x) str_split(x,';')[[1]][1]) %>% sapply(.,function(x)str_split(x,':')[[1]][2])
  Domain[is.na(Domain)]='unknown'
  Domain[Domain=='']='unknown'
  ### Geochip phylum
  Phylum =  apply(TAX[,'Lineage'],1,function(x) str_split(x,';')[[1]][2]) %>% sapply(.,function(x)str_split(x,':')[[1]][2])
  Phylum[is.na(Phylum)]='unknown'
  Phylum[Phylum=='']='unknown'
  ###
  TAX=cbind(TAX,Domain) %>% cbind(.,Phylum)
  data.frame(TAX)
}

TAX = data.frame(TAX)
summary(TAX$Gene_category)
summary(TAX$Subcategory1)

load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
TAX = tax_table(project)
TAX = TAX_reformat(TAX)
NA_full_index=paste(rownames(TAX),TAX$Gene,TAX$Domain,TAX$Phylum)
NA_genecate_index=paste(rownames(TAX),TAX$Gene,TAX$Gene_category)
NA_genename_index=paste(rownames(TAX),TAX$Gene)
NA_gene_ID = rownames(TAX)

load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
TAX = tax_table(project)
TAX = TAX_reformat(TAX)
SNNR_full_index=paste(rownames(TAX),TAX$Gene,TAX$Domain,TAX$Phylum)
SNNR_genecate_index=paste(rownames(TAX),TAX$Gene,TAX$Gene_category)
SNNR_genename_index=paste(rownames(TAX),TAX$Gene)
SNNR_gene_ID = rownames(TAX)


load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
TAX = tax_table(project)
TAX = TAX_reformat(TAX)
HKV_full_index=paste(rownames(TAX),TAX$Gene,TAX$Domain,TAX$Phylum)
HKV_genecate_index=paste(rownames(TAX),TAX$Gene,TAX$Gene_category)
HKV_genename_index=paste(rownames(TAX),TAX$Gene)
HKV_gene_ID = rownames(TAX)

### change mathced pair of different level into gene ID number
### find gene Name unmatch
gene_ID_match = intersect(SNNR_gene_ID,HKV_gene_ID)
gene_name_match = intersect(SNNR_genename_index,HKV_genename_index) %>% sapply(.,function(x)str_split(x,' ')[[1]][1])
gene_ID_name_unmatch = intersect(SNNR_gene_ID,HKV_gene_ID)[!gene_ID_match%in%gene_name_match]
# what they are
HKV_cate =HKV_genecate_index[HKV_gene_ID %in% gene_ID_name_unmatch] %>% sort(.)
HKV_name =HKV_genename_index[HKV_gene_ID %in% gene_ID_name_unmatch] %>% sort(.)
SNNR_cate=SNNR_genecate_index[SNNR_gene_ID %in% gene_ID_name_unmatch] %>% sort(.)
SNNR_name=SNNR_genename_index[SNNR_gene_ID %in% gene_ID_name_unmatch] %>% sort(.)
Check_result = data.frame(ID=sort(gene_ID_name_unmatch),HKV_name,HKV_cate,SNNR_name,SNNR_cate)

### find gene Cate unmatch
gene_ID_match = intersect(SNNR_gene_ID,HKV_gene_ID)
gene_cate_match = intersect(SNNR_genecate_index,HKV_genecate_index) %>% sapply(.,function(x)str_split(x,' ')[[1]][1])
gene_ID_cate_unmatch = intersect(SNNR_gene_ID,HKV_gene_ID)[!gene_ID_match%in%gene_cate_match]
# what they are
HKV_cate =HKV_genecate_index[HKV_gene_ID %in% gene_ID_cate_unmatch] %>% sort(.)
HKV_name =HKV_genename_index[HKV_gene_ID %in% gene_ID_cate_unmatch] %>% sort(.)
SNNR_cate=SNNR_genecate_index[SNNR_gene_ID %in% gene_ID_cate_unmatch] %>% sort(.)
SNNR_name=SNNR_genename_index[SNNR_gene_ID %in% gene_ID_cate_unmatch] %>% sort(.)
Check_result = data.frame(ID=sort(gene_ID_cate_unmatch),HKV_name,HKV_cate,SNNR_name,SNNR_cate)

### find gene Domain phlyum unmatch
gene_ID_match = intersect(SNNR_gene_ID,HKV_gene_ID)
gene_taxa_match = intersect(SNNR_full_index,HKV_full_index) %>% sapply(.,function(x)str_split(x,' ')[[1]][1])
gene_ID_taxa_unmatch = intersect(SNNR_gene_ID,HKV_gene_ID)[!gene_ID_match%in%gene_taxa_match]
# what they are
HKV_cate =HKV_genecate_index[HKV_gene_ID %in% gene_ID_taxa_unmatch] %>% sort(.)
HKV_taxa =HKV_full_index[HKV_gene_ID %in% gene_ID_taxa_unmatch] %>% sort(.)
SNNR_cate=SNNR_genecate_index[SNNR_gene_ID %in% gene_ID_taxa_unmatch] %>% sort(.)
SNNR_taxa=SNNR_full_index[SNNR_gene_ID %in% gene_ID_taxa_unmatch] %>% sort(.)
Check_result = data.frame(ID=sort(gene_ID_taxa_unmatch),HKV_taxa,HKV_cate,SNNR_name,SNNR_taxa)

### find gene Cate unmatch in GeoChip5
gene_ID_match = intersect(NA_gene_ID,HKV_gene_ID)
gene_cate_match = intersect(NA_genecate_index,HKV_genecate_index) %>% sapply(.,function(x)str_split(x,' ')[[1]][1])
gene_ID_cate_unmatch = intersect(NA_gene_ID,HKV_gene_ID)[!gene_ID_match%in%gene_cate_match]
# what they are
HKV_cate =HKV_genecate_index[HKV_gene_ID %in% gene_ID_cate_unmatch] %>% sort(.)
HKV_name =HKV_genename_index[HKV_gene_ID %in% gene_ID_cate_unmatch] %>% sort(.)
NA_cate=NA_genecate_index[NA_gene_ID %in% gene_ID_cate_unmatch] %>% sort(.)
NA_name=NA_genename_index[NA_gene_ID %in% gene_ID_cate_unmatch] %>% sort(.)
Check_result = data.frame(ID=sort(gene_ID_cate_unmatch),HKV_name,HKV_cate,NA_name,NA_cate)
