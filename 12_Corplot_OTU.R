source('0_function.R')
library(ggthemes)
library(pheatmap)
library(viridis)
# ,"Plant_richness"
ENV_var = c("NH4.N","NO3.N","TN","TC","CN_ratio","pH","Rainfall","Temperature")
set.seed(0430)

load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))

load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))

load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))

################################################
ENV <- sample_data(project) %>% data.frame(.) 
# define OA levels
Site = ENV$Site %>% as.factor(.)
####################################################
OTU <- otu_table(project)
#OA_info = group_count(table = OTU, group = Site) 
#OTU = OTU[rowSums(OA_info)>=5,] 
seed = sample(nrow(OTU),2000)
OTU = OTU[seed,]

#########################################################
ENV = ENV[,ENV_var]

Cor.OTU=c()
Cor.dt=data.frame(matrix(nrow=nrow(OTU),ncol=ncol(ENV),dimnames = list(rownames(OTU),colnames(ENV))))
for (j in 1:ncol(ENV)) {
  for (i in 1:nrow(OTU)) {
    index = OTU[i,]!=0
    Env_factor = ENV[index,j] %>% as.numeric(.)
    if (sum(index)>=6&all(!is.na(Env_factor))) {
      Abundance=OTU[i,index] %>% as.numeric(.)
      Cor = cor.test(Abundance,Env_factor,method = 'spearman',exact = F)
      Cof = Cor$estimate
    } else {Cof=NA}
    Cor.OTU[i]=Cof
  }
  if (sum(!is.na(Cor.OTU))>0.5*nrow(OTU)) {
    Cor.dt[[ENV_var[j]]] = Cor.OTU
  }
}
Cor.dt = Cor.dt[!rowSums(!is.na(Cor.dt))==0,!colSums(!is.na(Cor.dt))==0]
Cor.dt[is.na(Cor.dt)]=0
#Cor.dt = Cor.dt[,c("NH4.N","NO3.N","TN","TC","CN_ratio","Plant_richness","pH","Rainfall","Temperature")] 

#colnames(Cor.dt) =  c('Ammonium','Nitrate','Total nitrogen', 'Total carbon', 'C/N ratio', 'Plant richness', 'pH', 'Rainfall', 'Temperature')

######
cluster.k=4
#cluster_tree <- hclust(dist(Cor.dt), method = "complete")
########
#cluster_info = cutree(tree = as.dendrogram(cluster_tree),k=4) 
#cluster_info = cluster_info[rownames(Cor.dt)] 
#cluster_info = cluster_info[order(rownames(Cor.dt))] %>% paste('Cluster',.)
###############
taxon_info = tax_table(project)[rownames(Cor.dt),'phylum'] %>% .[order(rownames(Cor.dt))]
taxon_info[!taxon_info%in%c('Acidobacteria','Actinobacteria','Bacteroidetes','Chloroflexi','Protobacteria','Planctomycetes','Gemmatimonadetes','Verrucomicrobia')]='Others'
#################### cluster_info,
OA_info = group_count(table = OTU, group = Site)%>% .[rownames(Cor.dt),] %>% apply(., 1, function(x) sum(x>=1))
summary(OA_info)

OA_info[OA_info==nlevels(Site)] = "Ubiquitous"
OA_info[OA_info==1] = "Endemic"
OA_info[!OA_info%in%c("Endemic","Ubiquitous")] = "Occasional"
OA_info=OA_info[order(rownames(Cor.dt))]
OA_info%>% as.factor(.) %>% summary(.)
####################################
`Environmental Gradients`=c("Latitude (North America)",'Altitude (SNNR)','Altitude (HKV)')
Gradient_info=c()
Gradient_info[1:nrow(Cor.dt)]=`Environmental Gradients`[`Environmental Gradients`==prefix]
########################################

anno_row_dt = data.frame(taxon_info,OA_info,Gradient_info)
anno_col_dt = data.frame('Environmental_factor'=colnames(Cor.dt))
rownames(anno_col_dt)=colnames(Cor.dt)

#cividis(10)
#viridis(10)
#magma(10)
#library(ggsci)


pheatmap(Cor.dt,
         color =magma(10)[2:10],
         #breaks = c(-1,-0.6,-0.2,0.2,0.6,1),
         annotation_colors = anno_colour,
         annotation_row = anno_row_dt,
         annotation_col = anno_col_dt,
         annotation_names_row = F,
         annotation_names_col = F,
         cutree_rows =cluster.k,
         cluster_cols = F,
         #cluster_rows = F,
         show_rownames = F,
         show_colnames = F,
         width = 5,
         legend = T,
         filename = paste0('Corplot_',prefix,'_',note,'.pdf'))

##############################
row_col = c()
row_col['NH4.N']='#ef3b2c'
row_col['NO3.N']='#fb6a4a'
row_col['TN']='#fc9272'
row_col['TC']='#fcbba1'
row_col['CN_ratio']='#fee0d2'
row_col['Plant_richness']='#fff5f0'
row_col['Temperature']='#bdbdbd'
row_col['Rainfall']='#969696'
row_col['pH']='#737373'

Taxa_color = c('orangered','#43a2ca','#018571','#3C5488','#f4a582','slategray','#80cdc1','#ca0020','#b2abd2','#5e3c99')
ggthemes_data$wsj$palettes$colors6
OA_color = c(Endemic= '#c72e29',Occasional ='#016392',Ubiquitous='#be9c2e')
Gradient_color = c(`Latitude (North America)`='#08519c',`Altitude (SNNR)`='#3182bd',`Altitude (HKV)`='#6baed6')

anno_colour = list(Environmental_factor=row_col,OA_info=OA_color,Gradient_info=Gradient_color)

