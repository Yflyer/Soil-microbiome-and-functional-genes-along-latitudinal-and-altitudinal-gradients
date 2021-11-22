source('0_function.R')
library(ggthemes)
library(pheatmap)
library(viridis)
# ,"Plant_richness"
ENV_var = c("NH4.N","NO3.N","TN","TC","CN_ratio","pH","Rainfall","Temperature")
set.seed(123)

load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))

load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))

load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))


OTU <- otu_table(project) 
################################################
ENV <- sample_data(project) %>% data.frame(.) 
Site = ENV$Site %>% as.factor(.)
#########################################################
ENV = ENV[,ENV_var]

Cor.OTU=c()
P.OTU=c()
Cor.dt=data.frame(matrix(nrow=nrow(OTU),ncol=ncol(ENV),dimnames = list(rownames(OTU),colnames(ENV))))
P.dt=data.frame(matrix(nrow=nrow(OTU),ncol=ncol(ENV),dimnames = list(rownames(OTU),colnames(ENV))))
for (j in 1:ncol(ENV)) {
  for (i in 1:nrow(OTU)) {
    index = OTU[i,]!=0
    Env_factor = ENV[index,j] %>% as.numeric(.)
    if (sum(index)>6&all(!is.na(Env_factor))) {
      Abundance=OTU[i,index] %>% as.numeric(.)
      Cor = cor.test(Abundance,Env_factor,method = 'spearman',exact = F)
      Cof = Cor$estimate
      P = Cor$p.value
    } else {Cof=P=NA}
    Cor.OTU[i]=Cof
    P.OTU[i]=P
  }
  P.OTU=p.adjust(P.OTU, method = "BH")
  #if (sum(!is.na(Cor.OTU))>0.5*nrow(OTU)) {
    Cor.dt[[ENV_var[j]]] = Cor.OTU
    P.dt[[ENV_var[j]]] = P.OTU
  #}
}
Cor.dt = Cor.dt[!rowSums(!is.na(Cor.dt))==0,!colSums(!is.na(Cor.dt))==0]
Cor.dt[is.na(Cor.dt)]=0


cluster.k=4
########

###############
taxon_info = tax_table(project)[rownames(Cor.dt),'Gene_category'] %>% .[order(rownames(Cor.dt))]
taxon_info[taxon_info=='nitrogen']='Nitrogen'
taxon_info[taxon_info=='Carbon cycling']='Carbon Cycling'
taxon_info[!taxon_info%in%c('Carbon Cycling','Nitrogen')]='Others'
#################### cluster_info,
type_var = c("Fun","Fun","Fun","Fun","Fun","Env","Env","Env")
######
type_info = group_count(table = P.dt, group = type_var)
type_info$type='Both uncorrelated'
type_info$type[type_info$Env>0]='Envonmental correlation'
type_info$type[type_info$Fun>0]='Soil-functional correlation'
type_info$type[intersect(type_info$Env>0,type_info$Fun>0)]='Both correlated'
type_info=type_info$type

OA_info = group_count(table = OTU, group = Site)
OA_info$OA_info = "Occasional"
OA_info$OA_info[OA_info==nlevels(Site)] = "Ubiquitous"
OA_info$OA_info[OA_info==1] = "Endemic"
OA_info=OA_info$OA_info

OA_info%>% as.factor(.) %>% summary(.)
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

ggthemes_data$wsj$palettes$colors6
OA_color = c(Endemic= '#c72e29',Occasional ='#016392',Ubiquitous='#be9c2e')
Gene_color = c(`Carbon Cycling`= '#9970ab',Nitrogen ='#fdb863',Others='#bababa')
Gradient_color = c(`Latitude (North America)`='#08519c',`Altitude (SNNR)`='#3182bd',`Altitude (HKV)`='#6baed6')

anno_colour = list(Environmental_factor=row_col,OA_info=OA_color,Gene_category = Gene_color,Gradient_info=Gradient_color)

