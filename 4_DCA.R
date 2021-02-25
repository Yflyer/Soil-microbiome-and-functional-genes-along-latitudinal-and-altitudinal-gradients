source('0_function.R')

#################### save dir #######################
folder = paste0('4_DCA')
dir.create(folder)
# This code is to investigate the cluster of OTUs or Genes
######################################
dca_table <- function(project,prefix,note,Site='Site',Gradient){
  OTU <- otu_table(project) 
  ENV <- sample_data(project)
  dis = vegdist(t(OTU),method = 'bray')
  dca <- decorana(dis)
  dt_dca = dca$rproj %>% as.data.frame(.)
  dt_dca$Gradient=ENV[[Gradient]]
  dt_dca$Site=ENV[[Site]]
  dt_dca$Study=prefix
  dt_dca$Data=note
  dt_dca
}

######################################
method = 'bray'
###################################    OTU   ############################
load('project_Latitude (North America)_OTU')
project = get(paste0('project_',prefix,'_',note))
dt_dca=dca_table(project,prefix,note,Gradient='Lat')
assign(paste0('dt_dca',prefix,'_',note),dt_dca)
###############################################################################
load('project_Altitude (SNNR)_OTU')
project = get(paste0('project_',prefix,'_',note))
dt_dca=dca_table(project,prefix,note,Gradient='Elevation')
dt_dca$Gradient=round((dt_dca$Gradient)/100,0.1)
Site=dt_dca$Site%>%summary(.)
dt_dca=dca_table(project,prefix,note,Gradient='Elevation')
dt_dca$Gradient=c(rep(2800,Site[1]),rep(2200,Site[2]),rep(1000,Site[3]),rep(1800,Site[4]),rep(1700,Site[5]),rep(2500,Site[6]))
assign(paste0('dt_dca',prefix,'_',note),dt_dca)
###############################################################################
load('project_Altitude (HKV)_OTU')
project = get(paste0('project_',prefix,'_',note))
dt_dca=dca_table(project,prefix,note,Gradient='Elevation')
assign(paste0('dt_dca',prefix,'_',note),dt_dca)
###############################################################################
############################### Geochip ################################
load('project_Latitude (North America)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
dt_dca=dca_table(project,prefix,note,Gradient='Lat')
assign(paste0('dt_dca',prefix,'_',note),dt_dca)
###############################################################################
load('project_Altitude (SNNR)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
dt_dca=dca_table(project,prefix,note,Gradient='Elevation')
dt_dca$Gradient=c(rep(2800,Site[1]),rep(2200,Site[2]),rep(1000,Site[3]),rep(1800,Site[4]),rep(1700,Site[5]),rep(2500,Site[6]))
assign(paste0('dt_dca',prefix,'_',note),dt_dca)
###############################################################################
load('project_Altitude (HKV)_GeoChip')
project = get(paste0('project_',prefix,'_',note))
dt_dca=dca_table(project,prefix,note,Gradient='Elevation')
assign(paste0('dt_dca',prefix,'_',note),dt_dca)
###############################################################################
dt_dca = Reduce(rbind,list(`dt_dcaAltitude (HKV)_GeoChip`,`dt_dcaAltitude (HKV)_OTU`,`dt_dcaAltitude (SNNR)_GeoChip`,`dt_dcaAltitude (SNNR)_OTU`,`dt_dcaLatitude (North America)_GeoChip`,`dt_dcaLatitude (North America)_OTU`))
dt_dca$Study = factor(dt_dca$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
dt_dca$Gradient = as.factor(dt_dca$Gradient)
###############################################################################

#### plot
PTs <- theme(panel.background=element_blank(),
             text = element_text(size=12),
             axis.title=element_blank(),
             legend.position="top",
             panel.border = element_rect(colour = "grey", fill=NA, size=1))
cycle = stat_ellipse(geom = "polygon",aes(fill= Gradient),level = 0.9,alpha=0.2)
round = stat_ellipse(aes(color= Gradient),level = 0.9)

plot = dt_dca %>% group_by(Study,Data) %>% 
  do(gg ={ ggplot(.,aes(DCA1, DCA2)) + 
      geom_point(shape = 16, aes( color = Gradient),size=2.5,alpha=0.5)+
      PTs+
      round+cycle+facet_grid(~Study,scales = 'free_x')
    ggsave(filename = paste0(folder,'/4_DCA_',unique(.$Study),'_',unique(.$Data),'_',method,'.jpg'),dpi = 900,width=6,height=5.5)
  }) 

