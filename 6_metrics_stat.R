source('0_function.R')

#################### save dir #######################
folder = paste0('6_metrics_effect')
dir.create(folder)
####### 
# hkv
dt_hkv = read.csv(file = '6_metrics_effect/Three-metrics_comparision_data_Altitude (HKV).csv',row.names = 1)

stat_summary = dt_hkv %>% group_by(Method,Comparison,Metric) %>% summarise(mean.Sim = label_percent(accuracy=0.001)(mean(Dist)),sd.Sim = label_percent(accuracy=0.001)(sd(Dist)))

write.csv(stat_summary,paste0(folder,'/stat_summary_hkv.csv'),row.names = T)
# snnr
dt_snnr = read.csv(file = '6_metrics_effect/Three-metrics_comparision_data_Altitude (SNNR).csv',row.names = 1)

stat_summary = dt_snnr %>% group_by(Method,Comparison,Metric) %>% summarise(mean.Sim = label_percent(accuracy=0.001)(mean(Dist)),sd.Sim = label_percent(accuracy=0.001)(sd(Dist)))

write.csv(stat_summary,paste0(folder,'/stat_summary_snnr.csv'),row.names = T)
# NA
dt_na = read.csv(file = '6_metrics_effect/Three-metrics_comparision_data_Latitude (North America).csv',row.names = 1)

stat_summary = dt_na %>% group_by(Method,Comparison,Metric) %>% summarise(mean.Sim = label_percent(accuracy=0.001)(mean(Dist)),sd.Sim = label_percent(accuracy=0.001)(sd(Dist)))

write.csv(stat_summary,paste0(folder,'/stat_summary_NA.csv'),row.names = T)

######## process merge dt ####
dt = Reduce(rbind,list(dt_hkv,dt_snnr,dt_na))
dt$Pair=paste(dt$row.Site,dt$col.Site)
dt$Log.Dist=log10(dt$Dist+1)
dt$Study = factor(dt$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
dt$Metric = factor(dt$Metric,levels = c('Composition','Abundance','Diversity'))
write.csv(dt,paste0(folder,'/stat_data.csv'),row.names = T)

### stat exam
stat_test = dt %>% group_by(Method,Comparison,Study) %>% do(
  ab_vs_dv_t=get_ttest_t(.,m1='Abundance',m2='Diversity'),
  ab_vs_dv_p=get_ttest_p(.,m1='Abundance',m2='Diversity'),
  dv_vs_cp_t=get_ttest_t(.,m1='Diversity',m2='Composition'),
  dv_vs_cp_p=get_ttest_p(.,m1='Diversity',m2='Composition'),
  ab_vs_cp_t=get_ttest_t(.,m1='Abundance',m2='Composition'),
  ab_vs_cp_p=get_ttest_p(.,m1='Abundance',m2='Composition'),
  ab_vs_dv_cs=get_kruskal_cs(.,m1='Abundance',m2='Diversity'),
  ab_vs_dv_kp=get_kruskal_p(.,m1='Abundance',m2='Diversity'),
  dv_vs_cp_cs=get_kruskal_cs(.,m1='Diversity',m2='Composition'),
  dv_vs_cp_kp=get_kruskal_p(.,m1='Diversity',m2='Composition'),
  ab_vs_cp_cs=get_kruskal_cs(.,m1='Abundance',m2='Composition'),
  ab_vs_cp_kp=get_kruskal_p(.,m1='Abundance',m2='Composition'),
  )
stat_test = apply(stat_test,2,as.character)
write.csv(as.data.frame(stat_test),paste0(folder,'/stat_test.csv'),row.names = F)

### specifically check:
### jaccard and bray have the same chi-square, check it
test1 = filter(dt,Study=='Altitude (HKV)'&Comparison=='Within-group'&Method=='Jaccard')
test1$Dist
kruskal.test(Dist ~ Metric, data = test1)
test2 = filter(dt,Study=='Altitude (HKV)'&Comparison=='Within-group'&Method=='Bray-Curits')
test2$Dist
kruskal.test(Dist ~ Metric, data = test2)
             
### plot box
PTs = theme(axis.text.x=element_blank(),
            axis.title=element_blank(),
            legend.position="bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            panel.background=element_blank(),
            panel.border = element_rect(colour = "grey", fill=NA, size=1))

plot_dt = dt %>% group_by(Method,Comparison,Study,Pair,Metric) %>% summarise(Dist=mean(Dist),Log.Dist=mean(Log.Dist)) %>% group_by(Comparison) %>% do (
  plot_point(.)
)

### function test ###
# get test data
test = filter(dt,Method=='Bray-Curits'&Study=='Altitude (HKV)') %>% group_by(Method,Comparison,Study,Pair,Metric) %>% summarise(Dist=mean(Dist),Log.Dist=mean(Log.Dist))
# test plot
data=test
plot_point <- function(data) {
  st_mst <- ggplot(aes(x=Metric,y=Dist),data = data)+geom_boxplot(aes(color=Metric),lwd=0.25)#+geom_jitter(aes(color=Metric),shape=16,alpha=0.3,position=position_jitter(0.2),size=1)
  st_mst+PTs+facet_grid(Method~Study,scales = "free") # , space = "free"
  #ggsave(filename = paste0(folder,'/simbox_',unique(data$Comparison),'.pdf'),dpi=600,width=17,height=20,units='cm')
}
### add signif
library(ggsignif)
st_mst <- ggplot(aes(x=Metric,y=Dist),data = data)+geom_boxplot(aes(color=Metric),lwd=0.25)#+geom_jitter(aes(color=Metric),shape=16,alpha=0.3,position=position_jitter(0.2),size=1)
st_mst+PTs+facet_grid(Method~Study,scales = "free")+geom_signif(
  comparisons = list(c("Composition", "Abundance",'Diversity')),
  map_signif_level = TRUE
)
plot_point(test)
# get test data
test = filter(dt,Method=='Sorensen'&Comparison=='Within-group')
# test stat 
subdt=test
get_ttest_p <- function(subdt,m1='Abundance',m2='Composition'){
  x = subdt$Dist[subdt$Metric==m1]
  y = subdt$Dist[subdt$Metric==m2]
  tmodel = t.test(x,y)
  tmodel$p.value %>% specify_decimal(.,3)
}
get_ttest_t <- function(subdt,m1='Abundance',m2='Composition'){
  x = subdt$Dist[subdt$Metric==m1]
  y = subdt$Dist[subdt$Metric==m2]
  tmodel = t.test(x,y)
  tmodel$statistic %>% specify_decimal(.,3)
}
get_kruskal_p <- function(subdt,m1='Abundance',m2='Composition'){
  test = filter(subdt,Metric%in%c(m1,m2))
  kmodel = kruskal.test(Dist ~ Metric, data = test) #%>% summary()
  kmodel$p.value %>% specify_decimal(.,3)
}
get_kruskal_cs <- function(subdt,m1='Abundance',m2='Composition'){
  test = filter(subdt,Metric%in%c(m1,m2))
  kmodel = kruskal.test(Dist ~ Metric, data = test) #%>% summary()
  kmodel$statistic %>% specify_decimal(.,3)
}
# 
get_ttest_p(test)
get_ttest_t(test)
get_kruskal_p(test)
get_kruskal_cs(test)

