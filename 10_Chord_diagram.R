# devtools::install_github("jokergoo/circlize")
library('circlize')
library('ggsci')
library(scales)
library(dplyr)
mypal <- pal_npg("nrc", alpha = 0.7)(9)
show_col(mypal)

chord_data <- read.csv("10_RF/Total_Gene.csv",row.names = 1) %>% as.matrix(.)


##############################
row_col = c('#636363','#08519c' ,'#3182bd','#6baed6','#d95f0e','#fec44f','#fff7bc','#f03b20','#fb6a4a','#fcae91')
row.names(chord_data)
names(row_col)=row.names(chord_data)
row_col['NH4.N']='#ef3b2c'
row_col['NO3.N']='#fb6a4a'
row_col['TN']='#fc9272'
row_col['TC']='#fcbba1'
row_col['CN_ratio']='#fee0d2'
row_col['Plant_richness']='#fff5f0'
row_col['Temperature']='#bdbdbd'
row_col['Rainfall']='#969696'
row_col['pH']='#737373'
row_col['Distance']='#525252'
##############################
col_col= c('#08519c','#3182bd' ,'#6baed6')
colnames(chord_data)
names(col_col)=colnames(chord_data)
##########################

circos.par(gap.after = row_gap)
chordDiagram(chord_data, transparency = 0.5,row.col = row_col,grid.col = c(col_col,row_col),
             link.sort = T, link.decreasing = T)

#png("my_plot.pdf") 
row_gap = c(rep(1, nrow(chord_data)-1), 15, rep(6, ncol(chord_data)-1), 15)
circos.par(gap.after = row_gap)
chordDiagram(chord_data, transparency = 0.5,row.col = row_col,grid.col = c(col_col,row_col),
             link.sort = T, link.decreasing = T,
             directional = 1)
circos.clear()
#dev.off()

chord_data <- read.csv("10_RF/Total_OTU.csv",row.names = 1) %>% as.matrix(.)


##############################
row_col = c('#636363','#08519c' ,'#3182bd','#6baed6','#d95f0e','#fec44f','#fff7bc','#f03b20','#fb6a4a','#fcae91')
row.names(chord_data)
names(row_col)=row.names(chord_data)
row_col['NH4.N']='#ef3b2c'
row_col['NO3.N']='#fb6a4a'
row_col['TN']='#fc9272'
row_col['TC']='#fcbba1'
row_col['CN_ratio']='#fee0d2'
row_col['Plant_richness']='#fff5f0'
row_col['Temperature']='#bdbdbd'
row_col['Rainfall']='#969696'
row_col['pH']='#737373'
row_col['Distance']='#525252'
##############################
col_col= c('#08519c','#3182bd' ,'#6baed6')
colnames(chord_data)
names(col_col)=colnames(chord_data)
##########################

circos.par(gap.after = row_gap)
chordDiagram(chord_data, transparency = 0.5,row.col = row_col,grid.col = c(col_col,row_col),
             link.sort = T, link.decreasing = T)

#png("my_plot.pdf") 
row_gap = c(rep(1, nrow(chord_data)-1), 15, rep(6, ncol(chord_data)-1), 15)
circos.par(gap.after = row_gap)
chordDiagram(chord_data, transparency = 0.5,row.col = row_col,grid.col = c(col_col,row_col),
             link.sort = T, link.decreasing = T,
             directional = 1)
circos.clear()

