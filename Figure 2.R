
### Figure 2  ####   




# feature plots of some key marker genes

FeaturePlot(Merged_myeloid_data, features = c("Itgam","P2ry12", "Ly6c2", "Mrc1"))


#2b

#first per datasets
dittoBarPlot(Merged_myeloid_data,"Cell_Type", group.by = "dataset", 
             var.labels.reorder = c(1,3,2),
             x.reorder = c(2,1,4,3), 
             color.panel = c("BAM" = "red" ,'Mo/MG' = 'blue', 'MG' = '#4daf4a')) +
  theme(text = element_text(size = 16))

##subset only NORLUX (PDOX and Yolanda) datasets and show stack chart
norlux_datasets <- subset(Merged_myeloid_data, subset=dataset%in%c("PDOX", 'Yolanda'))
table(norlux_datasets$dataset)
#PDOX  Bozena Yolanda   Pombo 
#5239       0    1336       0 
table(norlux_datasets$orig.ident)

dittoBarPlot(norlux_datasets,"Cell_Type", group.by = "orig.ident", 
             var.labels.reorder = c(1,3,2),
             x.reorder = c(5,8,9,13,14,11,12,7,10,6,4,3,1,2), 
             color.panel = c("BAM" = "red" ,'Mo/MG' = 'blue', 'MG' = '#4daf4a')) +
  theme(text = element_text(size = 16))

