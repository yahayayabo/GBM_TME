
### Figure 2  ####   
Merged_myeloid_data$Final_Annot <- as.numeric(Merged_myeloid_data$Cluster) -1
Merged_myeloid_data$Final_Annot[Merged_myeloid_data$Cluster%in% c(2,10,14,15)] <- "Mo-TAMs"
Merged_myeloid_data$Final_Annot[Merged_myeloid_data$Cluster%in% c(0,1,3,4,5,6,7,8,11,12)] <- "Mg-TAMs"
Merged_myeloid_data$Final_Annot[Merged_myeloid_data$Cluster%in% c(9,13,16)] <- "BAMs"



# feature plots of some key marker genes
DimPlot(Merged_myeloid_data,group.by="Final_Annot",label=T)+NoLegend()
features <- c("Itgam","P2ry12","Ly6c2","Mrc1")
FeaturePlot(Seurat,features= features,order=T)

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

