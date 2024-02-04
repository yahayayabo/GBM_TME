
####  Figure 3  #####

library(Seurat)
library(dittoSeq)

#Fig 3C
DimPlot(Merged_myeloid_data, group.by = "Cell_Type", split.by = "dataset", cols = c("BAM" = "red" , 'MG' = '#4daf4a', 'Mo/MG' = 'blue'))

# feature plots of some key marker genes

FeaturePlot(Merged_myeloid_data, features = c("Itgam","P2ry12", "Ly6c2", "Mrc1"))


#Fig 3D

#first per datasets
dittoBarPlot(Merged_myeloid_data,"CellType", group.by = "dataset", 
             var.labels.reorder = c(1,3,2),
             x.reorder = c(2,1,4,3), 
             color.panel = c("BAM" = "red" ,'Mo/MG' = 'blue', 'MG' = '#4daf4a')) +
  theme(text = element_text(size = 16))
