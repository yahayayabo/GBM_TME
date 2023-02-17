
#4a

library(Seurat)
library(dittoSeq)
library(ggplot2)

# Dotplot of signature scores in myeloid subsets

#load myeloid_object.rds

myeloid_obj <- readRDS("myeloid_object.rds")

#load signature lists
Phago_gene_list <- read.csv("~/Phago_gene_list.csv", row.names=1)
APC_gene_list <- read.csv("~/APC_gene_list.csv", row.names=1)
Sensome_gene_list <- read.csv("~/Sensome_gene_list.csv", row.names=1) 
Migration_gene_list <- read.csv("~/Migration_gene_list.csv", row.names=1)


#score cell 
myeloid_obj <- AddModuleScore(myeloid_obj, list(APC_gene_list), name = "APC_score", search = T)
myeloid_obj <- AddModuleScore(myeloid_obj, list(Phago_gene_list), name = "Phagocytosis_score", search = T)
myeloid_obj <- AddModuleScore(myeloid_obj, list(Sensome_gene_list), name = "Sensome_score", search = T)
myeloid_obj <- AddModuleScore(myeloid_obj, list(Migration_gene_list), name = "Migration_score", search = T)

#plot dotplot
DotPlot(myeloid_obj, features = c("APCs_score", "Phagocytosis_score", "Sensome_score", "Migration_score"), 
        group.by = "Clusters",  scale = TRUE, dot.scale = 8, cols = c("blue", "red")) +   
  RotatedAxis()+coord_flip() +
  theme(legend.key.height = unit(0.4, 'cm'), legend.text = element_text(family = 'Helvetica', size = 10), 
        legend.title = element_text(family = 'Helvetica', size = 10))


### Calculate the significance of the signatures between myeloid clusters
library(escape)

### Significance Testing using Kruskal wallis test

ES <- data.frame(myeloid_obj[[]], Idents(myeloid_obj))

output_KW <- getSignificance(ES, 
                             group = "Clusters", 
                             fit = "KW")




