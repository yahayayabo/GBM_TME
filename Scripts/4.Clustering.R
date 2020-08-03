# ---
# title: "Epi-Analysis"
# author: "Dimitrios Kyriakis"
# date: "09/07/2020"
# ---
options(future.globals.maxSize= 20122317824)
# ============================== Libraries ========================================
library(HDF5Array)
library(BiocParallel)
library(tidyverse)
library(tictoc)
library(Seurat)
library(RColorBrewer)
library(tictoc)
library(crayon)
library(NMF)
library(cowplot)
library(ggpubr)
library(org.Hs.eg.db)
library("cluster")
library("clustree")

# ----------------------------------------------------------------------------------


# ================================ Color Pelet =====================================
color_cond <- c( "magenta4", "#007A87",brewer.pal(6,"Dark2")[-1],"#FF5A5F","black")
color_clust <- c("#FF5A5F", "#FFB400", "#007A87", "#8CE071", "#7B0051",
                 "#00D1C1", "#FFAA91", "#B4A76C", "#9CA299", "#565A5C", "#00A04B", "#E54C20",brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust,Oligo_Pop=color_clust)
# ----------------------------------------------------------------------------------


# =============================== Load Functions ====================================
script_libs <- list.files("/home/users/dkyriakis/PhD/Projects/epi-scRNA/Script_Library",full.names = TRUE)
lapply(script_libs,source)
# ----------------------------------------------------------------------------------

# ============================= PATH OF FILES =======================================
path <- "/mnt/irisgpfs/users/dkyriakis/PhD/Projects/epi-scRNA"


filenames<- c("/work/projects/esprit/epi-scRNA/Control/outs/filtered_feature_bc_matrix/",
"/work/projects/esprit/epi-scRNA/Epilepsy/outs/filtered_feature_bc_matrix/")
condition_names <- c("control","epilepsy")
# -----------------------------------------------------------------------------------
setwd(path)



test <- readRDS("/home/users/dkyriakis/PhD/Projects/epi-scRNA/3.Projection/SCT_Merged.rds")




setwd("4.Clustering")
# =================== Optimal Clusters ======================================
test <- FindClusters(test, resolution = seq(0,1.1,0.1))
clustree <-clustree(test,prefix="integrated_snn_res.")
test$Cluster <- test$integrated_snn_res.0.4
pdf("Clustering.pdf")
plot(clustree)
DimPlot(test,group.by = c("Cluster"),reduction="umap",cols=color_clust)
DimPlot(test,group.by = c("condition"),reduction="umap")
dev.off()





# ================================================== MARKERS ===================================================
excit_neurons <- c("SYT1",
"RBFOX3",
"CUX2",
"SATB2",
"RORB",
"TLE4")

neurons <-c("GAD1",
"GAD2",
"PVALB",
"SST",
"VIP",
"SV2C")

Glia<-c("NRGN",
"THY1",
"SLC1A2",
"GFAP",
"PLP1",
"PDGFRA",
"PTPRC",
"CLDN5")




genelist<-c(excit_neurons,neurons,Glia)
DefaultAssay(test)<-"RNA"
pdf("PAPER_Markers_RNA.pdf")
for (gene in genelist){
  plot(FeaturePlot(test,gene ,order=TRUE,cols = c("lightgrey","#FDBB84","#EF6548","#D7301F","#B30000","#7F0000"),reduction="umap"))
}
dev.off()


genelist<-c(excit_neurons,neurons,Glia)
DefaultAssay(test)<-"SCT"
pdf("PAPER_Markers_SCT.pdf")
for (gene in genelist){
  plot(FeaturePlot(test,gene ,order=TRUE,cols = c("lightgrey","#FDBB84","#EF6548","#D7301F","#B30000","#7F0000"),reduction="umap"))
}
dev.off()


# ================================= ASSIGN CLUSTERS TO CELL TYPES ====================================

test$Cell_Type <- as.vector(test$Cluster)
test$Cell_Type[test$Cell_Type %in% c(11)] <- "AST-PP"
test$Cell_Type[test$Cell_Type %in% c(15)] <- "AST-FB"

test$Cell_Type[test$Cell_Type %in% c(14)] <- "OPC"
test$Cell_Type[test$Cell_Type %in% c(10)] <- "Oligo"

test$Cell_Type[test$Cell_Type %in% c(12)] <- "Endothelial"
test$Cell_Type[test$Cell_Type %in% c(2,4,7,13,16)] <- "Interneurons"
test$Cell_Type[test$Cell_Type %in% c(8)] <-"L5/6"
test$Cell_Type[test$Cell_Type %in% c(0,1,3,5,6,9,17)] <-"Neu"

pdf("Split_Cond.pdf")
plot(DimPlot(test, reduction = "umap", group.by = "Cell_Type", split.by = "condition"))
plot(DimPlot(test, reduction = "umap", group.by = "condition"))
plot(DimPlot(test, reduction = "umap", group.by = "Cluster"))
plot(DimPlot(test, reduction = "umap", group.by = "Cell_Type"))
dev.off()
# -----------------------------------------------------------------------------------------------------



DefaultAssay(test) <- "RNA"
Idents(test)<-test$Condition
nk.markers <- FindConservedMarkers(test, ident.1 = "control", grouping.var = "Cell_Type", verbose = FALSE)
write.table(nk.markers,"Conserved.tsv")
nk.markers_ord <- nk.markers[order( nk.markers$minimump_p_val),]
library(patchwork)
pdf("Conserved.pdf",width=12)
plots <- VlnPlot(test, features =rownames(nk.markers_ord)[1:3], split.by = "condition", group.by = "Cell_Type", 
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)
Idents(test)<-test$Cell_Type
plot(DotPlot(test, features = rownames(nk.markers_ord)[1:40], cols = c("blue", "red"), dot.scale = 8, split.by = "condition")+ 
    RotatedAxis())
test2 <- ScaleData(test,features=rownames(nk.markers)[1:40])
plot(DoHeatmap(test2, features = rownames(nk.markers)[1:40]) + NoLegend())
dev.off()



saveRDS(test,"/home/users/dkyriakis/PhD/Projects/epi-scRNA/4.Clustering/SCT_Merged_Cl.rds")

setwd("../")