
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
library(useful)
library("cluster")
library("clustree")
# ----------------------------------------------------------------------------------


# ================================ Color Pelet =====================================
color_cond <- c( "magenta4", "#007A87",brewer.pal(6,"Dark2")[-1],"#FF5A5F","black",
    "#FFB400", "#007A87", "#8CE071", "#7B0051",
                 "#00D1C1")
color_clust <- c("#FF5A5F", "#FFB400", "#007A87", "#8CE071", "#7B0051",
                 "#00D1C1", "#FFAA91", "#B4A76C", "#9CA299", "#565A5C", "#00A04B", "#E54C20",brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust,Oligo_Pop=color_clust)
# ----------------------------------------------------------------------------------


# =============================== Load Functions ====================================
script_libs <- list.files("/home/users/dkyriakis/PhD/Projects/Yahaya/Script_Library",full.names = TRUE)
lapply(script_libs,source)
# ----------------------------------------------------------------------------------

# ============================= PATH OF FILES =======================================
# Produced by Kamil (GBMST)
path = "/work/projects/esprit/Dropseq_projects/GBMST/mouse/counts"
filenames_old <-list.files(path,full.names = TRUE)
filenames_DGE_old <- filenames_old[grep("Summary|readcounts",filenames_old,invert=TRUE)][c(1,3,5)]
cond_names_old <- c("P13S","T16S","T192S")

# Produced by Kamil New (LIH2)
path_new = "/work/projects/esprit/Dropseq_projects/LIH2/new/mouse/counts"
filenames_new <-list.files(path_new,full.names = TRUE)
filenames_DGE_new <- filenames_new[grep("Summary",filenames_new,invert=TRUE)][c(1,2)]
cond_names_new <- c("T470S","T101S")

filenames <- c(filenames_DGE_old,filenames_DGE_new)
condition_names <- c(cond_names_old,cond_names_new)


# New data 15/07/2020 from Kamil (LIHII)
path_new<-"/work/projects/esprit/Dropseq_projects/LIHII/mouse/counts"
filenames_new <-list.files(path_new,full.names = TRUE)
filenames_DGE_new <- filenames_new[grep("Summary|readcounts",filenames_new,invert=TRUE)]
cond_names_new <- c("T347S","T233S","P3TMZ","P3Con")

filenames <- c(filenames,filenames_DGE_new)
condition_names <- c(condition_names,cond_names_new)


old_suresh <- list.files("/home/users/dkyriakis/PhD/Projects/Yahaya/DATA_old",full.names=T)
filenames <- c(filenames,old_suresh)
condition_names <- c(condition_names,c("Control","p13_old","p3_old","p8"))
# -----------------------------------------------------------------------------------


setwd("/mnt/irisgpfs/projects/esprit/scAnalysis/Yahaya/")


mouse_data <- readRDS("5.Seperation/Mouse.rds")



dir.create("6.Mouse_Analysis")
setwd("6.Mouse_Analysis")
DefaultAssay(mouse_data) <-"RNA"
mouse_data2 <-liger_projection(mouse_data,split.by="condition",k=20,lamda=5,n_neighbors=10)
mouse_data2$Cluster <- mouse_data2$Liger_Clusters



mouse_data2 <- FindNeighbors(mouse_data2, dims = 1:20,reduction ="liger_inmf")
mouse_data2 <- FindClusters(mouse_data2, resolution = seq(0,1.1,0.1))
#Liger_withbatch <- FindClusters(Liger_withbatch, resolution = 0.3)
clustree <-clustree(mouse_data2,prefix="RNA_snn_res.")
mouse_data2$Cluster <- as.factor(as.numeric(mouse_data2$RNA_snn_res.0.3))
pdf("Clustree.pdf")
plot(clustree)
dev.off()


pdf("Markers.pdf")
FeaturePlot(mouse_data2, c("PLP1","P2RY12","GJA1","OLIG1","GAD2","IGF2","CCDC153"),
    order=TRUE,cols = c("lightgrey","#FDBB84","#EF6548","#D7301F","#B30000","#7F0000"),reduction="liger_umap")
dev.off()

pdf("Cells.pdf")
DimPlot(mouse_data2,group.by = "condition",cols=color_cond,reduction="liger_umap")
DimPlot(mouse_data2,group.by = "Cluster",cols=color_clust,reduction="liger_umap")
DimPlot(mouse_data2,group.by = "Phase",reduction="liger_umap")
dev.off()


dat <- data.frame(table(mouse_data2$Cluster,mouse_data2$Condition))
names(dat) <- c("Cluster","Condition","Count")

pdf("Barplot.pdf",width=12)
ggplot(data=dat, aes(x=Cluster, y=Count, fill=Condition)) + geom_bar(stat="identity")+theme_cowplot()+ scale_fill_manual(values=color_cond)
dev.off()

setwd("../")






dir.create("7.Mouse")
setwd("7.Mouse")
Seurat <- mouse_data2
Idents(Seurat)<-Seurat$Cluster
DefaultAssay(Seurat) <- "RNA"
Seurat <-NormalizeData(Seurat)
Seurat <-ScaleData(Seurat)
markers <- FindAllMarkers(Seurat,test.use="MAST",only.pos =T,latent.vars="nCount_RNA",logfc.threshold = 0.1)

markers$gene <- rownames(markers)

pdf("DoHeatmap.pdf",width=12,height=10)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = abs(avg_logFC))
DoHeatmap(object = Seurat,features = unique(top10$gene),raster = F)+ 
    theme(text = element_text(size = 6))
top3 <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = abs(avg_logFC))
DotPlot(Seurat, features = unique(top3$gene),  dot.scale = 8) + RotatedAxis()+ 
    theme(text = element_text(size = 9))
dev.off()
write.table(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.1,],"sig_avatar.txt")



#saveRDS(Seurat,"/home/users/dkyriakis/PhD/Projects/Yahaya/6.Mouse_Analysis/Mouse.rds")


#Seurat <- readRDS("/home/users/dkyriakis/PhD/Projects/Yahaya/6.Mouse_Analysis/Mouse.rds")


Seurat$Cell_Type <- as.vector(Seurat$Cluster)
Seurat$Cell_Type[Seurat$Cluster %in% c(1)] <- "Microglia"
Seurat$Cell_Type[Seurat$Cluster %in% c(2)] <- "Microglia"
Seurat$Cell_Type[Seurat$Cluster %in% c(3,5,8)] <- "Astrocytes"
Seurat$Cell_Type[Seurat$Cluster %in% c(16)] <- "Ependymal"
Seurat$Cell_Type[Seurat$Cluster %in% c(4,9,10,13,14)] <- "Endothelial"
Seurat$Cell_Type[Seurat$Cluster %in% c(6,7)] <- "Oligo"
Seurat$Cell_Type[Seurat$Cluster %in% c(11)] <- "OPC/Astro"
Seurat$Cell_Type[Seurat$Cluster %in% c(12)] <- "Microglia"
Seurat$Cell_Type[Seurat$Cluster %in% c(15)] <- "OPC"
Seurat$Cell_Type[Seurat$Cluster %in% c(17)] <- "Neurons"

pdf("Cells.pdf")
DimPlot(mouse_data2,group.by = "Cell_Type",reduction="liger_umap")
DimPlot(mouse_data2,group.by = "condition",cols=color_cond,reduction="liger_umap")
DimPlot(mouse_data2,group.by = "Cluster",cols=color_clust,reduction="liger_umap")
DimPlot(mouse_data2,group.by = "Phase",reduction="liger_umap")
dev.off()


