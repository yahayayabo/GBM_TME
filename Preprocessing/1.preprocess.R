# ---
# title: "GBM_TME"
# author: "Dimitrios Kyriakis"
# date: "26/6/2020"
# ---

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


# ============================ Setting Params =====================================
feat.mad <- 4
count.mad <- 4
mito.mad <- 4
n_cores <- 3
if (get_os()=="windows"){
    bpparam <- BiocParallel::SnowParam(workers = n_cores, type = "SOCK")
}else{
    bpparam <- BiocParallel::MulticoreParam(workers = n_cores)
}
# ---------------------------------------------------------------------------------

# ================================== QC ===========================================
dir.create("1.QC")
setwd("1.QC")
tm_list <- create_seurat(filenames = filenames ,conditions = condition_names,
                         elbow = FALSE,data_10x = FALSE,n_cores=1)
setwd("../")
# ---------------------------------------------------------------------------------

# ================================= Filtering ======================================
dir.create("2.Filtering")
setwd("2.Filtering")
tm_list_fl <- filtering_diy(tm_list,feat.mad = 3,count.mad = 3,mito.mad = 3,n_cores = 1)
setwd("../")
# ----------------------------------------------------------------------------------
save.image("Filtered.RData")

# ================================= Projection ======================================
dir.create("3.Projection_nobatch")
setwd("3.Projection_nobatch")
test <- nobatch_projection(tm_list_fl)
pdf("Projections.pdf")
FeaturePlot(test, c("PLP1","P2RY12","GJA1","OLIG1","GAD2","IGF2","CCDC153"),order=TRUE,cols = c("lightgrey","#FDBB84","#EF6548","#D7301F","#B30000","#7F0000"))
dev.off()
pdf("Cells.pdf",width=12)
DimPlot(test,group.by = "condition",cols=color_cond)
dev.off()
setwd("../")



dir.create("3.Projection")
setwd("3.Projection")

require(liger)
require(SeuratWrappers)
test <-liger_projection(test,split.by="condition",k=20,lamda=5,n_neighbors=10)
# test <- SCT_projection(tm_list_fl)
pdf("Markers.pdf")
FeaturePlot(test, c("PLP1","P2RY12","GJA1","OLIG1","GAD2","IGF2","CCDC153"),
    order=TRUE,cols = c("lightgrey","#FDBB84","#EF6548","#D7301F","#B30000","#7F0000"),reduction="liger_umap")
dev.off()

pdf("Cells.pdf",width=12)
DimPlot(test,group.by = "condition",cols=color_cond,reduction="liger_umap")
dev.off()
setwd("../")
# ----------------------------------------------------------------------------------

save.image("Liger.RData")

# ================================= Clustering ======================================

dir.create("4.Clustering")
setwd("4.Clustering")
source('/home/users/dkyriakis/PhD/Projects/Yahaya/Script_Library/Clustering.R')
library("cluster")
library("clustree")
# =================== Optimal Clusters ======================================
# optimal_output<- optimal_clusters(test,k.max=10,save=TRUE,resolution=FALSE)
# test <- optimal_output$object
# opt_num <- optimal_output$opt_num
# sil_scor <- optimal_output$sil_scor
# test <- FindClusters(test, resolution = seq(0,1.1,0.1))


# clustree <-clustree(test,prefix="integrated_snn_res.")
# clustree <-clustree(test,prefix="SCT_snn_res.")
test$Cluster <- test$Liger_Clusters
pdf("Clustering.pdf")
# plot(clustree)
DimPlot(test,group.by = c("Cluster"),reduction="liger_umap",cols=color_clust)
DimPlot(test,group.by = c("Phase"),reduction="liger_umap")
dev.off()
# ---------------------------------------------------------------------------
graphics.off()
setwd("../")

saveRDS(test,"Seurat_obj.rds")
save.image("All_data.RData")