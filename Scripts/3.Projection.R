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




cat("SCT")
tm_list_fl <- readRDS("/home/users/dkyriakis/PhD/Projects/epi-scRNA/2.Preprocess/Fl_Object_list.rds")

dir.create("3.Projection")
setwd("3.Projection")
SCT_d <- SCT_projection(tm_list_fl)
pdf("Cells.pdf")
DimPlot(SCT_d,group.by = "condition",cols=color_cond)
DimPlot(SCT_d,group.by = "Phase")
dev.off()

saveRDS(SCT_d,"/home/users/dkyriakis/PhD/Projects/epi-scRNA/3.Projection/SCT_Merged.rds")
setwd("../")