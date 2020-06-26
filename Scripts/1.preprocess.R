---
title: "GBM_TME"
author: "Dimitrios Kyriakis"
date: "26/6/2020"
---

# ============================== Libraries ========================================
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
color_cond <- c( "magenta4", "#007A87",brewer.pal(6,"Dark2")[-1],"#FF5A5F","black")
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
# Produced by Suresh
path = "/work/projects/esprit/Dropseq_projects/GBMST/mouse/counts"

# Produced by Kamil
path_new = "/work/projects/esprit/Dropseq_projects/LIH2/new/mouse/counts"
filenames_new <- grep("DGE",list.files(path,full.names = TRUE))[1:2]
condition_names <- c("T470S","T101S")
# -----------------------------------------------------------------------------------


# ============================ Setting Params =====================================
feat.mad <- 4
count.mad <- 4
mito.mad <- 4
n_cores <- 4
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
                         elbow = TRUE,data_10x = FALSE)
setwd("../")
# ---------------------------------------------------------------------------------

# ================================= Filtering ======================================
dir.create("2.Filtering")
setwd("2.Filtering")
tm_list_fl <- filtering_diy(tm_list,feat.mad = 2,count.mad = 2,mito.mad = 1.5,n_cores = 4)
setwd("../")
# ----------------------------------------------------------------------------------

