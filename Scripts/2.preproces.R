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


setwd("/home/users/dkyriakis/PhD/Projects/Yahaya/")

tm_list <- readRDS("/home/users/dkyriakis/PhD/Projects/Yahaya/1.QC/QC_list_data.rds")
# ================================= Filtering ======================================
dir.create("2.Filtering")
setwd("2.Filtering")
tm_list_fl <- filtering_diy(tm_list,feat.mad = 3,count.mad = 3,mito.mad = 3,n_cores = 1)
saveRDS(tm_list_fl,"/home/users/dkyriakis/PhD/Projects/Yahaya/2.Preprocess/Fl_QC_list_data.rds")
setwd("../")
# ----------------------------------------------------------------------------------
