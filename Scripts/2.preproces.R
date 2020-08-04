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




tm_list_fl <- lapply(c(1:length(condition_names)), FUN=function(x) { 

        # =============================== Load Functions ====================================
        script_libs <- list.files("/home/users/dkyriakis/PhD/Projects/epi-scRNA/Script_Library",full.names = TRUE)
        lapply(script_libs,source)
        # ----------------------------------------------------------------------------------
        
        object <- list.object[[x]]
        condition <-condition_names[x]
        cell_data  <- object@meta.data

        out_indexes <- outlier_indexes(cell_data,feat.mad = feat.mad,count.mad = count.mad,mito.mad = mito.mad)
        out_indexes_n <- rep(TRUE,dim(cell_data)[1])
        out_indexes_n[out_indexes]<-FALSE
        object$Keep2 <- out_indexes_n

        df_qc <- as.data.frame(object@meta.data)
        colnames(df_qc)[1] <- "Cond"
        out_indexes <- outlier_indexes(cell_data,feat.mad = feat.mad,count.mad = count.mad,mito.mad = mito.mad)
        out_indexes_n <- rep(TRUE,dim(cell_data)[1])
        out_indexes_n[out_indexes]<-FALSE

        
        pdf(paste("2.Filtering",condition,".pdf"))
        Vln_QC2(df_qc,condition=condition,title=condition,outliers=out_indexes_n)
        
        object$Keep <- out_indexes_n
        object_aft <- subset(object,subset = Keep == TRUE)
        df_qc <- as.data.frame(object_aft@meta.data)
        colnames(df_qc)[1] <- "Cond"
        Vln_QC2(df_qc,condition=condition,title=condition)
        dev.off()
        
        object_aft$condition <- object_aft$Condition
        
        object_aft$Joint <- (object_aft$nFeatures / object_aft$nCounts)
        
        x_data <- object_aft@meta.data[, "Joint"]
        med    <- median(x_data)
        MAD    <- mad(x_data, center = med, na.rm = TRUE)
        lower  <- med -  2.5* MAD
        higher <- med + 2.5 * MAD
        n_low  <- which(x_data < lower)
        n_high <- which(x_data > higher)
        
        object_aft$Out <- rep(1,length(object_aft$Joint))
        object_aft$Out[c(n_low,n_high)] <- 2

        
        
        plot2 <- FeatureScatter(object_aft, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "Out")
        pdf(paste0(condition[1],"_Joint.pdf"),width=10)
        plot(outlierHistogram(object_aft@meta.data,x="Joint",mad = 1:4))
        plot(plot2)
        dev.off()
        
        # ==================================== Normalization ================================== #
        object_aft <- SCTransform(object = object_aft, verbose = FALSE)
        # -----------------------------------------------------------------------------------------------------
        return(object_aft)
})


saveRDS(tm_list_fl,"Fl_Object_list.rds")

setwd("../")