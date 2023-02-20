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

library(EnhancedVolcano)
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



test <- readRDS("/home/users/dkyriakis/PhD/Projects/epi-scRNA/4.Clustering/SCT_Merged_Cl.rds")




setwd("5.DF")

DefaultAssay(test) <- "RNA"
Idents(test) <- test$Cell_Type
neu.cells <- subset(test, subset= Cell_Type %in% c("Neu"))
int.cells <- subset(test, subset= Cell_Type %in% c("Interneurons"))
oligo.cells <- subset(test, subset= Cell_Type %in% c("Oligo","OPC"))
astro.cells <- subset(test, subset= Cell_Type %in% c("AST-FB","AST-PP"))

subs_names <- c("Neurons","Interneurons","Oligodendrocytes","Astrocytes")
subs_data <- list(neu.cells,int.cells,oligo.cells,astro.cells)

for (iter in 1:length(subs_names)){
  name_sub <- subs_names[iter]
  data_sub <- subs_data[[iter]]
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  Idents(data_sub) <- data_sub$condition
  pbmc.markers <- FindAllMarkers(data_sub, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0,return.thresh=1)
   
  top15 <- pbmc.markers[pbmc.markers$p_val_adj <0.01,] %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
  
  write.table(pbmc.markers,paste0(name_sub,"_DF.tsv"))
  write.table(top15,paste0("top_",name_sub,"_DF.tsv"))
  
  data_sub <- ScaleData(data_sub,features=top15$gene)
  pdf(paste0(name_sub,"_DF.pdf"))
  plot(DoHeatmap(data_sub, features = top15$gene) + NoLegend())
  fc_alt <- pbmc.markers$avg_logFC
  fc_alt[pbmc.markers$cluster=="epilepsy"] <- -fc_alt[pbmc.markers$cluster=="epilepsy"] 
  res1 <- data.frame("log2FoldChange"=fc_alt,"adj_pvalue"=pbmc.markers$p_val_adj)
  rownames(res1) <- pbmc.markers$gene
  plot(EnhancedVolcano(res1,
      lab = rownames(res1),
      x = 'log2FoldChange',
      y = 'adj_pvalue',
      pCutoff = 0.01,
      FCcutoff = 0.6))
  dev.off()
  saveRDS(data_sub,paste0(name_sub,".rds"))
}


setwd("../")

