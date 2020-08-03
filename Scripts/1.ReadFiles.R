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





# ================================== 1.ReadFiles ===========================================
dir.create("/home/users/dkyriakis/PhD/Projects/epi-scRNA/2.Preprocess")
setwd("/home/users/dkyriakis/PhD/Projects/epi-scRNA/2.Preprocess")

tm_list <-lapply(c(1:length(condition_names)), FUN=function(x) {
    print(filenames[x])
    gene <- Read10X(filenames[x], gene.column=2)
    batch <- condition_names[x]
    temp <- CreateSeuratObject(counts = gene, project = batch)

    cat(head(rownames(temp)))
    
    temp <- NormalizeData(temp)
    all.genes <- rownames(temp)
    temp <- ScaleData(temp, features = all.genes)
    # ======================================== CELL Cycle ==========================================
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    temp<- CellCycleScoring(temp, s.features = s.genes, g2m.features = g2m.genes)
    # ----------------------------------------------------------------------------------------------


    # ======================== HOUSE KEEPING GENES ==============================================
    hkgenes <- read.table("/home/users/dkyriakis/PhD/Projects/Yahaya/data/tirosh_house_keeping.txt", skip = 2)
    hkgenes <- as.vector(hkgenes$V1)    
    # remove hkgenes that were not found
    hkgenes.found <- which(toupper(rownames(temp@assays$RNA@data)) %in% hkgenes)
    # Add_number_of_house_keeping_genes
    n.expressed.hkgenes <- Matrix::colSums(temp@assays$RNA@data[hkgenes.found, ] > 0)
    temp <- AddMetaData(object = temp, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")
    # ----------------------------------------------------------------------------------------------


    # ======================================== Dissociation Genes ==========================================
    # Genes upregulated during dissociation of tissue into single cells.
    genes.dissoc <- c("ATF3", "BTG2", "CEBPB", "CEBPD", "CXCL3", "CXCL2", "CXCL1", "DNAJA1", "DNAJB1", "DUSP1", "EGR1", "FOS", "FOSB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "HSPE1", "HSPH1", "ID3", "IER2", "JUN", "JUNB", "JUND", "MT1X", "NFKBIA", "NR4A1", "PPP1R15A", "SOCS3", "ZFP36")
    #### seurat <- AddModuleScore(?, genes.list = list(?), ctrl.size = 20, enrich.name = "genes_dissoc")
    Seurat <- AddModuleScore(temp, features = list(genes.dissoc), ctrl.size = 20, enrich.name = "genes_dissoc",name="genes_dissoc")
    # ----------------------------------------------------------------------------------------------


    temp[["percent.mito"]]<-PercentageFeatureSet(temp,pattern="^MT-")
    #temp[["percent.rb"]]<-PercentageFeatureSet(temp,pattern="^RPL-|^RPS-")
    temp$nFeatures <- temp$nFeature_RNA
    temp$nCounts <- temp$nCount_RNA
    temp$Condition<- temp$orig.ident
    temp$condition<- temp$orig.ident
    temp$log10_total_Features <- log10(temp$nFeature_RNA)
    temp$log10_total_counts <- log10(temp$nCount_RNA)

    pdf(paste0(Sys.Date(),"_",batch,"_QC.pdf"))
    plot(VlnPlot(temp,c("nCount_RNA","nFeature_RNA","percent.mito")))
    dev.off()

    cell_data <- temp@meta.data
    p1h <- outlierHistogram(cell_data, "log10_total_Features", mads = 1:4)
    p2h <-outlierHistogram(cell_data, "log10_total_counts", mads = 1:4)
    p3h <-outlierHistogram(cell_data, "percent.mito", mads = 1:4)
    qc2 <- ggpubr::ggarrange(plotlist = list(p1h,p2h,p3h),nrow = 1)
    pdf(paste("1.QC",batch,".pdf"),width=12)
    plot(qc2)
    dev.off()
    return(temp)
})

saveRDS(tm_list,"Object_list.rds")

setwd("../")
# ---------------------------------------------------------------------------------


