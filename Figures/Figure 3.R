
#### Figure 3     ####

library(Seurat)
library(harmony)
library(dittoSeq)
library(SCENIC)
library(ggplot2)



#3a

# subset myeloid cluster from the PDOX_GL261_TME object

Myeloid_subset <- subset(PDOX_GL261_TME, subset=Celltype%in%'Myeloid')

Myeloid_subset <- ScaleData(Myeloid_subset, verbose = FALSE)
Myeloid_subset <- RunPCA(Myeloid_subset, features = VariableFeatures(object = Myeloid_subset))

print(Myeloid_subset[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Myeloid_subset, dims = 1:2, reduction = "pca")
DimPlot(Myeloid_subset, reduction = "pca", group.by = "sample")
DimHeatmap(Myeloid_subset, dims = 1, cells = 500, balanced = TRUE)

JackStrawPlot(Myeloid_subset, dims = 1:10)
ElbowPlot(Myeloid_subset)
Myeloid_subset <- Myeloid_subset %>% 
  RunHarmony("sample", plot_convergence = TRUE)



Myeloid_subset <- Myeloid_subset %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  RunTSNE(reduction = "harmony", dims = 1:10) %>%
  FindNeighbors(reduction = "harmony", dims = 1:6) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(Myeloid_subset, reduction = "umap", label = TRUE, pt.size = 1)


DimPlot(Myeloid_subset, reduction = "umap", group.by = "clusters", label = F, pt.size = 1, 
        cols = c('CL0'= "#22830B", 'CL1'= "#2CAA0E", 'CL2' = "#9EF57B", 
                 'CL3'= "#99B953",'CL4'= "#FFD966", 'CL5'= "#F08C06", 
                 'CL6' = "#7F6000", 'CL7'= "blue", 'CL8'= "#e31a1c"))



#3b

dittoBarPlot(Myeloid_subset,"clusters", group.by = "sample", 
             var.labels.reorder = c(9, 8, 7, 6, 5, 4, 3, 2, 1),
             x.reorder = c(5, 8, 9, 13, 14, 11, 12, 7, 10, 6, 4,2,1,3), 
             color.panel = c('CL8'= "#e31a1c",'CL7'= "blue",'CL6' = "#7F6000", 
                             'CL5'= "#F08C06",'CL4'= "#FFD966", 'CL3'= "#99B953", 
                             'CL2' = "#9EF57B",'CL1'= "#2CAA0E",'CL0'= "#22830B")) +
  theme(text = element_text(size = 16))


#3c  

#heatmap of z-score of expression of key marker genes for myeloid clusters 

myeloid_markers <- FindAllMarkers(Myeloid_subset, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
write.csv(myeloid_markers, file = "Table_S4.csv")

myeloid_obj <- readRDS("Myeloid_subset.rds")
subtype_name <- "Cluster_markers"

# list of key representative genes from each cluster obtained from the list of myeloid_markers after applying the threshold: FDR <=0.01 and log2FC>=0.5

genes <- c("Egr1",  "Fos",  "Ier5",  "Btg2","Zfp36","Jun","Hspa1a", "Siglech","Jund", "Rhob", "P2ry12", "Tmem119","Gpr34", 
           "Frmd4a", "Malat1", "Pnn", "Akap9", "Atrx","Dock10", "Cxcr1","Ssh2","Prpf4b", "Jmjd1c", "Sfrs18","Smarca2", 
           "mt.Nd1","Cxcl13", "Cst7",  "Ccl6", "Apoe", "Spp1", "Tyrobp","Cd52","Cd9","Trem2",  
           "Ctsd","Ctsz","Ctsl","B2m", "Ctsb", "Clec7a","Ccl2", "Ccl12","Cst7", "Ccl3", "Ccl4",
           "Sparcl1","Gfap", "Clu", "Slc1a2", "Ttr","Ttyh1", "Atp1a2", "Aldoc","Cpe", "Gja1","Ptn", "Ptprz1", 
           "Rps12", "Rps29","Rps27","Rpl17", "Gnb2l1", "Rps28","Rps23", "Ptplb","Apoc1", "Fau","Cd74", "Sepp1", 
           "Esam","Ly6c1","Pecam1", "Cldn5","Itm2a", "Igfbp7", "Flt1", "Ptprb","Spock2", "Slco1c1","Ctla2a","Ly6a","Crip1",
           "Ccr2","Ly6c2", "Plac8", "H2.Eb1","Il1b", "Thbs1", "Lgals3","S100a11","S100a10","Lyz2", "Tgfbi", 
           "Mrc1", "F13a1","Dab2", "Pf4","Ms4a7","Cd163","Ccl7", "Ccl8", "Cp", "Stab1","Clec10a")



DefaultAssay(myeloid_obj) <- "RNA" 
Idents(myeloid_obj) <- "clusters"

markers.df <- genes
markers.df <- unique(markers.df) %>% as.data.frame()
markers.df <- markers.df %>% tibble::rownames_to_column(var = "GeneOrder")
colnames(markers.df) <- c("GeneOrder", "Gene")
markers.df$Gene <- as.character(markers.df$Gene)
markers <- markers.df$Gene

markers <- markers[markers %in% rownames(myeloid_obj)]
expr.df <- AverageExpression(myeloid_obj, assays = "RNA", features = markers, slot = "data")
expr.df <- as.data.frame(expr.df)

expr.long <- expr.df %>%
  tibble::rownames_to_column(var = "Gene") %>%
  left_join(markers.df, by = "Gene")
expr.long <- expr.long %>% 
  tidyr::gather(key = cell.type , value = expr, -Gene,-GeneOrder) %>%
  group_by(Gene) %>%
  dplyr::mutate(z_score = (expr - mean(expr))/(sd(expr) + 0.01))

expr.long$GeneOrder <- as.numeric(expr.long$GeneOrder)

myeloid_heatmap <- ggplot(data = expr.long) +
  geom_tile(aes(x = cell.type, y = reorder(Gene,-GeneOrder),fill = z_score), color = "white", size = 0.1) + 
  scale_fill_gradient2(low = "#3182bd",mid= "white", high = "darkred") +
  scale_y_discrete(position = "right") +
  xlab("") + ylab("") +
  theme_grey(base_size = 10) + 
  ggtitle(subtype_name) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10, colour = "gray50"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 10))


ggplot(expr.long)+
  geom_tile(aes(x = cell.type, y = reorder(Gene,-GeneOrder),fill = z_score), color = "white", size = 0.1)+
  scale_fill_gradient2(low = "#3182bd",mid= "white", high = "darkred") +
  scale_y_discrete(position = "right") +
  plot_theme +
  theme(axis.text.y.right = element_text(angle = 0, hjust=0, size=10),
        legend.position="right")



#3e

# List of transcription factors (TF) was generated as described in SCENIC vignette http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html

### Initialize settings
library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
library(GENIE3)
library(doSNOW)
library(doParallel) 
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(matrixTests)
library(ggpubr)
library(ggrepel)

# upload mouse DB 
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")


exprMat <- rawData
cellInfo <- data.frame(cellInfo)
cellTypeColumn <- "Annotation"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cbind(table(cellInfo$CellType))


colVars <- list(CellType=c("CL0" = "#e5f5e0",
                           "CL1" = "#c7e9c0",
                           "CL2" = "#a1d99b",
                           "CL3" = "#74c476",
                           "CL4" = "#41ab5d",
                           "CL5" = "#238b45",
                           "CL6" = "#005a32",
                           "CL7" = "blue",
                           "CL8" = "#e31a1c"))


colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]

plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))



org="mgi" # or hgnc, or dmel
dbDir="C:/Users/yabubakaryabo/Documents/cisTarget_databases" # RcisTarget databases location
myDatasetTitle="SCENIC myeloid cells" # choose a name for your analysis
data(defaultDbNames)
dbs=c("mm9-500bp-upstream-7species.mc9nr.feather", "mm9-tss-centered-10kb-7species.mc9nr.feather")

### Initialize settings
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=1)


scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions)

interestingGenes <- c("Sparc", "Tgfbi", "Hexb", "Cd74", "Ccl5")
interestingGenes[which(!interestingGenes %in% genesKept)]
exprMat_filtered <- exprMat[genesKept, ]


# Correlation 

runCorrelation(exprMat_filtered, scenicOptions) 
exprMat_filtered <- log2(exprMat_filtered+1) 

runGenie3(exprMat_filtered, scenicOptions) 

#Build and score the GRN
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- 123

logMat <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) 
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, logMat)


# Binarize the network activity (regulon on/off)  

aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# heatmap of enriched TF was done with modifications as described here https://github.com/TheJacksonLaboratory/singlecellglioma-verhaaklab/blob/master/analysis/Fig3c-scenic-IDHmut.R

#load data

myeloid_obj <- readRDS("Myeloid_subset.rds")
myeloid_obj_DGEM  <- myeloid_obj@assays$RNA@counts

metadata <-  as.data.frame(myeloid_obj[["umap"]]@cell.embeddings)
myeloid_obj <- AddMetaData(myeloid_obj, metadata$UMAP_1, col.name = "UMAP_1")
myeloid_obj <- AddMetaData(myeloid_obj, metadata$UMAP_2, col.name = "UMAP_2")

umap_coords_2d <- read.csv("Subset_annotation.csv", header = TRUE)
View(umap_coords_2d)

# Limit to sample of interest
umap_data_mut <- umap_coords_2d %>%
  filter(case_barcode %in% c("MGBM11", "MGBM13", "MGBM15", "MGBM17", "Normal","P13", "P3CTR", "P8", "T101", "T16","T192", "T233", "T347", "T470"))


## Keep only cell type of interest
cells_to_keep = which(umap_coords_2d$cell_state %in%c("CL0", "CL1", "CL2", "CL3", "CL4", "CL5", "CL6", "CL7", "CL8"))
clust_annot <- umap_coords_2d[cells_to_keep, ]

# Extract the exact name of cells.
cell_names_keep <- clust_annot$cell_barcode


expr_norm_data <- myeloid_obj_DGEM[ ,colnames(myeloid_obj_DGEM)%in%cell_names_keep]
all(colnames(expr_norm_data)==clust_annot$cell_barcode)
expr_norm_data_sample <- expr_norm_data
clust_annot_sample <- clust_annot

# Make sure that the same cells are being subsetted.
all(colnames(expr_norm_data_sample)==clust_annot_sample$cell_barcode)

## Use broad tumor cell classification.
myeloid_obj@meta.data$sample <- clust_annot_sample$sample
myeloid_obj@meta.data$New_names <- clust_annot_sample$New_names

# Re-load SCENIC results 
auc_rankings <- readRDS(".../3.3_aucellRankings.Rds")
regulonAUC <- readRDS(".../3.4_regulonAUC.Rds")

## Create a data.frame with the gene sets/TFs and cells.
regulonAUC_df = as.data.frame(getAUC(regulonAUC))

## The annotation files we have match the regulonAUC data.
all(clust_annot_sample$cell_barcode==colnames(regulonAUC_df))

## generate z-scores for variable A using the scale() function
## scale(A, center = TRUE, scale = TRUE). These are the defaults.
regulonAUC_scaled = t(apply(as.matrix(regulonAUC_df), 1, scale))

## Provide the case_barcode, cell_state, and subtype annotations for each cell.
annot_df = data.frame(clust_annot_sample$cell_barcode, rep("Myeloid", dim(clust_annot_sample)[1]), clust_annot_sample$case_barcode, clust_annot_sample$cell_state)
colnames(annot_df) <- c("barcode","New_names", "case_barcode", "cell_state")
epimut_cols <- colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1.0), c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"))

## Define the annotation colors:
ha = HeatmapAnnotation(df = annot_df,
                       col = list(case_barcode = c("MGBM11" = "#a6cee3",
                                                   "MGBM13" = "#1f78b4",
                                                   "MGBM15" = "#b2df8a",
                                                   "MGBM17" = "#33a02c",
                                                   "Normal" = "#fb9a99",
                                                   "P13" = "#e31a1c",
                                                   "P3CTR" = "#fdbf6f",
                                                   "P8" = "#ff7f00",
                                                   "T101" = "#cab2d6",
                                                   "T16" = "#6a3d9a",
                                                   "T192" = "#ffff99",
                                                   "T233" = "#b15928",
                                                   "T347" = "#ffed6f",
                                                   "T470" = "#fb8072"),
                                  cell_state = c("CL0" = "#e5f5e0",
                                                 "CL1" = "#c7e9c0",
                                                 "CL2" = "#a1d99b",
                                                 "CL3" = "#74c476",
                                                 "CL4" = "#41ab5d",
                                                 "CL5" = "#238b45",
                                                 "CL6" = "#005a32",
                                                 "CL7" = "blue",
                                                 "CL8" = "#e31a1c"),
                                  subtype = c("Myeloid" = "#AF8DC3")))

## What are the TFs considered in this analysis?
scenic_tfs = data.frame(tf = sapply(strsplit(rownames(regulonAUC_scaled), "_| "), "[[", 1))


cell_state <- clust_annot_sample$cell_state

## Apply kruskal.wallis across three cell types. Select X number of TFs
kw_test_results = row_kruskalwallis(regulonAUC_scaled, clust_annot_sample$cell_state)
kw_test_results$adj_pvalue = p.adjust(kw_test_results$pvalue, method = "fdr", n = length(kw_test_results$pvalue))
kw_test_results_diff = kw_test_results[which(kw_test_results$pvalue< 1e-127), ]
tfs_to_keep = rownames(kw_test_results_diff)

## Remove any "_extended" TFs that are also represented.
tmp = sapply(strsplit(rownames(regulonAUC_scaled), "_| "), "[[", 1) %>% as.data.frame()
colnames(tmp) <- "tf"
duplicated_tf = tmp %>%
  group_by(tf) %>%
  summarise(tf_counts = n()) %>%
  filter(tf_counts >= 2)
dup_tfs_remove = rownames(regulonAUC_scaled)[which(sapply(strsplit(rownames(regulonAUC_scaled), "_| "), "[[", 1)%in%duplicated_tf$tf & grepl("_extended", rownames(regulonAUC_scaled)))]
tfs_to_keep_uniq = rownames(regulonAUC_scaled)[!rownames(regulonAUC_scaled)%in%dup_tfs_remove]

## Take the 15 most enriched TFs per cell state. Some may overlap leading there to be fewer than 45 across the two cell states
CL0_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="CL0"], 1, median), decreasing = TRUE)
CL1_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="CL1"], 1, median), decreasing = TRUE)
CL2_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="CL2"], 1, median), decreasing = TRUE)
CL3_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="CL3"], 1, median), decreasing = TRUE)
CL4_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="CL4"], 1, median), decreasing = TRUE)
CL5_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="CL5"], 1, median), decreasing = TRUE)
CL6_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="CL6"], 1, median), decreasing = TRUE)
CL7_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="CL7"], 1, median), decreasing = TRUE)
CL8_rank = sort(apply(regulonAUC_scaled[tfs_to_keep_uniq, cell_state=="CL8"], 1, median), decreasing = TRUE)

ranked_tfs_to_keep = unique(c( names(CL0_rank[1:15]), names(CL1_rank[1:15]), names(CL2_rank[1:15]),
                               names(CL3_rank[1:15]), names(CL4_rank[1:15]),names(CL5_rank[1:15]),
                               names(CL6_rank[1:15]), names(CL7_rank[1:15]), names(CL8_rank[1:15])))

## Filter the regulonAUC plot.
regulonAUC_scaled_filt = regulonAUC_scaled[ranked_tfs_to_keep, ]
rownames(regulonAUC_scaled_filt) <- gsub("_extended", "", rownames(regulonAUC_scaled_filt))

## Gather activity by cell state.
Cl_0_df <- as.data.frame(CL0_rank)
Cl_0_df$tf <- sapply(strsplit(rownames(Cl_0_df), " "), "[[", 1)
Cl_0_df$tf  <- reorder(Cl_0_df$tf , Cl_0_df$CL0_rank)

Cl_1_df <- as.data.frame(CL1_rank)
Cl_1_df$tf <- sapply(strsplit(rownames(Cl_1_df), " "), "[[", 1)
Cl_1_df$tf  <- reorder(Cl_1_df$tf , Cl_1_df$CL1_rank)

Cl_2_df <- as.data.frame(CL2_rank)
Cl_2_df$tf <- sapply(strsplit(rownames(Cl_2_df), " "), "[[", 1)
Cl_2_df$tf  <- reorder(Cl_2_df$tf , Cl_2_df$CL2_rank)

Cl_3_df <- as.data.frame(CL3_rank)
Cl_3_df$tf <- sapply(strsplit(rownames(Cl_3_df), " "), "[[", 1)
Cl_3_df$tf  <- reorder(Cl_3_df$tf , Cl_3_df$CL3_rank)

Cl_4_df <- as.data.frame(CL4_rank)
Cl_4_df$tf <- sapply(strsplit(rownames(Cl_4_df), " "), "[[", 1)
Cl_4_df$tf  <- reorder(Cl_4_df$tf , Cl_4_df$CL4_rank)

Cl_5_df <- as.data.frame(CL5_rank)
Cl_5_df$tf <- sapply(strsplit(rownames(Cl_5_df), " "), "[[", 1)
Cl_5_df$tf  <- reorder(Cl_5_df$tf , Cl_5_df$CL5_rank)

Cl_6_df <- as.data.frame(CL6_rank)
Cl_6_df$tf <- sapply(strsplit(rownames(Cl_6_df), " "), "[[", 1)
Cl_6_df$tf  <- reorder(Cl_6_df$tf , Cl_6_df$CL6_rank)

Cl_7_df <- as.data.frame(CL7_rank)
Cl_7_df$tf <- sapply(strsplit(rownames(Cl_7_df), " "), "[[", 1)
Cl_7_df$tf  <- reorder(Cl_7_df$tf , Cl_7_df$CL7_rank)

Cl_8_df <- as.data.frame(CL8_rank)
Cl_8_df$tf <- sapply(strsplit(rownames(Cl_8_df), " "), "[[", 1)
Cl_8_df$tf  <- reorder(Cl_8_df$tf , Cl_8_df$CL8_rank)

### Create a matrix of summary values:
median_activity_df <- Cl_0_df %>%
  inner_join(Cl_1_df, by="tf") %>%
  inner_join(Cl_2_df, by="tf") %>%
  inner_join(Cl_3_df, by="tf") %>%
  inner_join(Cl_4_df, by="tf") %>%
  inner_join(Cl_5_df, by="tf") %>%
  inner_join(Cl_6_df, by="tf") %>%
  inner_join(Cl_7_df, by="tf") %>%
  inner_join(Cl_8_df, by="tf") %>%
  select(tf, CL0_rank, CL1_rank, CL2_rank, CL3_rank, CL4_rank, CL5_rank, CL6_rank, CL7_rank, CL8_rank)
median_activity_df_filt <- median_activity_df %>%
  pivot_longer(cols= c(CL0_rank:CL8_rank),
               names_to = "cell_state",
               values_to = "activity") %>%
  mutate(cell_state = recode(cell_state, 
                             `CL0_rank` = "CL0",
                             `CL1_rank` = "CL1",
                             `CL2_rank` = "CL2",
                             `CL3_rank` = "CL3",
                             `CL4_rank` = "CL4",
                             `CL5_rank` = "CL5",
                             `CL6_rank` = "CL6",
                             `CL7_rank` = "CL7",
                             `CL8_rank` = "CL8")) %>%
  mutate(tf = gsub("_extended", "", tf)) %>%
  filter(tf%in%c(sapply(strsplit(rownames(regulonAUC_scaled_filt), "_| "), "[[", 1)))

## Set the TF order to be ranked for each cell state.
tf_order <- sapply(strsplit(rownames(regulonAUC_scaled_filt), "_| "), "[[", 1)
tf_order
#tf_reverse <- reverse(tf_order)
tf_reverse <- rev(tf_order)
tf_reverse

#tf_levels = tf_order
tf_levels = tf_reverse
state_levels = c("CL8","CL7","CL6","CL5","CL4","CL3","CL2", "CL1", "CL0")
median_activity_df_filt$tf <-  factor(median_activity_df_filt$tf, levels = tf_levels)
median_activity_df_filt$cell_state <-  factor(median_activity_df_filt$cell_state, levels = rev(state_levels))


##Heatmap
palette     <- c( "#8CB9DA" , "white", "#A43838", "#8B0000")
ggplot(data = median_activity_df_filt) +
  geom_tile(aes(x = cell_state, y = tf, fill=activity), color = "white", size = 0.1) + # bigger size means bigger spacer
  scale_fill_gradientn(colours=c(palette))  +
  scale_y_discrete(position = "right") +
  labs(x="", y= "", fill="Relative TF activity\n(Z-score)")  +
  theme_grey(base_size = 10) + # text size
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10, colour = "gray50"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 10))







