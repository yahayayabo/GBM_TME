
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



