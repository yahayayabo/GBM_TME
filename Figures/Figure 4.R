
####   Figure 4     #####

library(Seurat)
library(dittoSeq)
library(tidyverse)
library(ggplot2)

#load data
seurat_obj <- readRDS("myeloid_object.rds")


#Fig 4A
DimPlot(seurat_obj, reduction = "umap", group.by = "New_Clusters", label = F, pt.size = 1, 
        cols = c('CL0'= "#22830B", 'CL1'= "#2CAA0E", 'CL2' = "#9EF57B", 'CL3'= "#99B953",'CL4'= "#FFD966", 
                 'CL5' = "#7F6000", 'CL6'= "#F08C06" , 'CL7'= "blue", 'CL8'= "#e31a1c"))

#Fig 4B
library(dittoseq)
dittoBarPlot(seurat_obj,"New_Clusters", group.by = "sample", 
             var.labels.reorder = c(9, 8, 7, 6, 5, 4, 3, 2, 1),
             x.reorder = c(5, 8, 9, 13, 14, 11, 12, 7, 10, 6, 4,2,1,3), 
             color.panel = c('CL8'= "#e31a1c",'CL7'= "blue",'CL6' = "#F08C06", 
                             'CL5'= "#7F6000",'CL4'= "#FFD966", 'CL3'= "#99B953", 
                             'CL2' = "#9EF57B",'CL1'= "#2CAA0E",'CL0'= "#22830B")) +
  theme(text = element_text(size = 16))


#Fig 4C
#find cluster markers
Idents(seurat_obj) <- "New_Clusters"
New_cluster_markers_ordered <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(New_cluster_markers_ordered, "New_cluster_markers_ordered_DEGS.csv", quote = FALSE,  row.names = TRUE)

markers_for_heatmap<- New_cluster_markers_ordered %>% 
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05)

object <-  Myeloid_subset
subtype_name <- "New_Clusters"

genes <- c( "Egr1","Ier5","Btg2", "Fos",  "P2ry12", "Gpr34", "Tmem119", "Rhob",  "Jund","Jun",  #CL0
            "Prpf4b", "Pnn","Luc7l3","Luc7l2","Mycbp2","Malat1", "Arid4a", "Rbm5","Zfhx3","Mertk", #CL1
            "Ctsd", "Apoe", "Ctsb", "Cd63", "Lpl", "Ccl6", "Creg1", #CL2
            "Cst7","Tyrobp","Ch25h","Cd63","Cd52","Apoc1","Spp1","Cxcl13","Ccl3","Ccl4",#CL3
            "Sparcl1", "Gfap", "Clu", "Plpp3", "Atp1a2", "Scg3","Mt2","Cpe","Slc1a2","Aldoc", "Ptn", #CL4
            "Igfbp7", "Ly6c1", "Cldn5", "Itm2a", "Abcb1a", "Ptprb", "Adgrf5", "Spock2", "Bsg", "Cxcl12", #CL5
            "Mki67", "Top2a", "Cenpf", "Cenpe", "Prc1", "Hmmr", "Smc2", "Smc4","2810417H13Rik","Kif20b",  #CL6
            "Ly6c2", "Plac8", "Lgals3","Ms4a4c","Il1b","H2.Eb1","Ly6a","Tgfbi", #CL7
            "Mrc1", "F13a1", "Pf4", "Ms4a7", "Dab2", "Cp", "Wfdc17", "Wwp1", "Stab1", "Maf")#Cl8

# HEATMAP ----
DefaultAssay(object) <- "RNA" # you will need to use SCT if the RNA is not normalised and scaled
Idents(object) <- "New_Clusters"

markers.df <- genes
markers.df <- unique(markers.df) %>% as.data.frame()
markers.df <- markers.df %>% tibble::rownames_to_column(var = "GeneOrder")
colnames(markers.df) <- c("GeneOrder", "Gene")
markers.df$Gene <- as.character(markers.df$Gene)
markers <- markers.df$Gene

markers <- markers[markers %in% rownames(object)]
expr.df <- AverageExpression(object, assays = "RNA", features = markers, slot = "data")
expr.df <- as.data.frame(expr.df)

expr.long <- expr.df %>%
  tibble::rownames_to_column(var = "Gene") %>%
  left_join(markers.df, by = "Gene")
expr.long <- expr.long %>% 
  tidyr::gather(key = cell.type , value = expr, -Gene,-GeneOrder) %>%
  group_by(Gene) %>%
  dplyr::mutate(z_score = (expr - mean(expr))/(sd(expr) + 0.01))

expr.long$GeneOrder <- as.numeric(expr.long$GeneOrder)
# expr.long$cell.type <- gsub("RNA.", "", expr.long$cell.type) # don't know if you need this 

heatmap.s <- ggplot(data = expr.long) +
  geom_tile(aes(x = cell.type, y = reorder(Gene,-GeneOrder),fill = z_score), color = "white", size = 0.1) + # bigger size means bigger spacer
  scale_fill_gradient2(low = "#3182bd",mid= "white", high = "darkred") +
  scale_y_discrete(position = "right") +
  xlab("") + ylab("") +
  theme_grey(base_size = 10) + # text size
  ggtitle(subtype_name) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10, colour = "gray50"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) 

heatmap.s

Cairo::Cairo(2600, 5200, file = "heatmap.s.png", type = "png", #tiff
             bg = "transparent", #white or transparent depending on your requirement 
             dpi = 500)
plot(heatmap.s)
dev.off()

#Fig 4E

#############################################################################
########################     SCENIC Analysis  ################################
#############################################################################

projectFolder = "path"

#load seurat object
setwd("path")
seurat_obj <- readRDS("seurat_obj_final.rds")
dim(seurat_obj)
#[1] 24067  6977
table(seurat_obj$New_Clusters)
#CL3  CL1  CL0  CL2  CL7  CL4  CL6  CL8  CL5 
#1438 1288 1171  863  792  731  353  184  157 

#the files needed
seurat_obj_DGEM  <- seurat_obj@assays$RNA@data
cell_metadata <- seurat_obj@meta.data
gene_metadata <- rownames(seurat_obj@assays$RNA@data)

#write.table(cell_metadata, file = "cell_metadata.txt", sep ="\t")
#write.table(gene_metadata, file = "gene_metadata.txt", sep ="\t")

#### upload mouse DB #######
#For mouse:
library(SCENIC)

dbFiles <- c("scenic_feather_files/mm9-500bp-upstream-7species.mc9nr.feather",
             "scenic_feather_files/mm9-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs

# Tutorial
vignette("SCENIC_Running") # open

if(!dir.exists("SCENIC_MouseMyeloidCells")){
  
  dir.create("SCENIC_MouseMyeloidCells")
  setwd("SCENIC_MouseMyeloidCells")
  
}else{
  
  setwd("SCENIC_MouseMyeloidCells")
  
}


rawData <- as.matrix(seurat_obj_DGEM)

cellInfo <- cell_metadata

dim(rawData)
#[1] 24067  6977

head(cellInfo)

exprMat <- rawData
cellInfo <- data.frame(cellInfo)
cellTypeColumn <- "New_Clusters"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cbind(table(cellInfo$CellType))


if(!dir.exists("int")){
  
  dir.create("int")
  
}


saveRDS(cellInfo, file="int/cellInfo.Rds")


colVars <- list(CellType=c("CL0" = "#22830B",
                           "CL1" = "#2CAA0E",
                           "CL2" = "#9EF57B",
                           "CL3" = "#99B953",
                           "CL4" = "#FFD966",
                           "CL5" = "#7F6000",
                           "CL6" = "#F08C06",
                           "CL7" = "blue",
                           "CL8" = "#e31a1c"))


colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]

plot.new();legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

saveRDS(colVars, file="int/colVars.Rds")


################################ Running your sample #############################################

### Initialize settings
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
library(GENIE3)

setwd(projectFolder)

org="mgi" # mgi or hgnc, or dmel
dbDir= "scenic_feather_files" # I have sent you the files
myDatasetTitle="SCENIC myeloid cells" # choose a name for your analysis
data(defaultDbNames)
dbs=c("mm9-500bp-upstream-7species.mc9nr.feather", 
      "mm9-tss-centered-10kb-7species.mc9nr.feather")


## https://github.com/aertslab/SCENIC/issues/364
motifAnnotations_mgi <- motifAnnotations


### Initialize settings
scenicOptions <- initializeScenic(org=org, 
                                  dbDir=dbDir, 
                                  dbs=dbs, 
                                  datasetTitle=myDatasetTitle, 
                                  nCores=4)
#Motif databases selected: 
#  mm9-500bp-upstream-7species.mc9nr.feather 
# mm9-tss-centered-10kb-7species.mc9nr.feather


scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

###############        Co-expression network      #######################################
#########################################################################################

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions)


interestingGenes <- c("Sparc", "Tgfbi", "Hexb", "Cd74", "Ccl5")

# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]
#0  means all interesting genes are present

exprMat_filtered <- exprMat[genesKept, ]

dim(exprMat_filtered)


#######################       Correlation        ###########################################
############################################################################################

runCorrelation(exprMat_filtered, scenicOptions) # takes a few minutes

########################     Run GENIE3         ###########################################

# GENIE3 will typically take several hours (or days)

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 


runGenie3(exprMat_filtered, scenicOptions) # takes many hours or days to finish running:   took little more than 3 days from Wed ~4pm to Sat ~23.30 pm (24.06.2023)
#Finished running GENIE3.
#Warning message:
#  In runGenie3(exprMat_filtered, scenicOptions) :
#  Only 500 (35%) of the 1412 TFs in the database were found in the dataset. Do they use the same gene IDs?


### saving the GENIE3 output results as RDS
saveRDS(exprMat_filtered, file="output/exprMat_filtered_scenicResults_after_GENIE3_run.Rds")
saveRDS(scenicOptions, file="output/scenicOptions_scenicResults_after_GENIE3_run.Rds")


########## Build and score the GRN      ############################################

#library(SCENIC)
#scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
#scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123

logMat <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) 



scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) #takes up to 4-5 hrs 
#Warning message:
#  In RcisTarget::addLogo(tableSubset, motifCol = motifCol, dbVersion = dbVersion,  :
#  There is no annotation version attribute in the input table (it has probably been loaded with an 
#  older version of the package).'v9' will be used as it was the old default,but we recommend to 
#  re-load the annotations and/or re-run the enrichment to make sure everything is consistent.

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, logMat)




#================================================================================
#Error in foreach(param = allParams) %dopar% { : 
#    could not find function "%dopar%"

# quick fix is to re-install these packages

library("doSNOW")
library("doParallel") 
#library("doMPI")
#=================================================================================
library(tidyverse)


saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status

####################    Binarize the network activity (regulon on/off)   ###############
# to create the "Binarized regulon activity matrix", which can be used for upstream analysis (e.g. clustering)
#To determine in which cells each regulon is active, we will use an AUC threshold

aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

runSCENIC_4_aucell_binarize(scenicOptions)
#Warning message:
#  In AUCell_plotTSNE(tSNE = tSNE, exprMat = exprMat, cellsAUC = regulonAUC,  :
#                       Expression plot was requested, but no expression matrix provided.
                     
saveRDS(scenicOptions, file="int/scenicOptions.Rds") # 

### from this stage you send me the last saved scenicOptions.Rds file:  Final results was sent to Yahaya on 26/06.


#### Output analysis:
##################################
library(Seurat)
library(SCENIC)
library(tidyverse)
library(SCopeLoomR)

projectFolder = "path"
setwd(projectFolder)
scenicOptions <- readRDS(file.path(projectFolder, "int/scenicOptions.Rds"))

seurat_obj <- readRDS("path/seurat_obj_final.rds")

cellInfo <- data.frame(seuratCluster=seurat_obj$New_Clusters)
cellInfo$seuratCluster <- factor(cellInfo$seuratCluster, levels = paste("CL", 0:8, sep = ""))

###
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

tmp <- split(rownames(cellInfo), cellInfo$seuratCluster)
regulonActivity_byCellType <- sapply(tmp, function(cells) {
                                            
                                            rowMeans(AUCell::getAUC(regulonAUC)[,cells])
                                            
                                      })

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")


###
# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=AUCell::getAUC(regulonAUC), 
               cellAnnotation=cellInfo[colnames(regulonAUC), "seuratCluster"])
rss <- rss[,paste("CL", 0:8, sep = "")]

gene.names <- rownames(rss)

gene.names.final <- unlist(lapply(gene.names, function(x) strsplit(x, " ")[[1]][1]))
gene.names.final <- unlist(lapply(gene.names.final, function(x) strsplit(x, "_")[[1]][1]))

## extended genes
tmp.1 <- grep("extended", gene.names)

## keep these genes
tmp.11 <- setdiff(c(1:nrow(rss)), tmp.1)
tmp.111 <- unlist(lapply(gene.names[tmp.11], function(x) strsplit(x, " ")[[1]][1]))

## some genes in extended list are not tmp.111
tmp.2 <- unlist(lapply(gene.names[tmp.1], function(x) strsplit(x, "_")[[1]][1]))

## genes with only "_extended"
tmp.3 <- setdiff(tmp.2, tmp.111)
#maf, mafb, maff, mafk
tmp.33 <- unlist(lapply(tmp.3, function(x, y) grep(paste0(x, "$"), y, perl = T), gene.names.final))

## combine
tmp.final <- sort(c(tmp.11, tmp.33))

gene.names.final <- gene.names.final[tmp.final]

rss <- rss[tmp.final,]
rownames(rss) <- gene.names.final

## 
rssPlot <- plotRSS(rss, thr = 0.05)
rssPlot$plot$data$cellType <- factor(rssPlot$plot$data$cellType, levels = paste("CL", 0:8, sep = ""))

plotly::ggplotly(rssPlot$plot) + theme(axis.text.y = element_text(size = 14))

ggsave(file.path(projectFolder, "SCENIC_Z_and_RSS_scores_plot.png"),
       width = 5,
       height = 10,
       dpi = 200)

##
rssPlot <- plotRSS(rss, thr = 0.05, revCol = T)
rssPlot$plot$data$cellType <- factor(rssPlot$plot$data$cellType, levels = paste("CL", 0:8, sep = ""))

plotly::ggplotly(rssPlot$plot) + theme(axis.text.y = element_text(size = 14))

ggsave(file.path(projectFolder, "SCENIC_Z_and_RSS_scores_plot_reverseColor.png"),
       width = 5,
       height = 10,
       dpi = 200)



