
####  Figure 7  #####

library(Seurat)
library(BiocParallel)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(dittoSeq)
library(Matrix)
library(dplyr)
library(CellChat)



### Fig 7C
# load P3 CTR and TMZ treated tumor cells

P3_CTR_TMZ_tumor <- readRDS("~/P3_CTR_TMZ_tumor.rds")

DimPlot(P3_CTR_TMZ_tumor, group.by =  "orig.ident", pt.size = 1,  cols = c("P3CTR" = "black" , 'P3TMZ' = 'red'))



###Fig 7f

# reference mapping to add P3TMZ to the myeloid subset

myeloid.reference <- SCTransform(Myeloid_subset, verbose = FALSE)
myeloid.query <- SCTransform(P3TMZ_TME, verbose = FALSE)


anchors <- FindTransferAnchors(
  reference = myeloid.reference,
  query = myeloid.query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:20
)

predictions <- TransferData(anchorset = anchors, refdata = myeloid.reference$New_names, dims = 1:20)

myeloid.reference <- RunUMAP(myeloid.reference, dims = 1:20, reduction = "pca", return.model = TRUE)


Myeloid_mapped <- MapQuery(
  anchorset = anchors,
  query = myeloid.query,
  reference = myeloid.reference,
  refdata = list(celltype = "New_names"), reference.reduction = "pca", reduction.model = "umap")



CTR <- subset(Myeloid_subset, subset = sample == "P3CTR")

P3.merged <- merge(CTR, y= Myeloid_mapped, project = "TMZ")

### stack chart of myeloid clusters in P3 TMZ and CTR

dittoBarPlot(P3.merged,"New_names", group.by = "sample", 
             var.labels.reorder = c(8, 7, 6, 5, 4, 3, 2, 1),
             color.panel = c('CL8'= "#e31a1c",'CL7'= "blue",'CL6' = "#7F6000", 
                             'CL5'= "#F08C06",'CL4'= "#FFD966", 'CL3'= "#99B953", 
                             'CL2' = "#9EF57B",'CL1'= "#2CAA0E",'CL0'= "#22830B")) +
  theme(text = element_text(size = 16))


### Fig 7G

DotPlot(P3.merged, features = c("APCs_score", "Phagocytosis_score", "Sensome_score", "Migration_score"), col.min = -2,
        col.max = 2, split.by = "sample", group.by = "Clusters", scale = TRUE, dot.scale = 10, cols = c("#3182bd","red")) +   
  RotatedAxis()+coord_flip()


#Fig 7H, I, J and K

# this analysis was done following the CellChat vignette with some modifications https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html


library(NMF)
library(ggalluvial)
library(ComplexHeatmap)
library(nichenetr)
library(ggplot2)                  
library(patchwork)
library(igraph)

#convert human genes to mouse in P3 tumor object

#P3CTR_tumor
P3CTR_tumor_count_matrix <- as.matrix(P3CTR_tumor@assays$RNA@counts)
rownames(P3CTR_tumor_count_matrix) = rownames(P3CTR_tumor_count_matrix) %>% convert_human_to_mouse_symbols() 
P3CTR_tumor_count_matrix = P3CTR_tumor_count_matrix %>% .[!is.na(rownames(P3CTR_tumor_count_matrix)), !is.na(colnames(P3CTR_tumor_count_matrix))]
P3CTR_tumor <- CreateSeuratObject(counts = P3CTR_tumor_count_matrix, project = "P3CTR_tumor")

#P3TMZ_tumor
P3TMZ_tumor_count_matrix <- as.matrix(P3TMZ_tumor@assays$RNA@counts)
rownames(P3TMZ_tumor_count_matrix) = rownames(P3TMZ_tumor_count_matrix) %>% convert_human_to_mouse_symbols() 
P3TMZ_tumor_count_matrix = P3TMZ_tumor_count_matrix %>% .[!is.na(rownames(P3TMZ_tumor_count_matrix)), !is.na(colnames(P3TMZ_tumor_count_matrix))]
P3TMZ_tumor <- CreateSeuratObject(counts = P3TMZ_tumor_count_matrix, project = "P3TMZ_tumor")


#merge tumor and TME data

P3CTR.full <- merge(P3CTR_TME, y= P3CTR_tumor, add.cell.ids = c("P3CTR_TME", "P3CTR_tumor"), project = "P3CTR.combined", merge.data = TRUE)
P3TMZ.full <- merge(P3TMZ_TME, y= P3TMZ_tumor, add.cell.ids = c("P3TMZ_TME", "P3TMZ_tumor"), project = "P3TMZ.combined", merge.data = TRUE)

# P3CTR tumor and TME
### change gene names to standard names with "-" instead of with "."
VlnPlot(P3CTR.full, features = "H2.Aa")
which(P3CTR.full@assays$RNA@data@Dimnames[[1]] == "H2.Aa")
P3CTR.full@assays$RNA@data@Dimnames[[1]] <- sub('\\.','-', P3CTR.full@assays$RNA@data@Dimnames[[1]])
VlnPlot(P3CTR.full, features = "H2-Aa")
which(P3CTR.full@assays$RNA@data@Dimnames[[1]] == "H2-Aa")

#Extract the CellChat input files from a Seurat object

data.input <- GetAssayData(P3CTR.full, assay = "RNA", slot = "data") # normalized data matrix
Idents(P3CTR.full) <- "celltypes"
labels <- Idents(P3CTR.full)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#Extract the CellChat input files from a Seurat object
cellchat.obj <- createCellChat(object = P3CTR.full, meta = meta, group.by = "group")


#Set the ligand-receptor interaction database        
CellChatDB  <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat.obj@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat.obj)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.mouse)

#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 0)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
#Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram

pathways.show <- c("GAS") 
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


netAnalysis_contribution(cellchat, signaling = pathways.show)


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:10), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 10, targets.use = c(1:9), remove.isolate = FALSE)


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1

cellchat.CTR <- cellchat


# Circle plot
netVisual_individual(cellchat.CTR, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# P3TMZ tumor and TME
### change gene names to standard names with "-" instead of with "."
VlnPlot(P3TMZ.full, features = "H2.Aa")
which(P3TMZ.full@assays$RNA@data@Dimnames[[1]] == "H2.Aa")
P3TMZ.full@assays$RNA@data@Dimnames[[1]] <- sub('\\.','-', P3TMZ.full@assays$RNA@data@Dimnames[[1]])
VlnPlot(P3TMZ.full, features = "H2-Aa")
which(P3TMZ.full@assays$RNA@data@Dimnames[[1]] == "H2-Aa")

#Extract the CellChat input files from a Seurat object

data.input <- GetAssayData(P3TMZ.full, assay = "RNA", slot = "data") # normalized data matrix
Idents(P3TMZ.full) <- "celltypes"
labels <- Idents(P3TMZ.full)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels


#Create a CellChat object using data matrix as input

cellchat.obj <- createCellChat(object = P3TMZ.full, meta = meta, group.by = "group")

groupSize <- as.numeric(table(cellchat.obj@idents)) # number of cells in each cell group


#Set the ligand-receptor interaction database        
CellChatDB  <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat.obj@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat.obj)


cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.mouse)

#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 0)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
#Here we also control the parameter edge.weight.max so that we can compare edge weights between different networks.

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram

pathways.show <- c("CSF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

netAnalysis_contribution(cellchat, signaling = pathways.show)


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:9), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 9, targets.use = c(1:8), remove.isolate = FALSE)


#Compute and visualize the network centrality scores

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1


cellchat.TMZ <- cellchat


# merge the 2 cell chat objects for comparison

#Lift up CellChat object and merge together

# Define the cell labels to lift up

group.new = levels(cellchat.CTR@idents)
cellchat.TMZ <- liftCellChat(cellchat.TMZ, group.new)


cellchat.CTR <- netAnalysis_computeCentrality(cellchat.CTR)
cellchat.TMZ <- netAnalysis_computeCentrality(cellchat.TMZ)

object.list <- list(CTR = cellchat.CTR, TMZ = cellchat.TMZ)


cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)


#### Part I: Predict general principles of cell-cell communication ###

# Compare the total number of interactions and interaction strength  

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


#The differential number of interactions or interaction strength in the cell-cell communication network between two datasets
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)


gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Myeloid")

gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Tumor")

patchwork::wrap_plots(plots = list(gg1,gg2))

#Part II: Identify the conserved and context-specific signaling pathways

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")

cellchat <- netEmbedding(cellchat, type = "functional")

cellchat <- netClustering(cellchat, type = "functional")

netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

#Compare the overall information flow of each signaling pathway

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 9, height = 13)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 9, height = 13)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 9, height = 13, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 9, height = 13, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 9, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 9, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs


levels(cellchat@idents$joint)
#[1] "Endothelial"      "Myeloid"          "Ependymal"        "Cycling"          "Pericytes"        "Astrocytes"      
#[7] "Lymphocytes"      "Oligodendrocytes" "OPCs"             "Tumor"    


#The increased signaling means these signaling have higher communication probability (strength) in one dataset compared to the other dataset.
##myeloid vs others
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3,4,5,6,7,8,9,10),  comparison = c(1, 2), angle.x = 90)

##tumor vs others
netVisual_bubble(cellchat, sources.use = 10, targets.use = c(1,2,3,4,5,6,7,8,9),  comparison = c(1, 2), angle.x = 90)

##myeloid vs others
gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use =  c(1,3,4,5,6,7,8,9,10),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in TMZ", angle.x = 90, remove.isolate = F)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 2, targets.use =  c(1,3,4,5,6,7,8,9,10),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in TMZ", angle.x = 90, remove.isolate = F)
#> Comparing communications on a merged object
gg1 + gg2

##tumor vs others
gg1 <- netVisual_bubble(cellchat, sources.use = 10, targets.use =  c(1,2,3,4,5,6,7,8,9),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in TMZ", angle.x = 90, remove.isolate = F)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 10, targets.use =  c(1,2,3,4,5,6,7,8,9),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in TMZ", angle.x = 90, remove.isolate = F)
#> Comparing communications on a merged object
gg1 + gg2

#Identify dysfunctional signaling by using differential expression analysis

pos.dataset = "TMZ"

# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)

# Use the joint cell labels from the merged CellChat object
net <- netMappingDEG(cellchat, features.name = features.name)

# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "TMZ",ligand.logFC = 0.2, receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "CTR",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


##myeloid vs others
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 2, targets.use =  c(1,3,4,5,6,7,8,9,10), comparison = c(1, 2),  angle.x = 90, remove.isolate = F, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 2, targets.use =  c(1,3,4,5,6,7,8,9,10), comparison = c(1, 2),  angle.x = 90, remove.isolate = F, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#> Comparing communications on a merged object
gg1 + gg2


##tumor vs others
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 10, targets.use =  c(1,2,3,4,5,6,7,8,9), comparison = c(1, 2),  angle.x = 90, remove.isolate = F,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 10, targets.use =  c(1,2,3,4,5,6,7,8,9), comparison = c(1, 2),  angle.x = 90, remove.isolate = F,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#> Comparing communications on a merged object
gg1 + gg2

#Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 10, targets.use =  c(1,2,4,3,5,6,7,8,9), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 10, targets.use =  c(1,2,4,3,5,6,7,8,9), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


pathways.show <- c("GAS")

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}


#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

par(mfrow = c(1,2), xpd=TRUE)

ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netAnalysis_contribution(object.list[[i]], signaling = pathways.show, title = paste("Contribution of L-R pair to",pathways.show, "signaling -",names(object.list)[i]))
}
ht[[1]] + ht[[2]]

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}



