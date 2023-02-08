

#### Figure 1  #####

### check the pre-processing file for details on the different steps during pre-processing here https://github.com/yahayayabo/GBM_TME/blob/master/Pre-processing

library(Seurat)
library(harmony)
library(ggplot2)
library(tidyverse)



#1b - generating umap of merged TME cells PDOX and GL261

##### merge PDOX and GL261 data together in one seurat object ######

PDOX_GL261_TME <- merge(MGBM11, y= c(MGBM13, MGBM15, MGBM17, Normal, P8, P13,T16,T101,T192,T233,T347,T470,P3CTR,P3TMZ ), project = "GL_PDOX")


PDOX_GL261_TME <- NormalizeData(PDOX_GL261_TME, normalization.method = "LogNormalize", scale.factor = 10000)
PDOX_GL261_TME <- FindVariableFeatures(PDOX_GL261_TME, selection.method = "vst", nfeatures = 2000)
PDOX_GL261_TME <- ScaleData(PDOX_GL261_TME, verbose = FALSE)
PDOX_GL261_TME <- RunPCA(PDOX_GL261_TME, features = VariableFeatures(object = PDOX_GL261_TME))

# remove batch effect using harmony 

VlnPlot(object = PDOX_GL261_TME, features = "PC_1", group.by = "sample", pt.size = .1)
options(repr.plot.height = 2.5, repr.plot.width = 6)
PDOX_GL261_TME <- PDOX_GL261_TME %>% 
  RunHarmony("sample", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(PDOX_GL261_TME, 'harmony')
harmony_embeddings[1:5, 1:5]


# run umap and tSNE with harmony embeddings 

#group = Naive Nude, PDOX, Naive Black6, Gl261
#sample = each individual PDOX and GL261 plus controls


PDOX_GL261_TME <- PDOX_GL261_TME %>% 
  RunUMAP(reduction = "harmony", dims = 1:15, metric = 'euclidean') %>% 
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(PDOX_GL261_TME, reduction = "umap", label = TRUE, pt.size = .1)
DimPlot(PDOX_GL261_TME, reduction = "umap", group.by = "sample")

# to split by sample
DimPlot(PDOX_GL261_TME, reduction = "umap", group.by = "sample", pt.size = .1, split.by = 'sample')

# generate list of marker genes for each cluster
PDOX_GL261_TME_markers <- FindAllMarkers(PDOX_GL261_TME, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
write.csv(PDOX_GL261_TME_markers, file = "clusters_markers.csv")

# load the merged full TME seurat object (full_TME_sobj)

DimPlot(PDOX_GL261_TME, reduction = "umap", pt.size = .1, label.size = 5,
    cols = c( 'Astrocytes' = '#bd0026', 'Myeloid' = '#006d2c', 'Endothelial' = '#bf812d',
              'OPCs' = '#a6cee3','Oligodendrocytes' = '#1f78b4','Cycling' = '#1a1a1a',
              'Ependymal' = '#fdcc0d','Pericytes' = '#e3f376', 'Lymphocytes' = '#7aff33')) +
    tiff("Fig1b.tiff", units="in", width=7.8, height=5.5, res=400) # to generate high resolution image
dev.off()

# stack chart of clusters per sample

Celltype <- rep(c("Cycling", "Pericytes","Ependymal", "OPCs", 
                  "Oligodendrocytes", "Endothelial", "Astrocytes",  
                  "Lymphocytes", "Myeloid"), 14)


cellnumber <- c(  28,  2,  6,  20,  88,  52,  57,   17, 108,                         
                  24,  5,  1,  26,  98,  55, 113,   16, 142,                                                  
                  50,  9, 19,  76, 250, 129, 266,  221, 614,                                            
                  14, 40, 101, 24, 182, 374, 707,    2, 528,                                         
                  27,  8,  17, 58,  27,  92, 721,    8, 734,                                           
                  97,155,  64, 69,  29, 823, 706,   24, 789,                                
                  55,  8,  25,  3,   1, 265, 113,    2, 237,
                  38,  9,  52,352, 111, 332, 672,    9, 618,
                  58, 16,  12, 49,  13, 307, 453,   43, 363,
                  65, 66,  44,182,  48, 697, 408,    5, 410,
                  52, 19,  12, 24,  11, 331, 138,   23,1134,                               
                  17, 19,  41, 27,   4, 145,  70,    1, 548,
                  110,15,  41, 57,  26, 242,1095,    4, 527,                                              
                  107,10,   6,162,  58, 276, 891,    1, 225)


data = data.frame(samples, cellnumber, Celltype)
data 

g <- ggplot(data, aes(y=cellnumber, x= samples, fill= factor(Celltype, levels= c("Cycling", "Pericytes","Ependymal", "OPCs", 
                                                                                 "Oligodendrocytes", "Endothelial", "Astrocytes",  
                                                                                 "Lymphocytes", "Myeloid"))))


g + geom_bar(position = "fill", stat="identity") + 
  scale_fill_manual(values =  c( 'Cycling' = '#1a1a1a',
                                 'Pericytes' = '#e3f376',
                                 'Ependymal' = '#fdcc0d',
                                 'OPCs' = '#a6cee3',
                                 'Oligodendrocytes' = '#1f78b4',
                                 'Endothelial' = '#bf812d',
                                 'Astrocytes' = '#bd0026',
                                 'Lymphocytes' = '#7aff33',
                                 'Myeloid' = '#006d2c' )) +
  tiff("stack_all.tiff", units="in", width=12.2, height=5.2, res=400)
dev.off()



#1c - generating list of differentially expressed genes PDOX vs Naive nude 

PDOX_vs_NaiveNude_Endothelial <- FindMarkers(PDOX_GL261_TME, ident.1 = "PDOX", ident.2 = "NaiveNude",
                                             group.by = "group", subset.ident = "Endothelial")

PDOX_vs_NaiveNude_Myeloid <- FindMarkers(PDOX_GL261_TME, ident.1 = "PDOX", ident.2 = "NaiveNude",
                                         group.by = "group", subset.ident = "Myeloid")

PDOX_vs_NaiveNude_Astrocytes <- FindMarkers(PDOX_GL261_TME, ident.1 = "PDOX", ident.2 = "NaiveNude",
                                            group.by = "group", subset.ident = "Astrocytes")

PDOX_vs_NaiveNude_OPCs <- FindMarkers(PDOX_GL261_TME, ident.1 = "PDOX", ident.2 = "NaiveNude",
                                      group.by = "group", subset.ident = "OPCs")

PDOX_vs_NaiveNude_Oligodendrocytes <- FindMarkers(PDOX_GL261_TME, ident.1 = "PDOX", ident.2 = "NaiveNude",
                                                  group.by = "group", subset.ident = "Oligodendrocytes")

PDOX_vs_NaiveNude_Pericytes <- FindMarkers(PDOX_GL261_TME, ident.1 = "PDOX", ident.2 = "NaiveNude",
                                           group.by = "group", subset.ident = "Pericytes")

PDOX_vs_NaiveNude_Ependymal <- FindMarkers(PDOX_GL261_TME, ident.1 = "PDOX", ident.2 = "NaiveNude",
                                           group.by = "group", subset.ident = "Ependymal")


# generating list of differentially expressed genes  GL261 Vs Naive black6 per cell type

GL261_vs_Black6_Endothelial <- FindMarkers(PDOX_GL261_TME, ident.1 = "GL261", ident.2 = "Black6",
                                             group.by = "group", subset.ident = "Endothelial")

GL261_vs_NaiveNude_Myeloid <- FindMarkers(PDOX_GL261_TME, ident.1 = "GL261", ident.2 = "Black6",
                                         group.by = "group", subset.ident = "Myeloid")

GL261_vs_NaiveNude_Astrocytes <- FindMarkers(PDOX_GL261_TME, ident.1 = "GL261", ident.2 = "Black6",
                                            group.by = "group", subset.ident = "Astrocytes")

GL261_vs_NaiveNude_OPCs <- FindMarkers(PDOX_GL261_TME, ident.1 = "GL261", ident.2 = "Black6",
                                      group.by = "group", subset.ident = "OPCs")

GL261_vs_NaiveNude_Oligodendrocytes <- FindMarkers(PDOX_GL261_TME, ident.1 = "GL261", ident.2 = "Black6",
                                                  group.by = "group", subset.ident = "Oligodendrocytes")

GL261_vs_NaiveNude_Pericytes <- FindMarkers(PDOX_GL261_TME, ident.1 = "GL261", ident.2 = "Black6",
                                           group.by = "group", subset.ident = "Pericytes")

GL261_vs_NaiveNude_Ependymal <- FindMarkers(PDOX_GL261_TME, ident.1 = "GL261", ident.2 = "Black6",
                                           group.by = "group", subset.ident = "Ependymal")

write.csv(All_list, "Supplementary_Table_S3.csv", quote = FALSE,  row.names = TRUE)



#1e


#### DotPlot of key genes in different cell types

# subset each cell type from the PDOX_GL261_TME

# Myeloid
Myeloid <- subset(PDOX_GL261_TME, subset = Celltype == "Myeloid")
genes = c("Spp1", "Cst7", "Ch25h", "P2ry12", "Gpr34", "Tmem119")
DotPlot(Myeloid, features = genes, dot.min = 0, scale = TRUE, dot.scale = 8, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() 

#OPC 
OPCs <- subset(PDOX_GL261_TME, subset = Celltype == "OPCs")
genes = c("S100a1","Cspg5","Cacng4", "Pdgfra","Ccnd2", "Stmn1","Apod", "Plp1")
DotPlot(OPCs, features = genes, dot.min = 0, scale = TRUE, dot.scale = 8, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() 

#Astrocytes 
Astrocytes <- subset(PDOX_GL261_TME, subset = Celltype == "Astrocytes")
genes = c("H2.D1", "Gfap", "Vim","Slc1a2", "Slc1a3", "Slc38a3")
DotPlot(Astrocytes, features = genes, dot.min = 0, scale = TRUE, dot.scale = 8, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() 

#Endothelial 
Endothelial <- subset(PDOX_GL261_TME, subset = Celltype == "Endothelial")
genes = c("Fn1", "Angpt2","Ifitm1", "Spock2", "Igfbp5", "Ttyh1")
DotPlot(Endothelial, features = genes, dot.min = 0, scale = TRUE, dot.scale = 8, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() 


##1f


## Upload Darmanis seurat object (Darmanis_obj)

Darmanis_obj <- readRDS("~/Darmanis_obj.rds")

dim(Darmanis_obj)
#[1] 23368  3589

## see expression of activation markers per celltype

Idents(Darmanis_obj) <- "Celltype"
Darmanis_obj@active.ident <- factor(Darmanis_obj@active.ident, levels=c("Immune cell", "OPC", "Astocyte", "Vascular", "Oligodendrocyte", "Neuron", "Neoplastic"))

var=c( "ANGPT2","IFITM1", "FN1", "CRYAB", "VIM", "GFAP","CSPG5","CACNG4", "PDGFRA","CLEC7A" , "APOE","SPP1")

DotPlot(Darmanis_obj, features = var,   scale = FALSE, dot.scale = 8, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() 


