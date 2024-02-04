
####  Figure 2  #####

library(Seurat)

#Fig 2A - generating list of differentially expressed genes for each cell type in PDOX vs Naive nude 

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


# generating list of differentially expressed genes for each cell type in GL261 Vs Naive black6 per cell type

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



#Fig 2C - DotPlot of key marker genes for the identified cell types

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



#Fig 2D - DotPlot of key marker genes from the human GBM dataset

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




