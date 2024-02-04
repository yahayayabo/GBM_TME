
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


