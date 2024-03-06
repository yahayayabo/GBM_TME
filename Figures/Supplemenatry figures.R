
#Supplementary Figures

library(ggplot2)
library(patchwork)
library(Seurat)
library(dittoSeq)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(viridis)
library(egg)
library(grid)
library(lattice)
library(gtools)
library(Biobase)
library(RColorBrewer)
library(stringr)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(ggmin)

### Supplementary figure 2

#load data

PDOX_GL261_TME <- readRDS("pathtoyourfolder/PDOX_GL261_TME.rds")

#Fig.S2A
features <- c("Bcan","Slc1a2",  "Gja1",  #Astro
               "P2ry12", "Cx3cr1", "Lyz2",  #Myeloid
               "Cxcl12", "Ptprb", "Cldn5",   #Endothelial
               "Cacng4", "Pdgfra", "Cspg5", #OPC
               "Plp1",   "Mbp",    "Mag",    #Oligo
               "Top2a",  "Cenpf",  "Mki67", #Cycling
               "Ttr",    "Enpp2",  "Igf2",  #Ependymal
               "Ccdc153","Tmem212", "Enkur",   #Pericytes
               "Trac",   "Ccl5",  "Nkg7" ) #Lymphocytes


## stacked violin plot
VlnPlot(PDOX_GL261_TME,features=features, split.by="merged_celltypes", pt.size = 0, stack=T, flip=T, 
        cols = c('Astrocytes' = '#bd0026','Myeloid' = '#006d2c', 'Endothelial' = '#bf812d', 'OPCs' = '#a6cee3', 
                 'Oligodendrocytes' = '#1f78b4', 'Cycling' = '#1a1a1a', 'Ependymal' = '#fdcc0d', 
                 'Pericytes' = '#e3f376','Lymphocytes' = '#7aff33' ))




#Fig.S2B
## lymphocytes
object <-  lympho
subtype_name <- "Lymphocytes"

# some selected higly expressed genes across lymphocyte clusters
genes <- c("Cd8b1","Cd8a","Trbc2","Cd2", "Cd3d","Cd3e","Cd3g", "Ctla4", "Gimap3","Foxp3", "Tigit", "Lag3", "Pdcd1","Ccl5","Cytip","Trac", #T cells
           "Cd19", "H2.Aa","H2.Ab1","Bank1","Ly6d","Fcmr", "Klf2","Cd37","Ms4a1", "Cd72","H2.Eb1", "Ece1","Ran",  # B cell
           "Klrb1b", "Klrb1c","Klra8","Klra4", "Klra13.ps", "Nkg7", "Gzma", "Gzmb", "Klrd1", "Gzmc","Eomes","Fcer1g","Itgb3","Klrg9", "Klrg1","Klrg4")  #  NK cells




# HEATMAP ----
DefaultAssay(object) <- "RNA" # you will need to use SCT if the RNA is not normalised and scaled
Idents(object) <- "groups"
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



heatmap.markers <- ggplot(data = expr.long) +
  geom_tile(aes(x=factor(cell.type, level=c("RNA.Naive_Nude", "RNA.PDOX", "RNA.Naive_Black6","RNA.GL261")), y = reorder(Gene,-GeneOrder),fill = z_score), color = "white", size = 0.1) + # bigger size means bigger spacer
  scale_fill_gradient2(low = "#3182bd",mid= "white", high = "darkred") +
  scale_y_discrete(position = "right") +
  xlab("") + ylab("") +
  theme_grey(base_size = 13) + # text size
  ggtitle(subtype_name) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10, colour = "gray50"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

heatmap.markers

x <- 4 # you can change this to shape the figure
y <- 5 # you can change this to shape the figure
ggsave("lymphocytes.pdf", plot = heatmap.markers, width = x, height = y) 


#Fig.S2C 

#load Suva gene list
suva_list <- read.table("Suva_paper_subtypes.csv",header=TRUE,sep=";")
head(suva_list)

#Set column headers as id variables
suva_list <- melt(suva_list)

#create separate vectors containing the Suva signature genes
#MES
MES_genes <- as.vector(suva_list$MES[suva_list$MES!=""])

#AC
AC_genes <- as.vector(suva_list$AC[suva_list$AC!=""])

#OPC
OPC_genes <- as.vector(suva_list$OPC[suva_list$OPC!=""])

#NPC
NPC_genes <- as.vector(suva_list$NPC[suva_list$NPC!=""])

#get only the unique subtype genes that are present in the normalized and filtered dataset
MES.subtype.related.genes<-unique(MES_genes[MES_genes%in%rownames(P3CTR_gem)])

AC.subtype.related.genes<-unique(AC_genes[AC_genes%in%rownames(P3CTR_gem)])

NPC.subtype.related.genes<-unique(NPC_genes[NPC_genes%in%rownames(P3CTR_gem)])

OPC.subtype.related.genes<-unique(OPC_genes[OPC_genes%in%rownames(P3CTR_gem)])

#create a list containing all subtype related genes
vector <- list(MES.subtype.related.genes, AC.subtype.related.genes, NPC.subtype.related.genes, OPC.subtype.related.genes)

#check class (supposed to be a list; package requirement)
class(vector)

Suva_Scores <- scrabble::score(P3CTR_gem, vector, center = F, nbin = 30, n = 100, replace = F)
#check class (it is matrix)
class(Suva_Scores)

#assign "MES", "AC", "NPC", "OPC" as colnames
colnames(Suva_Scores) <- c("MES", "AC", "NPC", "OPC")


#First find the cells that have 1 subtype highly expressed
rows <- apply(Suva_Scores , 1 , function(x) any( x > 1 ))

#Subset on this vector

criteria_one <- apply(Suva_Scores,1,function(x) colnames(Suva_Scores)[which.max(x)])

table(criteria_one)



AC_sub <- Suva_Scores[criteria_one=="AC",]
MES_sub <- Suva_Scores[criteria_one=="MES",]
NPC_sub <- Suva_Scores[criteria_one=="NPC",]
OPC_sub <- Suva_Scores[criteria_one=="OPC",]

AC_Quant <- quantile(AC_sub[,"AC"], 0.9)
MES_Quant <- quantile(MES_sub[,"MES"], 0.9)
NPC_Quant <- quantile(NPC_sub[,"NPC"], 0.9)
OPC_Quant <- quantile(OPC_sub[,"OPC"], 0.9)

Quantiles_scores <- c(AC_Quant,MES_Quant,NPC_Quant,OPC_Quant)
names(Quantiles_scores) <- c("AC","MES","NPC","OPC")
criteria_two <- c()

for (cell in rownames(Suva_Scores)){
  order_scores <- Suva_Scores[cell,][order(Suva_Scores[cell,])]
  sec_pos_subtype <- names(order_scores[3])
  second <- Suva_Scores[,sec_pos_subtype] > Quantiles_scores[sec_pos_subtype]
  temp_cell <- Suva_Scores[cell,][order(Suva_Scores[cell,])]
  if(temp_cell[3]-temp_cell[2] > 0.3){
    criteria_two[cell] <- paste(criteria_one[cell],sec_pos_subtype,sep="/")
  }else{
    criteria_two[cell] <- criteria_one[cell]
  }
}

criteria_two




#Fig.S2D 
#feature plot
FeaturePlot(PDOX_GL261_TME, features = c("Mki67","Top2a"))

#Fig.S2E and S2F

library(Seurat)
library(UCell)
library(R.filesets)
library(dplyr)
library(ggplot2)

srt_final <- loadRDS("~/path/srt_final.rds")
marker.genes <- as.data.frame(readxl::read_xlsx("~/path/Yahaya_Marker_Genes_From_Paper_June2023.xlsx"))
colnames(marker.genes) <- c("Astrocytes", "Endothelial", "Ependymal", "Myeloid", "Oligodendrocytes", "OPCs", "Pericytes")

## extracting cycling cells as a separate seurat obj
cycling.object = srt_final[, srt_final@meta.data$merged_celltypes == "Cycling"]


cell.types <- lapply(colnames(marker.genes), function(x){
  
  tmp <- is.na(marker.genes[, x])
  tmp <- marker.genes[, x][!tmp]
  return(tmp)

})
names(cell.types) <- colnames(marker.genes)


cycling.object <- AddModuleScore_UCell(cycling.object, 
                                       features = cell.types, 
                                       assay = "RNA")

## getting UCELL cell type columns to a separate df
ucell.cell.types.df <- cycling.object@meta.data[,13:19]
colnames(ucell.cell.types.df) <- colnames(marker.genes)

## without scaling
pheatmap::pheatmap(ucell.cell.types.df, 
                   angle_col = "45", 
                   show_rownames = F, 
                   main = "Cycling cluster cells scored by other cell type signatures",
                   scale = "none")


## scores are scaled
pheatmap::pheatmap(ucell.cell.types.df, 
                   angle_col = "45", 
                   show_rownames = F, 
                   main = "Cycling cluster cells scored by other cell type signatures",
                   scale = "row",
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100))

## assinging the cell type with max UCELL score to each single cell
new.cell.types.assigned <- colnames(ucell.cell.types.df)[apply(ucell.cell.types.df, 1, function(x) which.max(x))]

##
cycling.object@meta.data$new.celtype <- new.cell.types.assigned
cycling.object@meta.data$new.celtype <- factor(cycling.object@meta.data$new.celtype, 
                                               levels = c("Astrocytes", 
                                                          "Myeloid", 
                                                          "Endothelial", 
                                                          "OPCs", 
                                                          "Oligodendrocytes", 
                                                          "Ependymal", 
                                                          "Pericytes"))

srt_final@meta.data$new.celltype <- srt_final@meta.data$merged_celltypes
srt_final@meta.data$new.celltype[srt_final@meta.data$merged_celltypes == "Cycling"] <- as.character(cycling.object@meta.data$new.celtype)


Idents(cycling.object) <- cycling.object@meta.data$new.celtype


tmp.2 <- tibble(
  cell_type = Idents(cycling.object)
  ) %>%
  group_by(cell_type) %>%
  count() %>%
  mutate(
    percent=(100*n)/sum(n)
  )

tmp.2$percent <- round(tmp.2$n/sum(tmp.2$n)*100,1)
tmp.2 <- tmp.2[order(-tmp.2$percent),]


xlsx::write.xlsx(as.data.frame(tmp.2), 
                 "~/myData/projects/Yahaya_project/For_Yahayas_Paper/Yahaya_Cycling_Subcelltypes_numbers_and_percentages_June2023.xlsx")


my_cols <- c('3'='#F68282','15'='#31C53F','5'='#1FA195','1'='#B95FBB','13'='#D4D915',
             '14'='#28CECA','9'='#ff9a36','8'='#2FF18B','11'='#aeadb3','6'='#faf4cf',
             '2'='#CCB1F1','12'='#25aff5','7'='#A4DFF2','4'='#4B4BF7','16'='#AC8F14',
             '10'='#E6C122')

my_cols.all <- c('Astrocytes'='#F68282',
             'Myeloid'='#E6C122', 
             'Endothelial'='#4B4BF7',
             'OPCs'='#A4DFF2',
             'Oligodendrocytes'='#CCB1F1', 
             'Cycling'='#25aff5', 
             'Ependymal'='#31C53F',
             'Pericytes'='#aeadb3',
             'Lymphocytes'='#D4D915')

my_cols.cycling <- c('Astrocytes'='#F68282',
             'Myeloid'='#E6C122', 
             'Endothelial'='#4B4BF7',
             'OPCs'='#A4DFF2',
             'Oligodendrocytes'='#CCB1F1', 
             'Cycling'='#25aff5', 
             'Ependymal'='#31C53F',
             'Pericytes'='#aeadb3')

yahaya.cols = c('Astrocytes' = '#bd0026', 
                'Myeloid' = '#006d2c', 
                'Endothelial' = '#bf812d',
                'OPCs' = '#a6cee3',
                'Oligodendrocytes' = '#1f78b4',
                'Cycling' = '#1a1a1a',
                'Ependymal' = '#fdcc0d',
                'Pericytes' = '#e3f376', 
                'Lymphocytes' = '#7aff33')

yahaya.cols.cycling = c('Astrocytes' = '#bd0026', 
                'Myeloid' = '#006d2c', 
                'Endothelial' = '#bf812d',
                'OPCs' = '#a6cee3',
                'Oligodendrocytes' = '#1f78b4',
                'Cycling' = '#1a1a1a',
                'Ependymal' = '#fdcc0d',
                'Pericytes' = '#e3f376')

Idents(srt_final) <- srt_final@meta.data$merged_celltypes
p1 <- DimPlot(srt_final, cols = yahaya.cols, order = T) + 
        ggtitle("Original UMAP") +
        theme(plot.title = element_text(hjust = 0.5, size = 18))
p2 <- DimPlot(cycling.object, cols = yahaya.cols.cycling, order = tmp.2, pt.size = 1.2) + 
        ggtitle("Subclustering of Cycling cells") + 
        theme(plot.title = element_text(hjust = 0.5, size = 18))

Idents(srt_final) <- srt_final@meta.data$new.celltype
p3 <- DimPlot(srt_final, cols = yahaya.cols, group.by = "new.celltype", order = T) + 
        ggtitle("Updated UMAP") +
        theme(plot.title = element_text(hjust = 0.5, size = 18))





p4 <- ggplot(tmp.2, 
             aes(x="",
                 y=percent, 
                 fill=cell_type)) +
        geom_bar(width=1, 
                 stat = "identity") +
        coord_polar("y") +
        scale_fill_manual(values = yahaya.cols.cycling) +  
        geom_text(aes(label = percent),
                  size = 9,
                  position = position_stack(vjust = 0.5)) +
        theme_classic() +
        theme(
              axis.text.x=element_blank(),
              plot.title = element_text(vjust = 0.5, 
                                        size = 18, 
                                        face = "bold"),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              text = element_text(size = 20),
              legend.title = element_text(colour = "white")) +
        labs(fill = "cell_type", 
             x = NULL, 
             y = NULL,
             title = "Cellular subtypes of Cycling cluster")
        


p5 <- patchwork::wrap_plots(list(p1, p3, p2, p4), ncol = 2)
ggsave("/path/image.png", 
       width = 20, 
       height = 16, 
       dpi = 150, 
       plot = p5)




### Supplementary figure 3

#Fig. S3A

## upset plot

#upload DEG lists
#make a list of DEGS per comparison
library(readxl)
df=read_excel("path/DEGS_all_list.xlsx")
df=lapply(df, function(x)x[!is.na(x)])

PDOX_Nu_NB_Astro <- df$PDOX_Nude_Astro
GL261_BL6_NB_Astro <- df$GL261_BL6_Astro
PDOX_Nu_NB_ECs <- df$PDOX_Nude_Endo
GL261_BL6_NB_ECs <- df$GL261_BL6_Endo
PDOX_Nu_NB_Myeloid <- df$PDOX_Nude_Myeloid
GL261_BL6_NB_Myeloid <- df$GL261_BL6_Myeloid
PDOX_Nu_NB_OPCs <- df$PDOX_Nude_OPCs
GL261_BL6_NB_OPCs <- df$GL261_BL6_OPCs


Final_list <- list(PDOX_Nu_NB_Astro, GL261_BL6_NB_Astro, PDOX_Nu_NB_ECs, GL261_BL6_NB_ECs, 
                   PDOX_Nu_NB_Myeloid, GL261_BL6_NB_Myeloid, PDOX_Nu_NB_OPCs, GL261_BL6_NB_OPCs)

named_list <- list(PDOX_Nu_NB_Astro=PDOX_Nu_NB_Astro, GL261_BL6_NB_Astro=GL261_BL6_NB_Astro, PDOX_Nu_NB_ECs=PDOX_Nu_NB_ECs, 
                   GL261_BL6_NB_ECs=GL261_BL6_NB_ECs,  PDOX_Nu_NB_Myeloid=PDOX_Nu_NB_Myeloid, GL261_BL6_NB_Myeloid=GL261_BL6_NB_Myeloid, 
                   PDOX_Nu_NB_OPCs=PDOX_Nu_NB_OPCs, GL261_BL6_NB_OPCs=GL261_BL6_NB_OPCs)

library(UpSetR)

require(ggplot2); require(plyr); require(gridExtra); require(grid);

upset(fromList(named_list), nsets = 8, keep.order = TRUE, order.by = "freq", sets.bar.color=c("red","black","red","black","red","black","red","black"), 
      sets = c("GL261_BL6_NB_Astro","PDOX_Nu_NB_Astro", "GL261_BL6_NB_OPCs", "PDOX_Nu_NB_OPCs",
               "GL261_BL6_NB_ECs","PDOX_Nu_NB_ECs","GL261_BL6_NB_Myeloid","PDOX_Nu_NB_Myeloid"))

#download list of upset intersections
str(named_list)
df2 <- data.frame(gene=unique(unlist(named_list)))

head(df2)

dim(df2)
#[1] 1398    1
df1 <- lapply(named_list,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "comparison")

head(df1)


dim(df1)
#1] 2172    2

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(comparison) %>% 
    pull("comparison") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

head(df_int,n=5)


dim(df_int)
#[1] 1398    2

df_int %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))



write.csv(df_int, "Upset_sets.csv", quote = FALSE,  row.names = TRUE)


upset(fromList(named_list),nsets = 10) 


### Supplementary figure 4

#Fig. S4A
### order split umap
Merged_myeloid_data$dataset <- factor(x = Merged_myeloid_data$dataset, levels = c("PDOX", "Bozena", "Yolanda", "Pombo"))
DimPlot(Merged_myeloid_data, group.by = "Cell_Type", split.by = "dataset", cols = c("BAM" = "red" , 'MG' = '#4daf4a', 'Mo/MG' = 'blue'))

#Fig. S4B

FeaturePlot(Merged_myeloid_data, features = c("Ptprc", "Gpr34", "Ccr2", "Lyve1"))

#Fig. S4C

#stack chart per datasets
dittoBarPlot(Merged_myeloid_data,"CellType", group.by = "dataset", 
             var.labels.reorder = c(1,3,2),
             x.reorder = c(2,1,4,3), 
             color.panel = c("BAM" = "red" ,'Mo/MG' = 'blue', 'MG' = '#4daf4a')) +
  theme(text = element_text(size = 16))



### Supplementary figure 5

#Fig. S5A

#load data
mydata <- readRDS("myeloid_object.rds")


genes_selected <-  c("Il3ra", "Tcf4", "Tpm2", #pDCs
                     "Ccr7","Lamp3", "Fscn1", "Ccl22", #migratory DCs
                     "Ccr9", "Cd300c", "Cox6a2", "Clec10a", #cDC2 - 2
                     "Xcr1","Clec9a", "Thbd", "Rab7b", #cDC1
                     "Cd5", "Cd2","Axl", #pre-DCs 
                     "Rps12", "Rps27a", "Rps29", "Rpl17", "Rpl23a", "Gnb2l1", "Fau",
                     "Ptprb", "Cxcl12", "Cldn5", "Ctla2a", #endothelial
                     "Bcan", "Slc1a2", "Slc1a3", "Aldoc", #astrocytes
                     "Fos",  "Jund","Ier5",  "Egr1", #activation markers
                     "Mcm4",  "Cenpe","Cdk1", "Pcna","Mki67", "Top2a", #cc
                     "Tgm2","Tgfb1", "Il1rn", "Il27ra","Clec7a","Klf4","Mrc1", "Il10","Cd163","F13a1","Sepp1", "Tgfb2", "Msr1", "Arg1", "Rnase1", "Il2ra", #M2
                     "Tlr2","Stat1","Ido1", "Il6", "Cd80", "Cxcl10","Cd86","Nfkbia","Il1b","Tlr4","Cd34","Irf5", "Tnf", "Cxcl13","Cd72","Ccl4","Ccl3", #M1
                     "Itgam", "Ptprc","Csf1r", "Hexb") #myeloid 


Idents(mydata) <- "New_Clusters"

mydata@active.ident <- factor(mydata@active.ident, levels=c("CL0", "CL1", "CL2", "CL3", "CL4", "CL5", "CL6", "CL7", "CL8", "Endothelial", "Astrocytes"))

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(mydata, features = genes_selected, scale = FALSE, dot.scale = 8, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() +
  tiff("dotplot.tiff", units="in", width=6.5, height=15, res=400)
dev.off()

library(Cairo)

Cairo::Cairo(
  30, #length
  30, #width
  file = paste("nameofplot", ".png", sep = ""),
  type = "png", #tiff
  bg = "transparent", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
plot(p) #p is your graph object 
dev.off()


#Fig. S5B

## doublet detection using DoubletFinder

library(DoubletFinder)
nExp <- round(ncol(seurat_obj) * 0.04)  # expect 4% doublets
seurat_obj <- doubletFinder_v3(seurat_obj, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

DF.name = colnames(seurat_obj@meta.data)[grepl("DF.classification", colnames(seurat_obj@meta.data))]



cowplot::plot_grid(ncol = 2, DimPlot(seurat_obj, group.by = "New_Clusters") + NoAxes(),
                   DimPlot(seurat_obj, group.by = DF.name) + NoAxes())

DimPlot(seurat_obj, group.by = DF.name) 
VlnPlot(seurat_obj, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)


saveRDS(seurat_obj, "seurat_obj.rds")

data.filt = seurat_obj[, seurat_obj@meta.data[, DF.name] == "Singlet"]

dim(data.filt)

DimPlot(data.filt, group.by = "DF.classifications_0.25_0.09_279")


#Fig. S5D


#Fig. S5E



### Supplementary figure 6

#Fig. S6A
DimPlot(Merged_myeloid_data, group.by = "Bozena_clusters", label = T,
        cols = c(  "BAM"= "#F8766D","Macrophage"= "#A3A500", "inMoM"= "#D89000", "Monocytes"= "#FF62BC", 
                   "MG8"= "#E76BF3", "MG7"= "#9590FF", "MG6"="#00B0F6", "MG3"= "green","MG2"= "#00BF7D", "MG1"= "#39B600"))


dittoBarPlot(Merged_myeloid_data,"Bozena_clusters", group.by = "sample", 
             var.labels.reorder = c(1,3,2,10,9,8,7,6,5,4),
             x.reorder = c(5,8,9,13,14,11,12,7,10,6,4,3,1,2), 
             color.panel = c( "BAM"= "#F8766D","Macrophage"= "#A3A500", "inMoM"= "#D89000", "Monocytes"= "#FF62BC", 
                              "MG8"= "#E76BF3", "MG7"= "#9590FF", "MG6"="#00B0F6", "MG3"= "green","MG2"= "#00BF7D", "MG1"= "#39B600")) +
  theme(text = element_text(size = 16))


#Fig. S6D
seurat_obj <- readRDS("seurat_obj_final.rds")

Idents(seurat_obj) <- "New_Clusters"
seurat_obj@active.ident <- factor(seurat_obj@active.ident, levels=c("CL0", "CL1", "CL2", 'CL3', "CL4", "CL5" , "CL6", "CL7", "CL8"))



#Extract Meta Data
CellInfo <- seurat_obj@meta.data

#Assign colors

ClusterColors=c(CL0 = "#22830B", CL1 = "#2CAA0E", CL2 = "#9EF57B", CL3 = "#99B953", CL4 = "#FFD966", 
                CL5 = "#7F6000", CL6 = "#F08C06", CL7 = "blue", CL8 = "#e31a1c")


#table(seurat_obj$groups)
#GL261 Naive_Black6   Naive_Nude         PDOX 
#864          528          734         4851 

GroupColors=c(GL261 = "darkred", Naive_Black6 = "#fc9272",  Naive_Nude = "black", PDOX = "grey")



markers = FindAllMarkers(seurat_obj,logfc.threshold = 0.5,test.use = "wilcox",only.pos = T )


top20.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)



Idents(seurat_obj) <- "New_Clusters"
seurat_obj@active.ident <- factor(seurat_obj@active.ident, levels=c("CL0", "CL1", "CL2", 'CL3', "CL4", "CL5" , "CL6", "CL7", "CL8"))

seurat_obj.sce=as.SingleCellExperiment(seurat_obj)
plot.data<-as.data.frame(assay(seurat_obj.sce, "logcounts"))
plot.data<-plot.data[top20.markers$gene,]
plot.data <- plot.data - rowMeans(plot.data)
plot.data=na.omit(plot.data)

CellInfoS=seurat_obj@meta.data


column_annot <-CellInfoS[,c("New_Clusters","groups"),drop=F]

group_order <- c("Naive_Black6", "GL261", "Naive_Nude", "PDOX")
cluster_order <- c("CL0", "CL1", "CL2", 'CL3', "CL4", "CL5" , "CL6", "CL7", "CL8")

column_annot$groups <- factor(column_annot$groups, levels = group_order)
column_annot$New_Clusters <- factor(column_annot$New_Clusters, levels = cluster_order)


column_annot$groups = as.factor(as.character(column_annot$groups))
column_annot=with(column_annot, column_annot[order(groups), , drop=F])
column_annot=with(column_annot, column_annot[order(New_Clusters), , drop=F])
plot.data<-plot.data[,row.names(column_annot)]




column.colors=list()
column.colors[["groups"]]<-GroupColors
column.colors[["New_Clusters"]]<-ClusterColors
Patient=as.matrix(column_annot[,c("groups"),drop=F])


Cluster=as.matrix(column_annot[,c("New_Clusters"),drop=F])

colanno <- columnAnnotation(df=column_annot,
                            show_annotation_name =T,show_legend = F,col=column.colors)
top3 <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
genes= top3$gene
rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#cb181d"))



HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)        
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)


lgd1=Legend(labels = levels(as.factor(column_annot$New_Clusters)),title="New_Clusters",legend_gp = gpar(fill=ClusterColors,fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$groups)),title="groups",legend_gp = gpar(fill=GroupColors,fontsize=5))

draw(HM,heatmap_legend_list = list( lgd,lgd1, lgd2), heatmap_legend_side = "right")
dev.off()




#Fig. S6E

mydata <- readRDS("myeloid_seurat_obj.rds")

Idents(mydata) <- "celltype"

MG <- subset(mydata, subset = celltype == "MG")
mg_genes <- c("Apoe","Cxcl13", "Spp1", "Ccl6", "Itgax", "H2.Eb1", "Ccl12", "Hspa8", "Plp1")

Macrophage <- subset(mydata, subset = celltype == "Macrophage")     
Mo_genes <- c("Cst3", "Hexb", "Ifi202b", "Cst7", "Siglech", "Il1b", "Tgfbi", "Thbs1", "Vegfa", "Arg1")


BAM <- subset(mydata, subset = celltype == "BAM")
bam_genes <- c("Cst3", "Jun", "Ifi202b", "Mef2c", "Cd14", "Tgfbi", "Ly6a", "Il1b", "Arg1", "Cxcl10")


# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(mydata, features = genes,  group.by = "groups", scale = TRUE, dot.scale = 8, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() 


Idents(MG) <- "groups"
MG@active.ident <- factor(MG@active.ident, levels=c("Naive_Nude", "PDOX", "Naive_Black6", "GL261"))

DotPlot(MG, features = mg_genes,  scale = TRUE, dot.scale = 8, col.min = -2.0, col.max = 2.0, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() 

Idents(Macrophage) <- "groups"
Macrophage@active.ident <- factor(Macrophage@active.ident, levels=c("Naive_Nude", "PDOX", "Naive_Black6", "GL261"))

DotPlot(Macrophage, features = Mo_genes,  scale = TRUE, dot.scale = 8, col.min = -2.0, col.max = 2.0, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() 

Idents(BAM) <- "groups"
BAM@active.ident <- factor(BAM@active.ident, levels=c("Naive_Nude", "PDOX", "Naive_Black6", "GL261"))

DotPlot(BAM, features = bam_genes, group.by = "groups", scale = TRUE, dot.scale = 8, col.min = -2.0, col.max = 2.0, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() 

dittoDotPlot(BAM, vars = bam_genes, group.by = "groups", min.color = "blue",  max.color = "red", 
             size = 10, min.percent = 0, max.percent = 1, y.reorder = c(3, 4, 2, 1)) +
  RotatedAxis()+coord_flip()



dittoDotPlot(Myeloid_subset_enrichIt, vars = genes_selected,  group.by = "groups", scale = FALSE, min.color = "blue",  max.color = "red", 
             size = 8,  min.percent = 0,max.percent = 1) +   RotatedAxis() + coord_flip() 

### Supplementary figure 7

#Fig. S7A

#heatmap of functional markers expression

genes <-  c( "Cd47","Sirpa", #do not eat me signal
             "Lag3","Pdcd1","Cd274","Havcr2", "Vsir", #check point inhibitors 
             "H2.Eb1", "H2.Ab1", "H2.Aa", "Cd86", "Igf1","Itgax", #antigen presentation
             "Tyrobp","Trem2", "Axl", #phargocytosis
             "Cxcl16","Clec7a","Cd52","Cd74", #sensome
             "Fn1","Ccl3", "Cxcl13") #migration


DoHeatmap(Myeloid_subset, features = genes, group.by = "New_Clusters",disp.min = -2.0, disp.max = 2.0, slot = "scale.data") +
  scale_fill_gradientn(colours = c("#43a2ca", "white", "#e31a1c"))

Idents(Myeloid_subset) <- "New_Clusters"

Myeloid_subset@active.ident <- factor(Myeloid_subset@active.ident, levels=c("CL0", "CL1", "CL2", "CL3", "CL4", "CL5", "CL6", "CL7", "CL8"))

DotPlot(Myeloid_subset, features = genes, scale = FALSE, dot.scale = 8, cols = c("blue", "red")) +
  RotatedAxis() + coord_flip() + 
  tiff("dotplot_functional.tiff", units="in", width=6.5, height=7, res=400)
dev.off()


#Fig. S7B

##convert seurat obj metadata to data frame
mydata.df <- data.frame(mydata[[]], Idents(mydata))
head(mydata.df)
mydata.df.subset <- mydata.df[, c(10,11,12,13,14,15,19)]
rownames(mydata.df.subset) <- NULL

library(corrplot)
library(Hmisc)

mydata.df.subset.cor = cor(mydata.df.subset) #default =pearson

mydata.cor = cor(mydata.df.subset, method = c("spearman"))

mydata.rcorr = rcorr(as.matrix(mydata.df.subset))

mydata.rcorr

mydata.coeff = mydata.rcorr$r
mydata.p = mydata.rcorr$P
write.csv(mydata.p, file = "functional_scores_corr.p.csv")

corrplot(mydata.cor)
corrplot(mydata.cor, order = 'AOE') # after 'AOE' reorder

paletteFunc <- colorRampPalette(c(low = "#3182bd",mid= "white", high = "darkred"));
palette     <- paletteFunc(10)
corrplot(mydata.cor, method = 'circle', order = 'AOE', type = 'lower', diag = TRUE, tl.pos="ld", tl.cex = 0.6,
         tl.col="black",col=palette)

corrplot(mydata.cor, method = 'ellipse', order = 'AOE', type = 'upper')

corrplot(mydata.cor, order = 'hclust', addrect = 2)

## add significant level stars
corrplot(mydata.cor, p.mat = mydata.p, method = 'color', diag = FALSE, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')



### Supplementary figure 8

#Fig. S8A

human_myeloid <- readRDS("human_myeloid.rds")

###Add signature scores
library(readxl)
Human_Mg <- read_excel("Human_MG.xlsx")
Human_Mo <- read_excel("Human_Mo.xlsx")

Human_mg = list(c("HEXB", "P2RY12", "P2RY13", "TMEM119", "GPR34", "CXCR1", "TGFBR1", "CST3", "OLFML3", "SALL1", "SELPLG", "TREM2"))
Human_Mo = Human_Mo$Macrophages


human_myeloid  <- AddModuleScore(human_myeloid , features = Human_mg, name = "Human_MG_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(Human_Mo), name = "Human_Mo_score", search = T)



DimPlot(human_myeloid, reduction = "umap", group.by = "Yabo_MG_Mo_signature_scores", cols = c('Mg.TAM' = '#4daf4a', 'Mo.TAM' = 'blue'))

FeaturePlot(human_myeloid, features = c("TMEM119","LYZ", "MRC1", "CX3CR1", "S100A9", "LYVE1"))

#Fig. S8B
dittoBarPlot(human_myeloid,"YaboScores", group.by = "sample", 
             var.labels.reorder = c(2, 1),
             x.reorder = c(12,29,8,27,28,23,25,11,9,30,20,6,22,24,35,26,1,7,2,18,
                           33,31,4,5,19,3,10,34,32,16,15,17,14,21,13), 
             color.panel = c('Mg.TAM' = '#4daf4a', 'Mo.TAM' = 'blue')) +
  theme(text = element_text(size = 16))



#show siganture scores on umap

#first add sinature scores

Phago <- read_excel("Human_Phagocytosis.xlsx")
APC <- read_excel("Human_APC.xlsx")
Sensome<-read_excel("Human_sensome.xlsx")
migration <- read_excel("Human_migration.xlsx")

Phago_expressed <- Phago$Phagocytosis
APC_expressed <- APC$APC
Sensome_expressed <- Sensome$Sensome
migration_expressed <- migration$Migration


human_myeloid  <- AddModuleScore(human_myeloid , list(Phago_expressed), name = "Phagocytosis_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(APC_expressed), name = "APCs_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(migration_expressed), name = "migration_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(Sensome_expressed), name = "Sensome_score", search = T)

Cluster_markers <- read_excel("Microglia cluster markers.xlsx")
Cluster_markers  <- lapply(Cluster_markers, function(x)x[!is.na(x)])

CL0 = Cluster_markers$CL0
CL1 = Cluster_markers$CL1
CL2 = Cluster_markers$CL2
CL3 = Cluster_markers$CL3
CL4 = Cluster_markers$CL4
CL5 = Cluster_markers$CL5
CL6 = Cluster_markers$CL6
CL7 = Cluster_markers$CL7
CL8 = Cluster_markers$CL8

# to add mouse cluster signatures first convert to human genes for each of the clusters


library(dplyr)
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

gene_list <- CL8
convert_mouse_to_human <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
}

CL8_human <- convert_mouse_to_human(gene_list)


human_myeloid  <- AddModuleScore(human_myeloid , list(CL0_human), name = "CL0_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(CL1_human), name = "CL1_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(CL2_human), name = "CL2_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(CL3_human), name = "CL3_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(CL4_human), name = "CL4_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(CL5_human), name = "CL5_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(CL6_human), name = "CL6_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(CL7_human), name = "CL7_score", search = T)
human_myeloid  <- AddModuleScore(human_myeloid , list(CL8_human), name = "CL8_score", search = T)



FeaturePlot(human_myeloid, features = c("CL0_score","CL1_score", "CL2_score", "CL3_score", "CL4_score", 
                                        "CL5_score","CL6_score", "CL7_score", "CL8_score","migration_score",
                                        "Phagocytosis_score", "APCs_score"))



#Fig. S8E
library(Seurat)
library(Azimuth)
library(patchwork)
library(tidyverse)
library(SPATA2)

source('modified_azimuth.R')

ref <- readRDS('~/azimuth_core_GBmap.rds')

#Check integration

ref@meta.data$annotation_level_4 <- ref@meta.data$annotation_level_4 %>% str_replace_all(., " ", "_") %>% str_replace_all(., "-", "_") %>% str_replace_all(., "/", "_")
ref@meta.data$annotation_level_4 %>% table() %>% as.data.frame() %>% arrange(Freq)
ref.new <- subset(ref, cells=rownames(ref@meta.data)[ref@meta.data$annotation_level_4!="Neuron"])
rm(ref)




### Load LUX data 
setwd("~/Desktop/Single_cell/Luxemburg")
myeloid_lux <- readRDS("myeloid_object.rds")
DimPlot(myeloid_lux, group.by = "New_names")

FeaturePlot(myeloid_lux, feature="B2m")

library(dplyr)
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

human.only <- 
  mouse_human_genes %>% 
  dplyr::select(DB.Class.Key,Common.Organism.Name, Symbol) %>% 
  filter(Common.Organism.Name=="human") %>% 
  rename("human":=Symbol)

mice.only <- 
  mouse_human_genes %>% 
  dplyr::select(DB.Class.Key,Common.Organism.Name, Symbol) %>% 
  filter(Common.Organism.Name=='mouse, laboratory') %>% 
  rename("mice":=Symbol)

genes_mice <-data.frame(mice=myeloid_lux@assays$RNA@scale.data %>% rownames())
genes_mice <- genes_mice %>% left_join(., mice.only, by="mice") %>% dplyr::select(mice, DB.Class.Key)

#add human
genes_mice <- genes_mice %>% left_join(.,human.only, by="DB.Class.Key" ) %>% dplyr::select(mice, human)

#Remove all genes that are not match
genes_mice <- genes_mice %>% filter(!is.na(human))

mat <- myeloid_lux@assays$RNA@counts %>% as.matrix()
mat <- mat[genes_mice$mice,]
rownames(mat) <- genes_mice$human
RNA <- Seurat::CreateAssayObject(mat)

new <- CreateSeuratObject(RNA)
new@meta.data <- myeloid_lux@meta.data

new[["percent.mt"]] <- PercentageFeatureSet(new, pattern = "^MT-")
new[["percent.RB"]] <- PercentageFeatureSet(new, pattern = "^RPS")
new <- Seurat::SCTransform(new,vars.to.regress = c("percent.mt"), return.only.var.genes = F)

new@reductions <- myeloid_lux@reductions
DimPlot(new, group.by = "New_names")


ds <- new
DimPlot(ds, size=0.01, group.by = "New_names")+scale_color_manual(values=colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(20))


query <- Seurat::DietSeurat(ds)
DefaultAssay(query) <- "RNA"
query@assays$SCT <- NULL
ds.ref <- 
  RunAzimuth(query,
             reference = ref.new, 
             annotation.levels = c('annotation_level_3', 'annotation_level_4'))


## Visualization
p1 <- DimPlot(ds.ref, group.by = 'predicted.annotation_level_3', reduction = 'ref.umap',
              label = TRUE, label.size = 3) + NoLegend()
p1+scale_color_manual(values=colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(17))

p2 <- DimPlot(ds.ref, group.by = 'predicted.annotation_level_4', reduction = 'ref.umap',
              label = TRUE, label.size = 3) + NoLegend()
p2+scale_color_manual(values=sub_color_space$colors)


sub_color_space <- ref.new@meta.data %>% as.data.frame() %>% SPATAwrappers::getSubColors(., pal="Set3",add_n=2, class = 'annotation_level_3',
                                                                                         subclass = "annotation_level_4")

names(sub_color_space)[2] <- "predicted.annotation_level_4"
saveRDS(sub_color_space, "~/Desktop/ImmunoSpatial/Paper/colors_cell_deconv2.RDS")

# Remove all colors not in annotations
sub_color_space <- sub_color_space %>% filter(predicted.annotation_level_4 %in% unique(ds.ref@meta.data$predicted.annotation_level_4))

sub_color_space <- SPATAwrappers::getSubColors(sub_color_space, class = 'annotation_level_3',
                                               subclass = "predicted.annotation_level_4")



plot.df <- 
  ds.ref@meta.data %>% 
  mutate(UMAP1 =ds.ref@reductions$ref.umap@cell.embeddings[,1],
         UMAP2 =ds.ref@reductions$ref.umap@cell.embeddings[,2]) %>% 
  left_join(., sub_color_space, by="predicted.annotation_level_4")

ggplot()+
  geom_point(mapping = aes(x=ref.new@reductions$umap@cell.embeddings[,1], 
                           y=ref.new@reductions$umap@cell.embeddings[,2]),
             size=0.1, color="#EDEDED", alpha=0.8)+
  geom_point(plot.df, mapping=aes(x=UMAP1, y=UMAP2), color=plot.df$colors, size=0.8)+
  theme_classic()+
  coord_fixed()

ggplot()+
  geom_point(mapping=aes(x=1:nrow(sub_color_space), y=1), color=sub_color_space$colors, size=4)+
  geom_text(mapping=aes(x=1:nrow(sub_color_space), y=1.2, label=sub_color_space$predicted.annotation_level_4),
            size=3, angle=90, nudge_y=0.1)+
  theme_void()


plot.df$predicted.annotation_level_4 <- factor(plot.df$predicted.annotation_level_4, 
                                               levels = sub_color_space$predicted.annotation_level_4)



SPATAImmune::plotBarCompare(plot.df, "predicted.annotation_level_4", "New_names")+scale_fill_manual(values=sub_color_space$colors)





setwd("~/Desktop/SpatialTranscriptomics/Visium/Visium/CancerCell_Revision/Deconvolution")
scDF <- readRDS("275_scDF.RDS")

scDF <- scDF %>% filter(celltypes %in% sub_color_space$predicted.annotation_level_4)

scDF$celltypes <- factor(scDF$celltypes , levels= sub_color_space$predicted.annotation_level_4)

myeloid_tye <- sub_color_space %>% filter(annotation_level_3 %in% c("TAM-BDM", "TAM-MG","Mono" ))


myeloid_tye <- sub_color_space %>% filter(annotation_level_3 %in% c("TAM-MG"))

ggplot()+
  geom_point(data=scDF %>% filter(celltypes %in% myeloid_tye$predicted.annotation_level_4), mapping=aes(x=x, y=y, color=celltypes), size=1)+
  theme_classic()+
  coord_fixed()+
  scale_color_manual(values=myeloid_tye$colors)


myeloid_tye <- sub_color_space %>% filter(annotation_level_3 %in% c("TAM-MG"))
spata.obj %>% SPATAwrappers::plotSurfaceMixed(., 
                                              normalize=F,
                                              Mixed = myeloid_tye$predicted.annotation_level_4,
                                              mixed.colors = c("Reds", "Blues", "Greens", "Greys"))

myeloid_tye <- sub_color_space %>% filter(annotation_level_3 %in% c("TAM-BDM"))
spata.obj %>% SPATAwrappers::plotSurfaceMixed(., 
                                              normalize=F,
                                              Mixed = myeloid_tye$predicted.annotation_level_4,
                                              mixed.colors = c("Reds", "Blues", "Greens", "Greys"))

myeloid_tye <- sub_color_space %>% filter(annotation_level_3 %in% c("Mono"))
spata.obj %>% SPATAwrappers::plotSurfaceMixed(., 
                                              normalize=F,
                                              Mixed = myeloid_tye$predicted.annotation_level_4,
                                              mixed.colors = c("Reds", "Blues", "Greens", "Greys"))



spatialSub <- c("Module_consensus_1",
                "Module_consensus_2",
                "Module_consensus_4",
                "Module_consensus_3",
                "Module_consensus_5")

spata.obj %>% addFeatures(., joinWith(., gene_sets = spatialSub, normalize = F, smooth = T) %>% 
                            dplyr::select(-sample, -x,-y,-row,-col)) %>% 
  SPATAwrappers::plotSurfaceMixed(., 
                                  smooth = T,
                                  display_image = T,
                                  pt_size = 2.4,
                                  Mixed = spatialSub,
                                  mixed.colors = c("Reds", "Blues", "Greens", "Greys", "Purples"))




### Supplementary figure 9

#Fig. S9A

rm(list=ls())
graphics.off()
# library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("~/projects/Yahaya/Figures-2022.09.29/")

#######@
# https://glioblastoma.alleninstitute.org/api/v2/well_known_file_download/305873915
path2ivya_exprs = "ivygap_all-2022.09.19/gene_expression_matrix_2014-11-25.zip"

# https://glioblastoma.alleninstitute.org/api/v2/gbm/rna_seq_samples_details.csv
path2ivya_rna_seq_samples_details = "ivygap_all-2022.09.19/rna_seq_samples_details.csv"

#
path2Supp_table = "yahaya/Yabo et al, Supplementary tables.xlsx"

####### This code needs code of https://github.com/livnatje/ImmuneResistance
# to calculate OE score
# Download the code and change the path below
path2immune_res_code = '~/pipeline/scRNASeq/ImmuneResistance/Code/'

#######
tmp = unzip(path2ivya_exprs)

stopifnot(file.exists("fpkm_table.csv"))
stopifnot(file.exists("columns-samples.csv"))
stopifnot(file.exists("rows-genes.csv"))

## read fpkm
ivygap_x = read.csv(
  "fpkm_table.csv",header = T,row.names = 1,stringsAsFactors = F,check.names = F
)
ivygap_x = as.matrix(ivygap_x)
dim(ivygap_x)
range(ivygap_x)
hist(rowMeans(log2(ivygap_x+1)))

## read first meta
ivygap_meta = read.csv(
  "columns-samples.csv",header = T,stringsAsFactors = F,check.names = F
)
all(colnames(ivygap_x) == ivygap_meta$rna_well_id)
rownames(ivygap_meta) = ivygap_meta$rna_well_id

## read meta
ivygap_meta2 = read.csv(
  path2ivya_rna_seq_samples_details,header = T,stringsAsFactors = F,check.names = F
)

## merge 2 sample info
dim(ivygap_meta2) #279  22
dim(ivygap_meta)  #270  12

table(ivygap_meta$rna_well_id %in% ivygap_meta2$sample_id)
table(colnames(ivygap_x) %in% ivygap_meta2$sample_id)
ivygap_meta2 = ivygap_meta2[ivygap_meta2$sample_id %in% ivygap_meta$rna_well_id,]
ivygap_meta2 = ivygap_meta2[match(ivygap_meta$rna_well_id,ivygap_meta2$sample_id),]
stopifnot(all(ivygap_meta$rna_well_id == ivygap_meta2$sample_id))

common.cols = intersect(colnames(ivygap_meta),colnames(ivygap_meta2))
names(common.cols) = common.cols
stopifnot(all(sapply(common.cols, function(x) all(ivygap_meta[,x] == ivygap_meta2[,x]))))

ivygap_meta = merge(
  ivygap_meta,
  ivygap_meta2[,setdiff(colnames(ivygap_meta2),common.cols)],
  by.x = "rna_well_id", by.y = "sample_id"
)
rm(ivygap_meta2)
rownames(ivygap_meta) = ivygap_meta$rna_well_id

## aggragate region and simplified name
ivygap_meta$structure_name_abbr = ivygap_meta$structure_name
ivygap_meta$structure_name_abbr[grepl("^Cellular Tumor",ivygap_meta$structure_name_abbr)] = "Cellular Tumor"
ivygap_meta$structure_name_abbr[grepl("^Hyperplastic blood vessels",ivygap_meta$structure_name_abbr)] = "HyBV"
ivygap_meta$structure_name_abbr[grepl("^Infiltrating Tumor",ivygap_meta$structure_name_abbr)] = "Inf. Tumor"
ivygap_meta$structure_name_abbr[grepl("^Leading Edge",ivygap_meta$structure_name_abbr)] = "LD Edge"
ivygap_meta$structure_name_abbr[grepl("^Microvascular proliferation",ivygap_meta$structure_name_abbr)] = "MvP"
ivygap_meta$structure_name_abbr[grepl("^Perinecrotic zone sampled",ivygap_meta$structure_name_abbr)] = "Perinecrotic"
ivygap_meta$structure_name_abbr[grepl("^Pseudopalisading",ivygap_meta$structure_name_abbr)] = "Pseudopalisading"
table(ivygap_meta$structure_name_abbr)

# color for each region
structure_name_abbr.cols = c(
  "LD Edge" = "blue",
  'Inf. Tumor'="purple",
  "Cellular Tumor"="darkgreen",
  Perinecrotic="lightblue",
  Pseudopalisading="lightgreen",
  HyBV ="orange",
  MvP="red"
)

# order sample according to region, tumor_id and tumor_name
ivygap_meta = ivygap_meta[
  order(
    match(ivygap_meta$structure_name_abbr,names(structure_name_abbr.cols)),
    ivygap_meta$tumor_id,
    ivygap_meta$tumor_name
  ),
]

# same sample and order as in metadata
ivygap_x = ivygap_x[,rownames(ivygap_meta)]
range(ivygap_x)

# simplify name of IVyGap study and set color for each
ivygap_meta$study = ifelse(ivygap_meta$study_name == "Anatomic Structures RNA Seq","Anat. Str.",'CSC')
study.cols = c("Anat. Str."="grey85",'CSC'="grey50")

# gene annotation
gns = read.csv(
  "rows-genes.csv",header = T,stringsAsFactors = F,check.names = F
)
stopifnot(all(gns$gene_id == rownames(ivygap_x)))
rownames(ivygap_x) = gns$gene_symbol

file.remove(c("rows-genes.csv","columns-samples.csv","fpkm_table.csv","README.txt"))
rm(gns)

#########
load("yahaya_gnsig.RData") # "unique_gene", "yahaya_gnsig_", "yahaya_gnsig" 

cluster_signature_from_mo = readxl::read_xlsx(
  path2Supp_table,
  sheet = "Table S4",
  skip  = 3,
  col_names = T
)
cluster_signature_from_mo = cluster_signature_from_mo[,-1]
cluster_signature_from_mo = cluster_signature_from_mo[-(1:5),]
cluster_signature_from_mo = lapply(cluster_signature_from_mo, function(x) x[!is.na(x)])

#### convert mouse to human gene names,
## use file instead of biomart to avoid
## 'unexpected server error'.
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
mm_hs = split(mouse_human_genes,mouse_human_genes$DB.Class.Key)
mm_hs = mm_hs[sapply(mm_hs,nrow)>1]
table(sapply(mm_hs,function(x) sum(x[,"Common.Organism.Name"] == "mouse, laboratory")))
mm_hs = mm_hs[sapply(mm_hs,function(x) sum(x[,"Common.Organism.Name"] == "mouse, laboratory"))==1]
length(mm_hs) # 20124
mm = sapply(mm_hs,function(x) x[x$Common.Organism.Name == "mouse, laboratory","Symbol"])
hs = lapply(mm_hs,function(x) x[x$Common.Organism.Name == "human","Symbol"])
mm_hs = data.frame(
  mm = rep(mm,sapply(hs,length)),
  hs = unlist(hs),
  stringsAsFactors = F
)
rm(hs,mm,mouse_human_genes)
# clean
mm_hs = mm_hs[
  mm_hs$hs %in% rownames(ivygap_x) & 
    mm_hs$mm %in% unlist(cluster_signature_from_mo),
]

# manage duplicate, privilege to same symbol while different 'cases'
ident = sapply(split(toupper(mm_hs$mm) == toupper(mm_hs$hs) ,mm_hs$mm),any)
mm_hs$symbol = ifelse(
  table(mm_hs$mm)[mm_hs$mm] > 1 & ident[mm_hs$mm],
  toupper(mm_hs$mm),
  mm_hs$hs
)

cluster_signature_from_mo_ = lapply(
  cluster_signature_from_mo,
  function(x){
    unique(mm_hs$hs[mm_hs$mm %in% x])
  }
)
cbind(
  Mm = sapply(cluster_signature_from_mo,length),
  Hs = sapply(cluster_signature_from_mo_,length)
)
# Mm  Hs
# CL0  59  59 <<<
# CL1  43  43 <<<
# CL2  24  21 <<<
# CL3 104  64 <<<
# CL4 225 215 <<<
# CL5 235 156 <<<
# CL6 254 228 <<<
# CL7 433 348 <<<
# CL8  65  63 <<<
cluster_signature_from_mo = cluster_signature_from_mo_
rm(cluster_signature_from_mo_,mm_hs,ident)

####### Calculate OE score
wd = getwd()
setwd(path2immune_res_code)
source('ImmRes_source.R')
setwd(wd)

all.x = log2(ivygap_x+1)
sc = list()
sc$tpm = all.x[rowSums(all.x)  > 0,]
sc$genes = rownames(sc$tpm)
oe_score = get.OE.sc(sc,gene.sign = c(yahaya_gnsig_,cluster_signature_from_mo),num.rounds = 500)
rownames(oe_score) = colnames(sc$tpm)
oe_score = data.frame(oe_score)
rm(sc)

######## set colors
molecular_subtype.cols = c(
  "Classical" = "red",
  "Classical, Mesenchymal"  = "blue",
  "Classical, Neural" = "black",
  "Proneural" = "purple",
  "Neural" = "darkgreen",
  "Neural, Proneural" = "grey",
  "Mesenchymal" = "yellow",
  "Mesenchymal, Neural" = "lightblue"
)

tumor.cols = rainbow(length(unique(ivygap_meta$tumor_name)))
names(tumor.cols) = unique(ivygap_meta$tumor_name)

########### Figure 5C
col_annot = data.frame(
  row.names = rownames(ivygap_meta),
  Features = ivygap_meta$structure_name_abbr,
  Study = ivygap_meta$study,
  MolSubT = ivygap_meta$molecular_subtype,
  stringsAsFactors = F
)

col_annot_col = list(
  Features = structure_name_abbr.cols,
  MolSubT  = molecular_subtype.cols,
  Study = study.cols
)

colAnno = HeatmapAnnotation(
  df = col_annot,
  which="col",
  col=col_annot_col,
  show_legend = F
)

ComplexHeatmap::Heatmap(
  t(
    oe_score[,
             c(
               "Human.Mg","Human.Mo","CL0","CL1","CL2","CL3",
               "CL4","CL5","CL6","CL7","CL8","Migration",
               "Sensome","Phagocytosis","APC"
             )
    ]
  ),
  show_row_names  = T,
  cluster_rows    = F,
  # row_split = factor(ifelse(
  #     rownames(all.x) %in% yahaya_gnsig$`Human Mg`,
  #     "Human Mg","Human Mo"
  # ),levels = c("Human Mg","Human Mo")),
  # row_gap = unit(2,"mm"),
  
  cluster_columns = F,
  show_column_names = F,
  column_split = factor(
    ivygap_meta$structure_name_abbr,
    levels = names(structure_name_abbr.cols)
  ),
  column_gap = unit(2,"mm"),
  
  # col = c(
  #     rev(brewer.pal(9,"Blues")[c(1:9)]),
  #     "white",
  #     brewer.pal(9,"Reds")[c(3:9)]
  # ),
  col = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#3182bd", "white", "darkred")),
  top_annotation = colAnno,
  heatmap_legend_param = list(title="Z-score")
)

########### figure S8A
par(mfrow=c(3,5),mar=c(1,5,3,1))
for(i in c(
  "Human.Mg","Human.Mo","CL0","CL1","CL2","CL3",
  "CL4","CL5","CL6","CL7","CL8","Migration",
  "Sensome","Phagocytosis","APC"
)){
  boxplot(
    oe_score[,i] ~ factor(
      ivygap_meta$structure_name_abbr,
      levels = names(structure_name_abbr.cols)
    ),
    xlab="",ylab="Over-expression Score",
    main=sub("[.]"," ",i),
    las=2,
    xaxt="n",
    col = structure_name_abbr.cols
  )
  # text(1:7,par("usr")[3]-.07,names(structure_name_abbr.cols),adj = 1,srt=45,xpd="n",cex=.8)
}




#Fig S9B


# Import functional signatures
func_sig <- read.csv("functional signatures.csv", sep=";")
gs_funk_sig <- map_dfr(.x=1:ncol(func_sig), .f=function(i){
  data.frame(ont=names(func_sig)[i], gene= func_sig %>% dplyr::select(!!sym(names(func_sig)[i])) %>% filter(!!sym(names(func_sig)[i])!="") %>% pull(!!sym(names(func_sig)[i])) )
})

path_visium=("~/Desktop/SpatialTranscriptomics/Visium/Visium")
setwd(path_visium)
meta.st <- read.csv( "feat.csv", sep=";")

setwd("~/Desktop/SpatialTranscriptomics/Visium/Visium/CancerCell_Revision/Modules_WT")
obj.wt <- meta.st %>% filter(Region=="T") %>% filter(Tumor=="IDH-WT") %>% pull(files) %>% unique()
obj.wt <-data.frame(files=obj.wt, 
                    path=paste0("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions/", "Revision_",obj.wt,"_SPATA_CNV_Pred.RDS"),
                    imh.hrs=paste0("~/Desktop/SpatialTranscriptomics/Visium/Visium/",obj.wt,"/outs/spatial/tissue_hires_image.png"))
obj.wt$Seuratstep1=NA

spata.obj <- readRDS(obj.wt$path[15])
spata.obj@used_genesets <- rbind(spata.obj@used_genesets, gs_funk_sig)

unique(gs_funk_sig$ont)

plotSurface(spata.obj, pt_alpha = 0)

spata.obj <- spata.obj %>% addFeatures(., joinWith(., gene_sets = c("Microglia", "Macrophages","Sensome", "Migration", "APC", "Phagocytosis" ), normalize = F, smooth = T) %>% 
                                         dplyr::select(-sample, -x,-y,-row,-col))

SPATAwrappers::plotSurfaceMixed(spata.obj, 
                                smooth = F,
                                display_image = F,
                                pt_size = 2.4,
                                Mixed = c("Microglia", "Macrophages" ),
                                mixed.colors = c("Reds", "Greens"))


#Fig S9C

SPATAwrappers::plotSurfaceMixed(spata.obj, 
                                smooth = F,
                                display_image = F,
                                pt_size = 2.4,
                                Mixed = c("Sensome", "Migration", "APC", "Phagocytosis" ),
                                mixed.colors = c("Reds", "Greens", "Blues", "Oranges"))







### Supplementary figure 10



