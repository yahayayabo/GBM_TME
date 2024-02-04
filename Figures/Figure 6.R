
####  Figure 6  #####

library(Seurat)
library(Azimuth)
library(patchwork)
library(dittoSeq)
library(tidyverse)
library(ggpubr)
library(Matrix)
library(readxl)
library(reshape2)
library(scrabble)


#load data
sobj <- readRDS("azimuth_core_GBmap.rds")

#Fig 6A

#subset only TAM-BDM, TAM-MG and Mono

sobj_myeloid <- subset(sobj, subset=annotation_level_3%in%c("Mono","TAM-BDM","TAM-MG"))
dim(sobj_myeloid)


DimPlot(sobj_myeloid, reduction = "umap")

#score cell base on our signatures

#score cells using scrabble
Human_GEM <- as.matrix(sobj@assays$RNA@counts)


#load Human Mg and Human Mo gene lists from Table S6
Mg_Mo_TAM_sig <- read.table("Human_Mg_Human_Mo_gene_lists.csv", header=TRUE,sep=";")


#Set column headers as id variables
Mg_Mo_TAM_sig <- melt(Mg_Mo_TAM_sig)

#create separate vectors containing the Suva signature genes
#Mg
MG_genes <- as.vector(Mg_Mo_TAM_sig$Mg.TAM[Mg_Mo_TAM_sig$Mg.TAM!=""])

#Mo
Mo_genes <- as.vector(Mg_Mo_TAM_sig$Mo.TAM[Mg_Mo_TAM_sig$Mo.TAM!=""])



#get only the unique signature genes that are present in the normalized and filtered dataset
Mg.related.genes<-unique(MG_genes[MG_genes%in%rownames(Human_GEM)])

Mo.related.genes<-unique(Mo_genes[Mo_genes%in%rownames(Human_GEM)])


#create a list containing all signature related genes
vector <- list(Mg.related.genes, Mo.related.genes)

#check class (supposed to be a list; package requirement)
class(vector)


TAM_Scores <- scrabble::score(Human_GEM, vector, center = F, nbin = 30, n = 100, replace = F)


#check class (it should be a matrix)
class(TAM_Scores)

colnames(TAM_Scores) <- c("Mg.TAM", "Mo.TAM")

#write output
write.table(TAM_Scores, file = "TAM_Scores_GBmap.txt", sep="\t")

rows <- apply(TAM_Scores, 1 , function(x) any( x > 1 ))

criteria_one <- apply(TAM_Scores,1,function(x) colnames(TAM_Scores)[which.max(x)])


metadata <- as.data.frame(criteria_one)

sobj_myeloid <- AddMetaData(sobj_myeloid, metadata, col.name = "TAMScores")


DimPlot(sobj_myeloid, reduction = "umap", group.by = "TAMScores", cols = c('Mg.TAM' = '#4daf4a', 'Mo.TAM' = 'blue'))

Idents(sobj_myeloid) <- "TAMScores"

dittoBarPlot(sobj_myeloid, "ident", group.by = "patient",
             color.panel = c('Mg.TAM' = '#4daf4a', 'Mo.TAM' = 'blue'))


dittoBarPlot(sobj_myeloid, "ident", group.by = "patient",
             var.labels.reorder = c(2, 1),
             x.reorder = c(23,65,91,8,66,64,47,31,33,16,
                           95,63,29,34,94,60,32,80,87,50,
                           96,2,71,103,69,30,42,19,61,81,
                           98,45,56,55,39,85,62,13,18,84,
                           101, 46, 5,25,1,40,97,58,4,89,
                           17, 99,83,38,11, 35,7,70,28,24,
                           15,78,88,27,59,14,57,10,
                           41,48,77,6,37,26,82,68,72,54,100,12,44,49,76,
                           36,90,74,93,86,51,102,9,67,79,
                           73,3,92,20,75,22,52,43,53,21),
             color.panel = c('Mo.TAM' = 'blue', 'Mg.TAM' = '#4daf4a'))





#Fig 6C
rm(list=ls())
graphics.off()
# library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)


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
# zip contains : 
#   fpkm_table.csv 
#   columns-samples.csv
#   rows-genes.csv
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

#Fig 6D and E


library(SPATA2)
library(tidyverse)
library(Seurat)
library(tidyverse)

# Import the myeloid signatures and the reference dataset from verhaak
pan_glioma <- readRDS("~/Desktop/ImmunoSpatial/PanGlioma/Seurat_panglioma.RDS")

# Isolate Myeloid cells and subcluster
pan_glioma@meta.data %>% as.data.frame() %>% pull(cell_state) %>% unique()
myeloid <- pan_glioma@meta.data %>% as.data.frame() %>% filter(cell_state=="Myeloid") %>% rownames()

#Isolate object
myeloid <- subset(pan_glioma, cells=myeloid)
myeloid <- myeloid %>% Seurat::DietSeurat()

# Horizontal integration by MNN
myeloid <- SeuratWrappers::RunFastMNN(Seurat::SplitObject(myeloid, split.by = "sampleid"))

#QC
myeloid[["percent.mt"]] <- PercentageFeatureSet(myeloid, pattern = "^MT-")
myeloid[["percent.RB"]] <- PercentageFeatureSet(myeloid, pattern = "^RPS")

VlnPlot(myeloid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Remove Mitochondrial and Stress genes
exclude=c(rownames(myeloid@assays$RNA)[grepl("^RPL", rownames(myeloid@assays$RNA))],
          rownames(myeloid@assays$RNA)[grepl("^RPS", rownames(myeloid@assays$RNA))],
          rownames(myeloid@assays$RNA)[grepl("^MT-", rownames(myeloid@assays$RNA))],
          c('JUN','FOS','ZFP36','ATF3','HSPA1A","HSPA1B','DUSP1','EGR1','MALAT1'))

# @Method and genes from Sankowski et al., Nat Neuroscience

feat_keep=rownames(myeloid@assays$RNA[!(rownames(myeloid@assays$RNA) %in% exclude), ])
myeloid=subset(myeloid, features=feat_keep)

# PostProcessing
myeloid <- 
  myeloid %>% 
  Seurat::SCTransform(vars.to.regress = c("percent.mt"), return.only.var.genes = F)

myeloid <- Seurat::FindNeighbors(myeloid, reduction="mnn")

# Find optimal clusters
myeloid <- SPATAwrappers::run.SNN.stability(myeloid,reduction = "mnn", assay="RNA")
myeloid <- Seurat::RunUMAP(myeloid, dims=1:20, reduction="mnn")
DimPlot(myeloid)


# Read in LUX myeliod file:
setwd("~/Desktop/Single_cell/Luxenburg")
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


#Integrate to human reference

myeloid@active.assay <- "SCT"
myeloid <- FindNeighbors(myeloid, reduction = "mnn")
myeloid <- myeloid %>% RunPCA()


anchors <- FindTransferAnchors(
  reference = myeloid,
  query = new,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30
)


myeloid <- myeloid %>% RunUMAP(reduction="mnn", dims=1:30,return.model=T)
DimPlot(myeloid, reduction = "umap")

map <- MapQuery(
  anchorset = anchors,
  query = new,
  reference = myeloid,
  refdata = list(
    new = "new"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(map, group.by = "New_names", reduction = "ref.umap")


# From Markers

marker_verh <- Seurat::FindAllMarkers(myeloid, 
                                      logfc.threshold = 0.5, 
                                      min.pct = 0.4, 
                                      min.cells.group = 100, 
                                      densify=T)
marker_verh_order <- marker_verh %>% group_by(cluster) %>% top_n(200, wt="avg_log2FC")

# Compare the Luxenburg cluster with the myeloid subclustering
setwd("~/Desktop/Single_cell/Luxenburg")
marker_lux <- read.csv("Myeloid_cluster markers.csv", sep=";")
marker_lux <- marker_lux %>% filter(human.genes!="")


# 1. Compare by jaccard-Index

LX <- paste0("LX_", unique(marker_lux$cluster))
VER <- paste0("Verhaak_", unique(marker_verh$cluster))
Cor_mat <- matrix(0, nrow=length(LX), ncol=length(VER))
rownames(Cor_mat) <- VER; colnames(Cor_mat) <- LX

for(i in 1:length(VER)){
  for(j in 1:length(LX)){
    
    Cor_mat[VER[i], LX[j]] <- 
      SPATAwrappers::jaccard(
        marker_lux %>% filter(cluster==unique(marker_lux$cluster)[j]) %>% pull(human.genes),
        marker_verh_order %>% filter(cluster==unique(marker_verh_order$cluster)[i]) %>% pull(gene)
      )
    
    
  }
}
corrplot::corrplot(Cor_mat, is.corr = F, order = "hclust")


# 2. Compare by Gene expression

LX <- paste0("LX_", unique(marker_lux$cluster))
VER <- paste0("Verhaak_", unique(marker_verh$cluster))
Cor_mat <- matrix(0, nrow=length(LX), ncol=length(VER))
rownames(Cor_mat) <- VER; colnames(Cor_mat) <- LX

for(i in 1:length(VER)){
  for(j in 1:length(LX)){
    
    genes <- marker_lux %>% filter(cluster==unique(marker_lux$cluster)[j]) %>% pull(human.genes)
    mat.exp <- Seurat::GetAssayData(myeloid)
    sub <- mat.exp[rownames(mat.exp) %in% genes, ] %>% colMeans() %>% as.data.frame()
    cells <- myeloid@meta.data %>% as.data.frame() %>% filter(new==unique(marker_verh$cluster)[i]) %>% rownames()
    Cor_mat[VER[i], LX[j]] <- sub[cells, ] %>% mean()
    
    
    
  }
}
corrplot::corrplot(Cor_mat, is.corr = F, order = "hclust")



# Add myeloid subgroups to pan glioma
names <- marker_verh %>% group_by(cluster) %>% top_n(20, wt=avg_log2FC)
cluster_names <- data.frame(new=unique(names$cluster),
                            names=c("CX3CR1_pos","Inflammatory_Microglia","Homeostasis", "Monocytes", "Stress_induced","M1","NOS_1", "Hypoxia_induced", "M1" ))


myeloid_def <- myeloid@meta.data %>% as.data.frame() %>% rownames_to_column("barcodes") %>% left_join(., cluster_names, "new") %>% dplyr::select(barcodes,new, names)
myeloid@meta.data$cell_type <- myeloid_def$names

DimPlot(myeloid, group.by = "cell_type", label = T, raster=T, reduction="umap")+scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(8))


pan_glioma@meta.data[myeloid_def$barcodes, ]$cell_state <- myeloid_def$names
DimPlot(pan_glioma, group.by = "cell_state", label = T, raster=T)+scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(20))



# Import the SPATA datasets

path_visium=("~/Desktop/SpatialTranscriptomics/Visium/Visium")
setwd(path_visium)
meta.st <- read.csv( "feat.csv", sep=";")

setwd("~/Desktop/SpatialTranscriptomics/Visium/Visium/CancerCell_Revision/Modules_WT")
obj.wt <- meta.st %>% filter(Region=="T") %>% filter(Tumor=="IDH-WT") %>% pull(files) %>% unique()
obj.wt <-data.frame(files=obj.wt, 
                    path=paste0("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions/", "Revision_",obj.wt,"_SPATA_CNV_Pred.RDS"),
                    imh.hrs=paste0("~/Desktop/SpatialTranscriptomics/Visium/Visium/",obj.wt,"/outs/spatial/tissue_hires_image.png"))
obj.wt$Seuratstep1=NA

setwd("~/Desktop/Single_cell/Luxenburg")

#Create GeneSet from marker genes (top 50)
marker_gs <- 
  marker_lux %>% 
  group_by(cluster) %>% 
  top_n(n=100, wt=avg_log2FC) %>% 
  ungroup() %>% 
  dplyr::select(cluster, human.genes)

names(marker_gs) <- c("ont", "gene")

# Filter to genes that are in SPATA
marker_gs %>% group_by(ont) %>% summarise(sum=length(ont))

# Import functional signatures
func_sig <- read.csv("functional signatures.csv", sep=";")
gs_funk_sig <- map_dfr(.x=1:ncol(func_sig), .f=function(i){
  data.frame(ont=names(func_sig)[i], gene= func_sig %>% dplyr::select(!!sym(names(func_sig)[i])) %>% filter(!!sym(names(func_sig)[i])!="") %>% pull(!!sym(names(func_sig)[i])) )
})



CCA_matrix <- 
  map(.x=1:5, .f=function(i){
    
    spata <- readRDS(obj.wt$path[i])
    samples=getSampleNames(spata)
    names(spata@fdata[[samples]])
    
    
    #spata@used_genesets <- rbind(gs.all , marker_gs, gs_funk_sig)
    
    
    
    #plot some clusters
    #col <- colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Spectral")))
    
    #color.pram <- unique(marker_gs$ont)[1]
    #plotSurface(spata, color_by = color.pram, alpha_by = color.pram, display_image = T, smooth = T) +
    #  scale_colour_gradientn(colours = col(50), oob = scales::squish, limit=c(0.6, 1))
    
    #plotSurfaceComparison(spata, color_by = unique(marker_gs$ont), display_image = F, smooth = T) +
    #  scale_colour_gradientn(colours = col(50), oob = scales::squish, limit=c(0.4, 1))
    
    
    
    
    
    #Run myeloid subtype analysis
    
    
    # Get CCA correlation of all clusters 
    gs <- c("MILO_A1","MILO_A2" ,"MILO_fetal", "MILO_adult",
            "Verhaak_Classical","Verhaak_Mesenchymal","Verhaak_Proneural","Verhaak_Neural",
            "Neftel_MESlike2","Neftel_MESlike1","Neftel_AClike","Neftel_OPClike","Neftel_NPClike1","Neftel_NPClike2",
            "Module_consensus_1","Module_consensus_2","Module_consensus_4","Module_consensus_3","Module_consensus_5",
            "Developmental", "Injury_Response")
    
    
    subgroups <- factor(c("Stem-like","Diff.-like","Stress_induced","Prolif. stem-like","Oligodendrocyte","Monocytes","Granulocyte","Inflammatory_Microglia",
                          "T cell","Endothelial","M1","Dendritic cell","CX3CR1_pos","Fibroblast","Homeostasis","Pericyte","B cell","Hypoxia_induced", "NOS_1") %>% 
                          str_replace_all(., "-", "_") %>% 
                          str_replace_all(., " ", "_"), levels = c("Prolif. stem-like","Stem-like", "Diff.-like",
                                                                   "T cell","B cell","Dendritic cell",
                                                                   "Homeostasis","Inflammatory_Microglia","CX3CR1_pos","Hypoxia_induced","Stress_induced","M1","Monocytes","NOS_1",
                                                                   "Granulocyte","Endothelial","Fibroblast","Oligodendrocyte","Pericyte") %>% 
                          str_replace_all(., "-", "_") %>% 
                          str_replace_all(., " ", "_"))
    samples=getSampleNames(spata)
    message(paste0("Run CCA Analysis: ", samples))
    
    feat <- c(as.character(subgroups),
              gs,
              unique(marker_gs$ont),
              unique(gs_funk_sig$ont))
    
    mat.cor <- SPATAImmune::getSpatialRegression(spata, 
                                                 features = feat,
                                                 model = "CCA")
    
    return(mat.cor)
    
    
    
  })



cca.matrix <- Reduce(`+`, CCA_matrix) / length(CCA_matrix)


corrplot::corrplot(cca.matrix[c("Homeostasis","Inflammatory_Microglia","CX3CR1_pos","Hypoxia_induced","Stress_induced","M1","Monocytes","NOS_1"), 
                              unique(marker_gs$ont)], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Reds")))(50))


pram=c("Prolif._stem_like","Stem_like", "Diff._like")

corrplot::corrplot(cca.matrix[pram, 
                              unique(marker_gs$ont)], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Greens")))(50))

corrplot::corrplot(cca.matrix[pram, 
                              c("Homeostasis","Inflammatory_Microglia","CX3CR1_pos","Hypoxia_induced","Stress_induced","M1","Monocytes","NOS_1")], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Greens")))(50))

pram=c("T_cell","B_cell","Dendritic_cell",
       "Granulocyte","Endothelial","Fibroblast","Oligodendrocyte","Pericyte")

corrplot::corrplot(cca.matrix[pram, 
                              unique(marker_gs$ont)], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Oranges")))(50))

corrplot::corrplot(cca.matrix[pram, 
                              c("Homeostasis","Inflammatory_Microglia","CX3CR1_pos","Hypoxia_induced","Stress_induced","M1","Monocytes","NOS_1")], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Oranges")))(50))



pram=c("Module_consensus_1","Module_consensus_2","Module_consensus_4","Module_consensus_3","Module_consensus_5")

corrplot::corrplot(cca.matrix[pram, 
                              unique(marker_gs$ont)], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Blues")))(50))

corrplot::corrplot(cca.matrix[pram, 
                              c("Homeostasis","Inflammatory_Microglia","CX3CR1_pos","Hypoxia_induced","Stress_induced","M1","Monocytes","NOS_1")], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Blues")))(50))


pram=c("Neftel_MESlike2","Neftel_MESlike1","Neftel_AClike","Neftel_OPClike","Neftel_NPClike1","Neftel_NPClike2")

corrplot::corrplot(cca.matrix[pram, 
                              unique(marker_gs$ont)], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Purples")))(50))

corrplot::corrplot(cca.matrix[pram, 
                              c("Homeostasis","Inflammatory_Microglia","CX3CR1_pos","Hypoxia_induced","Stress_induced","M1","Monocytes","NOS_1")], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Purples")))(50))

pram=c("Microglia", "Macrophages", "Homeostatic.microglia",    "Sensome",  "Migration","APC", "Phagocytosis")

corrplot::corrplot(cca.matrix[pram, 
                              unique(marker_gs$ont)], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Greys")))(50))

corrplot::corrplot(cca.matrix[pram, 
                              c("Homeostasis","Inflammatory_Microglia","CX3CR1_pos","Hypoxia_induced","Stress_induced","M1","Monocytes","NOS_1")], is.corr = F,
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Greys")))(50))



spata <- readRDS(obj.wt$path[1])
plotSurfaceInteractive(spata)

pram=c("Microglia", "Macrophages", "Homeostatic.microglia",    "Sensome",  "Migration","APC", "Phagocytosis")

plotSurfaceComparison(spata, color_by = pram)







