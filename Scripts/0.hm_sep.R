# ---
# title: "GBM_TME"
# author: "Dimitrios Kyriakis"
# date: "26/6/2020"
# ---

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
library(useful)
# ----------------------------------------------------------------------------------


# ================================ Color Pelet =====================================
color_cond <- c( "magenta4", "#007A87",brewer.pal(6,"Dark2")[-1],"#FF5A5F","black",
    "#FFB400", "#007A87", "#8CE071", "#7B0051",
                 "#00D1C1")
color_clust <- c("#FF5A5F", "#FFB400", "#007A87", "#8CE071", "#7B0051",
                 "#00D1C1", "#FFAA91", "#B4A76C", "#9CA299", "#565A5C", "#00A04B", "#E54C20",brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust,Oligo_Pop=color_clust)
# ----------------------------------------------------------------------------------


# =============================== Load Functions ====================================
script_libs <- list.files("/home/users/dkyriakis/PhD/Projects/Yahaya/Script_Library",full.names = TRUE)
lapply(script_libs,source)
# ----------------------------------------------------------------------------------

# ============================= PATH OF FILES =======================================
# Produced by Kamil (GBMST)
path = "/work/projects/esprit/Dropseq_projects/GBMST/mouse/counts"
filenames_old <-list.files(path,full.names = TRUE)
filenames_DGE_old <- filenames_old[grep("Summary|readcounts",filenames_old,invert=TRUE)][c(1,3,5)]
cond_names_old <- c("P13S","T16S","T192S")

# Produced by Kamil New (LIH2)
path_new = "/work/projects/esprit/Dropseq_projects/LIH2/new/mouse/counts"
filenames_new <-list.files(path_new,full.names = TRUE)
filenames_DGE_new <- filenames_new[grep("Summary",filenames_new,invert=TRUE)][c(1,2)]
cond_names_new <- c("T470S","T101S")

filenames <- c(filenames_DGE_old,filenames_DGE_new)
condition_names <- c(cond_names_old,cond_names_new)


# New data 15/07/2020 from Kamil (LIHII)
path_new<-"/work/projects/esprit/Dropseq_projects/LIHII/mouse/counts"
filenames_new <-list.files(path_new,full.names = TRUE)
filenames_DGE_new <- filenames_new[grep("Summary|readcounts",filenames_new,invert=TRUE)]
cond_names_new <- c("T347S","T233S","P3TMZ","P3Con")

filenames <- c(filenames,filenames_DGE_new)
condition_names <- c(condition_names,cond_names_new)


old_suresh <- list.files("/home/users/dkyriakis/PhD/Projects/Yahaya/DATA_old",full.names=T)
filenames <- c(filenames,old_suresh)
condition_names <- c(condition_names,c("Control","p13_old","p3_old","p8"))
# -----------------------------------------------------------------------------------

setwd("/home/users/dkyriakis/PhD/Projects/Yahaya/")


test <- readRDS("/home/users/dkyriakis/PhD/Projects/Yahaya/4.Clustering/Liger_Merged_CL.rds")
backup_obj<-test



dir.create("5.Seperation")
setwd("5.Seperation")
WORKDIR<-"/home/users/dkyriakis/PhD/Projects/Yahaya/"

human_dir <- paste0(WORKDIR,"/human/")
mouse_dir <- paste0(WORKDIR,"/mouse/")

list.files(human_dir,full.names = TRUE)


order_files<-c("4","5","6","7","8","9","1","2","3")
iter=1
human_file <- NULL
for (file in list.files(human_dir,full.names = TRUE)){
  print(iter)
  print(file)
  temp_file1 <- read.table(file,header = TRUE)
  if(dim(temp_file1)[2]==4){
    temp_file1[2]<-NULL
  }
  temp_file1[,1] <- paste0(temp_file1[,1],"_",order_files[iter]) 
  print(dim(temp_file1))
  human_file <- rbind(human_file,temp_file1)
  iter=iter+1
}

iter=1
mouse_file <- NULL
for (file in list.files(mouse_dir,full.names = TRUE)){
  temp_file1 <- read.table(file,header = TRUE)
  if(dim(temp_file1)[2]==4){
    temp_file1[2]<-NULL
  }
  temp_file1[,1] <- paste0(temp_file1[,1],"_",order_files[iter]) 
  print(dim(temp_file1))
  mouse_file <- rbind(mouse_file,temp_file1)
  iter=iter+1
}


corner(human_file)
corner(mouse_file)
human_df <- human_file[human_file[,1] %in% colnames(backup_obj),]
mouse_df <- mouse_file[mouse_file[,1] %in% colnames(backup_obj),]
dim(mouse_df)
dim(human_df)
corner(human_df)

human_df[human_df$CELL_BARCODE=="AAAACGGTTGTT_1",]        
mouse_df[mouse_df$CELL_BARCODE=="AAAACGGTTGTT_1",]        
merged_df <- merge(mouse_df,human_df,by="CELL_BARCODE", all = TRUE)
merged_df[is.na(merged_df)] <- 0
merged_df_ord <- merged_df[match(colnames(backup_obj),merged_df$CELL_BARCODE),]
merged_df_ord$Cluster <- backup_obj$Cluster
colnames(merged_df_ord) <- c("CELL_BARCODE",
                             "Mouse_NUM_GENES",
                             "Mouse_NUM_TRANSCR",
                             "Human_NUM_GENES", 
                             "Human_NUM_TRANSCR", "Cluster")
merged_df_ord$Diff_tr <- merged_df_ord$Mouse_NUM_TRANSCR-merged_df_ord$Human_NUM_TRANSCR
merged_df_ord$Diff_gexp <- merged_df_ord$Mouse_NUM_GENES-merged_df_ord$Human_NUM_GENES
corner(merged_df_ord)

pdf("Seperation.pdf")
g1 <- ggplot(merged_df_ord,aes(x=Mouse_NUM_GENES,y=Human_NUM_GENES,fill=Cluster,col=Cluster))+
  geom_point()+theme_cowplot()+geom_abline(show.legend = "Y=X")+
  scale_color_manual(values = color_list[["Cluster"]], name = "Cluster",na.value = "gray")

g2 <- ggplot(merged_df_ord,aes(x=Mouse_NUM_TRANSCR,y=Human_NUM_TRANSCR,fill=Cluster,col=Cluster))+
  geom_point()+theme_cowplot()+geom_abline(show.legend = "Y=X")+
  scale_color_manual(values = color_list[["Cluster"]], name = "Cluster",na.value = "gray")

g3 <- ggplot(merged_df_ord,aes(x=Cluster,y=Diff_tr))+
  geom_violin()+geom_jitter(aes(fill=Cluster,col=Cluster))+
  scale_color_manual(values = color_list[["Cluster"]], name = "Cluster",na.value = "gray")+
  geom_abline(x=0)+theme_cowplot()

g4 <- ggplot(merged_df_ord,aes(x=Cluster,y=Diff_gexp))+
  geom_violin()+geom_jitter(aes(fill=Cluster,col=Cluster))+
  scale_color_manual(values = color_list[["Cluster"]], name = "Cluster",na.value = "gray")+
  geom_abline(x=0)+theme_cowplot()
g1
g2
g3
g4
DimPlot(test,group.by=c("condition"))
DimPlot(test,group.by=c("Cluster"))
dev.off()

pdf("Sep_2.pdf",width=15,height=10)
ggarrange(g1, g2, g3,g4,
                    labels = c("A", "B", "C","D"),
                    ncol = 2, nrow = 2)
dev.off()

pdf("Sep_3.pdf",width=15)
g1 <- DimPlot(test,group.by=c("condition"),cols=color_cond)
g2<- DimPlot(test,group.by=c("Cluster"),cols=color_clust)
g3 <- ggplot(merged_df_ord,aes(x=Cluster,y=Diff_gexp))+
  geom_violin()+geom_jitter(aes(fill=Cluster,col=Cluster))+
  scale_color_manual(values = color_list[["Cluster"]], name = "Cluster",na.value = "gray")+
  geom_abline(x=0)+theme_cowplot()
ggarrange(g1, g2, g3,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
dev.off()



pdf("Dim_wrap.pdf",width=15,height=15)
DimPlot(test,group.by=c("condition"),cols=color_cond,split.by="condition",ncol=3)
dev.off()

pdf("Sep2.pdf")
hist(merged_df_ord$Mouse_NUM_TRANSCR-merged_df_ord$Human_NUM_TRANSCR,breaks = 500)


hist(colSums(test@assays$RNA@counts==0),breaks = 200)
table(colSums(test@assays$RNA@counts==0) > 15000)


cells_rmv <- colnames(test)[colSums(test@assays$RNA@counts==0) > 15900]
DimPlot(test)
subset_com <- subset(test,cells=cells_rmv)
DimPlot(subset_com)
dev.off()






#mouse_data <- subset(test,subset=Cluster%in%c(1,2,3,4,5,6,7,9,10,13,14,15,16,17,18,20))
#human_data <- subset(test,subset=Cluster%in%c(8,11,12,19))



#mouse_data <- subset(test,subset=Cluster%in%c(1,2,3,4,5,6,8,10,13,14,15,16,17,18,20))
#human_data <- subset(test,subset=Cluster%in%c(7,9))

mouse_data <- subset(test,subset=Cluster%in%c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17))
human_data <- subset(test,subset=Cluster%in%c(2))



saveRDS(mouse_data,"/home/users/dkyriakis/PhD/Projects/Yahaya/5.Seperation/Mouse.rds")
saveRDS(human_data,"/home/users/dkyriakis/PhD/Projects/Yahaya/5.Seperation/Human.rds")
