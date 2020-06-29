---
title: "GBM_TME"
author: "Dimitrios Kyriakis"
date: "26/6/2020"
---


WORKDIR<-"/home/users/dkyriakis/PhD/Projects/Yahaya/"

human_dir <- paste0(WORKDIR,"/human/")
mouse_dir <- paste0(WORKDIR,"/mouse/")

list.files(human_dir,full.names = TRUE)[c(4,6,8,1,2)]


iter=1
human_file <- NULL
for (file in list.files(human_dir,full.names = TRUE)[c(4,6,8,1,2)]){
  print(iter)
  print(file)
  temp_file1 <- read.table(file,header = TRUE)
  if(dim(temp_file1)[2]==4){
    temp_file1[2]<-NULL
  }
  temp_file1[,1] <- paste0(temp_file1[,1],"_",iter) 
  print(dim(temp_file1))
  human_file <- rbind(human_file,temp_file1)
  iter=iter+1
}

iter=1
mouse_file <- NULL
for (file in list.files(mouse_dir,full.names = TRUE)[c(4,6,8,1,2)]){
  temp_file1 <- read.table(file,header = TRUE)
  if(dim(temp_file1)[2]==4){
    temp_file1[2]<-NULL
  }
  temp_file1[,1] <- paste0(temp_file1[,1],"_",iter) 
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

pdf("Cluster_4_Validation.pdf")
ggplot(merged_df_ord,aes(x=Mouse_NUM_GENES,y=Human_NUM_GENES,fill=Cluster,col=Cluster))+
  geom_point()+theme_cowplot()+geom_abline(show.legend = "Y=X")+
  scale_color_manual(values = color_list[["Cluster"]], name = "Cluster",na.value = "gray")

ggplot(merged_df_ord,aes(x=Mouse_NUM_TRANSCR,y=Human_NUM_TRANSCR,fill=Cluster,col=Cluster))+
  geom_point()+theme_cowplot()+geom_abline(show.legend = "Y=X")+
  scale_color_manual(values = color_list[["Cluster"]], name = "Cluster",na.value = "gray")

ggplot(merged_df_ord,aes(x=Cluster,y=Diff_tr))+
  geom_violin()+geom_jitter(aes(fill=Cluster,col=Cluster))+
  scale_color_manual(values = color_list[["Cluster"]], name = "Cluster",na.value = "gray")+
  geom_abline(x=0)+theme_cowplot()

ggplot(merged_df_ord,aes(x=Cluster,y=Diff_gexp))+
  geom_violin()+geom_jitter(aes(fill=Cluster,col=Cluster))+
  scale_color_manual(values = color_list[["Cluster"]], name = "Cluster",na.value = "gray")+
  geom_abline(x=0)+theme_cowplot()

dev.off()



hist(merged_df_ord$Mouse_NUM_TRANSCR-merged_df_ord$Human_NUM_TRANSCR,breaks = 500)


hist(colSums(Combined@assays$RNA@counts==0),breaks = 200)
table(colSums(Combined@assays$RNA@counts==0) > 15000)


cells_rmv <- colnames(Combined)[colSums(Combined@assays$RNA@counts==0) > 15900]
DimPlot(Combined)
subset_com <- subset(Combined,cells=cells_rmv)
DimPlot(subset_com)
