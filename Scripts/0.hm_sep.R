---
title: "GBM_TME"
author: "Dimitrios Kyriakis"
date: "26/6/2020"
---

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
backup_obj<-test
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


hist(colSums(Combined@assays$RNA@counts==0),breaks = 200)
table(colSums(Combined@assays$RNA@counts==0) > 15000)


cells_rmv <- colnames(Combined)[colSums(Combined@assays$RNA@counts==0) > 15900]
DimPlot(Combined)
subset_com <- subset(Combined,cells=cells_rmv)
DimPlot(subset_com)
dev.off()






mouse_data <- subset(test,subset=Cluster%in%c(1,2,3,4,5,6,7,9,10,13,14,15,16,17,18,20))
human_data <- subset(test,subset=Cluster%in%c(8,11,12,19))



dir.create("6.Projection_sep")
setwd("6.Projection_sep")
DefaultAssay(mouse_data) <-"RNA"
mouse_data2 <-liger_projection(mouse_data,split.by="condition",k=20,lamda=5,n_neighbors=10)
mouse_data2$Cluster <- mouse_data2$Liger_Clusters

pdf("Markers.pdf")
FeaturePlot(mouse_data2, c("PLP1","P2RY12","GJA1","OLIG1","GAD2","IGF2","CCDC153"),
    order=TRUE,cols = c("lightgrey","#FDBB84","#EF6548","#D7301F","#B30000","#7F0000"),reduction="liger_umap")
dev.off()

pdf("Cells.pdf")
DimPlot(mouse_data2,group.by = "condition",cols=color_cond,reduction="liger_umap")
DimPlot(mouse_data2,group.by = "Cluster",cols=color_clust,reduction="liger_umap")
DimPlot(mouse_data2,group.by = "Phase",reduction="liger_umap")
dev.off()


dat <- data.frame(table(mouse_data2$Cluster,mouse_data2$Condition))
names(dat) <- c("Cluster","Condition","Count")

pdf("Barplot.pdf",width=12)
ggplot(data=dat, aes(x=Cluster, y=Count, fill=Condition)) + geom_bar(stat="identity")+theme_cowplot()+ scale_fill_manual(values=color_cond)
dev.off()

setwd("../")






dir.create("7.Mouse")
setwd("7.Mouse")
Seurat <- mouse_data2
Idents(Seurat)<-Seurat$Cluster
DefaultAssay(Seurat) <- "RNA"
Seurat <-NormalizeData(Seurat)
Seurat <-ScaleData(Seurat)
markers <- FindAllMarkers(Seurat,test.use="MAST",only.pos =T,latent.vars="nCount_RNA",logfc.threshold = 0.1)

markers$gene <- rownames(markers)

pdf("DoHeatmap.pdf",width=12,height=10)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = abs(avg_logFC))
DoHeatmap(object = Seurat,features = unique(top10$gene),raster = F)+ 
    theme(text = element_text(size = 6))
top3 <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = abs(avg_logFC))
DotPlot(Seurat, features = unique(top3$gene),  dot.scale = 8) + RotatedAxis()+ 
    theme(text = element_text(size = 9))
dev.off()
write.table(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.1,],"sig_avatar.txt")


