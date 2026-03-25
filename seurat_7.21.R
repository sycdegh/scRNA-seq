
library(Seurat)
library(ggplot2)
args=commandArgs(T)

#5个10X数据的路径
Read10X_dir_M=args[1]
Read10X_dir_C=args[2]
Read10X_dir_T1=args[3]
Read10X_dir_T2=args[4]
Read10X_dir_T3=args[5]
out_dir=args[6]
setwd(out_dir)

###########################################################################
###########################################################load data
M.data <- Read10X(data.dir = Read10X_dir_M,gene.column = 2)
M<- CreateSeuratObject(counts = M.data, project = "M",min.cells = 3, min.features = 200)
M$group<-"M"
M$time<-"1week"
M$sample<-"treat"
M$gene<-"Td4F-c"

C.data <- Read10X(data.dir = Read10X_dir_C,gene.column = 2)
C<- CreateSeuratObject(counts = C.data, project = "C",min.cells = 3, min.features = 200)
C$group<-"C"
C$time<-"1week"
C$sample<-"control"
C$gene<-"None"


T1.data <- Read10X(data.dir = Read10X_dir_T1,gene.column = 2)
T1<- CreateSeuratObject(counts = T1.data, project = "T1",min.cells = 3, min.features = 200)
T1$group<-"T1"
T1$time<-"1week"
T1$sample<-"treat"
T1$gene<-"OSKM"

T2.data <- Read10X(data.dir = Read10X_dir_T2,gene.column = 2)
T2<- CreateSeuratObject(counts = T2.data, project = "T2",min.cells = 3, min.features = 200)
T2$time<-"1week"
T2$group<-"T2"
T2$sample<-"treat"
T2$gene<-"Td4F"

T3.data <- Read10X(data.dir = Read10X_dir_T3,gene.column = 2)
T3<- CreateSeuratObject(counts = T3.data, project = "T3",min.cells = 3, min.features = 200)
T3$time<-"1week"
T3$group<-"T3"
T3$sample<-"treat"
T3$gene<-"OSKM+Td4F"

W2C.data <- Read10X(data.dir = Read10X_dir_W2T2,gene.column = 2)
W2C<- CreateSeuratObject(counts = W2C.data, project = "W2C",min.cells = 3, min.features = 200)
W2C$group<-"W2C"
W2C$time<-"2week"
W2C$sample<-"Control"
W2C$gene<-"None"

W2T1.data <- Read10X(data.dir = Read10X_dir_W2T1a,gene.column = 2)
W2T1<- CreateSeuratObject(counts = W2T1.data, project = "W2T1",min.cells = 3, min.features = 200)
W2T1$group<-"W2T1"
W2T1$time<-"2week"
W2T1$sample<-"treat"
W2T1$gene<-"OSKM"


W2T2.data <- Read10X(data.dir = Read10X_dir_W2T2a,gene.column = 2)
W2T2<- CreateSeuratObject(counts = W2T2.data, project = "W2T2",min.cells = 3, min.features = 200)
W2T2$group<-"W2T2"
W2T2$time<-"2week"
W2T2$sample<-"treat"
W2T2$gene<-"Td4F"


W2T3.data <- Read10X(data.dir = Read10X_dir_W2T3a,gene.column = 2)
W2T3<- CreateSeuratObject(counts = W2T3.data, project = "W2T3",min.cells = 3, min.features = 200)
W2T3$group<-"W2T3"
W2T3$time<-"2week"
W2T3$sample<-"treat"
W2T3$gene<-"OSKM+Td4F"

pbmc <- merge(x = C, y = list(M,T1,T2,T3,W2C,W2T1,W2T2,W2T3))
######################################################################################################################
##############################################################################process
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

pdf("vln.pdf",width = 10,height = 8)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & nCount_RNA < 60000 & percent.mt < 25)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
pbmc@meta.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.4)
pbmc <- RunUMAP(pbmc, dims = 1:20)

###############################################################################################
pdf("DimPlot_10_8-unname.pdf",width = 10,height = 8)
DimPlot(pbmc, label=T,pt.size=0.6,label.size=6,reduction = "umap")
dev.off()

#################################################################################################rename
new.cluster.ids <- c("Tendon cell","Macrophage", "Endothelial","Proliferative Teno","Pericyte/SMC","Tendon cell",
                     "T/NK cell","Tendon cell","Neutrophil","Endothelial","Endothelial","Macrophage",
                     "Endothelial","Tendon cell","Tendon cell")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
#################################################################################################plot
pdf("DimPlot_10_8name.pdf",width = 10,height = 8)
DimPlot(pbmc, label=T,label.size=6,reduction = "umap")+NoLegend()
dev.off()

pdf("DimPlot_group_10_8.pdf",width = 10,height = 8)
DimPlot(pbmc, group.by="orig.ident",label=T,reduction = "umap")
dev.off()

pdf("DimPlot_group_gene.pdf",width = 10,height = 8)
DimPlot(pbmc, group.by="time",label=T,reduction = "umap")
dev.off()

pdf("DimPlot_split_20_4.pdf",width = 20,height = 4)
DimPlot(pbmc, split.by="orig.ident",label=T,reduction = "umap")
dev.off()

pdf("DimPlot_split_gene.pdf",width = 20,height = 4)
DimPlot(pbmc, split.by="time",label=T,reduction = "umap")
dev.off()

########################################################################findmarker

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, test.use = "wilcox",min.pct = 0.25, logfc.threshold = 0.25)

write.table(pbmc.markers,"pbmc-name.xlsx",sep="\t", row.names = T,quote=F)

#################################################################################################heatmap
library(ggplot2)
library(dplyr)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

pdf("top4-heatmap.pdf",width = 10,height = 8)
DoHeatmap(pbmc, features = top4$gene) + NoLegend()
dev.off()
##############################################################################subset
teno <- subset(pbmc,idents = c("Tendon cell","Proliferative Teno"))
saveRDS(teno,file = "teno.rds")

immune <- subset(pbmc,idents = c("Macrophage","Neutrophil","T/NK cell"))
saveRDS(immune,file = "immune.rds")







