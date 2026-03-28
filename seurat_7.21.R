


########################human data############################
library(dplyr)
library(Seurat)
library(patchwork)
memory.limit(10000000)
######################load data and name cluster#################

setwd(".....\\cardiomyocytes\\science\\allcell")

heart <- readRDS(".....\\cardiomyocytes\\science\\allcell\\humanheart.rds")
heart=UpdateSeuratObject(heart)

table(heart@active.ident)

H07 <- subset(heart, subset = heart@meta.data[[,"donor_id"]] == 'D1') 
H07 <- subset(x = heart, subset = Patient == "H07")
H07 <- NormalizeData(H07, normalization.method = "LogNormalize", scale.factor = 10000)
H07 <- NormalizeData(H07)
H07 <- FindVariableFeatures(H07, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(H07)
H07 <- ScaleData(H07, features = all.genes)
H07$type<-"PKP2"

H09 <- subset(x = heart, subset = Patient == "H09")
H09 <- NormalizeData(H09, normalization.method = "LogNormalize", scale.factor = 10000)
H09 <- NormalizeData(H09)
H09 <- FindVariableFeatures(H09, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(H09)
H09 <- ScaleData(H09, features = all.genes)
H09$type<-"PKP2"

H11 <- subset(x = heart, subset = Patient == "H11")
H11 <- NormalizeData(H11, normalization.method = "LogNormalize", scale.factor = 10000)
H11 <- NormalizeData(H11)
H11 <- FindVariableFeatures(H11, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(H11)
H11 <- ScaleData(H11, features = all.genes)
H11$type<-"PKP2"

H42 <- subset(x = heart, subset = Patient == "H42")
H42 <- NormalizeData(H42, normalization.method = "LogNormalize", scale.factor = 10000)
H42 <- NormalizeData(H42)
H42 <- FindVariableFeatures(H42, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(H42)
H42 <- ScaleData(H42, features = all.genes)
H42$type<-"PKP2"

D2 <- subset(x = heart, subset = Patient == "D2")
D2 <- NormalizeData(D2, normalization.method = "LogNormalize", scale.factor = 10000)
D2 <- NormalizeData(D2)
D2 <- FindVariableFeatures(D2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(D2)
D2 <- ScaleData(D2, features = all.genes)
D2$type<-"Normal"

D7 <- subset(x = heart, subset = Patient == "D7")
H07 <- NormalizeData(H07, normalization.method = "LogNormalize", scale.factor = 10000)
H07 <- NormalizeData(H07)
H07 <- FindVariableFeatures(H07, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(H07)
H07 <- ScaleData(H07, features = all.genes)
D7$type<-"Normal"

H51 <- subset(x = heart, subset = Patient == "H51")
H51 <- NormalizeData(H51, normalization.method = "LogNormalize", scale.factor = 10000)
H51 <- NormalizeData(H51)
H51 <- FindVariableFeatures(H51, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(H51)
H51 <- ScaleData(H51, features = all.genes)
H51$type<-"Normal"

H53 <- subset(x = heart, subset = Patient == "H53")
H53 <- NormalizeData(H53, normalization.method = "LogNormalize", scale.factor = 10000)
H53 <- NormalizeData(H53)
H53 <- FindVariableFeatures(H53, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(H53)
H53 <- ScaleData(H53, features = all.genes)
H53$type<-"Normal"

pbmc <- merge(x = H51, y = list(H53, D2, D7, H07, H09, H11, H42))

# split the dataset into a list of seurat objects
pbmc.list <- SplitObject(pbmc, split.by = "donor_id")

# normalize and identify variable features for each dataset independently
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = pbmc.list)

heartpkp2 <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features,dims = 1:20,k.anchor = 5,k.filter = 30)

# this command creates an 'integrated' data assay
heart <- IntegrateData(anchorset = heartpkp2)
###########################################################################
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(heart) <- "RNA"

# Run the standard workflow for visualization and clustering
heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "^MT-")
VlnPlot(heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "donor_id", ncol = 3)
heart <- subset(heart, subset = nFeature_RNA > 250 & nFeature_RNA < 7500 & nCount_RNA < 75000&percent.mt<25)

DefaultAssay(heart) <- "integrated"

heart <- ScaleData(heart, verbose = FALSE)
heart <- RunPCA(heart, npcs = 20, verbose = FALSE)
heart <- RunUMAP(heart, reduction = "pca", dims = 1:20)
heart <- FindNeighbors(heart, reduction = "pca", dims = 1:20)
heart <- FindClusters(heart, resolution = 0.2)
#######################################################################################################
my36colors <- c('#FF3366', '#33FF33', '#FF9900','#339933', '#FF33FF', '#3399FF', '#FF6666', '#009999',
                '#33FFFF', '#E59CC4', '#AB3282', '#968175', '#FF0000', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
##############################################################################################
new.cluster.ids <- c("Cardiomyocyte", "Fibroblast","Pericyte","Endothelial","Immune","Cardiomyocyte","Neuronal")
names(new.cluster.ids) <- levels(heart)
heart <- RenameIdents(heart, new.cluster.ids)

DimPlot(heart, label=T,label.size=6, reduction = "umap",raster=FALSE)+ NoLegend()
DimPlot(heart, group.by="Patient",label=T,reduction = "umap")
DimPlot(heart, group.by="type",label=F,pt.size=0.1,reduction = "umap",raster=FALSE,cols=my36colors)
VlnPlot(heart, features = c("PROM1"))
FeaturePlot(heart, features = c("MYH6"))

############################subset cardiomyocyte and find DEG####################################33
Cardiomyocytes <- subset(heart,idents = c("Cardiomyocyte"))

plot1 <- FeatureScatter(cardiomyocytes, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cardiomyocytes, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

cardiomyocytes <- subset(cardiomyocytes, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 30000 & percent.mt < 20)
cardiomyocytes <- NormalizeData(cardiomyocytes, normalization.method = "LogNormalize", scale.factor = 10000)
cardiomyocytes <- NormalizeData(cardiomyocytes)
cardiomyocytes <- FindVariableFeatures(cardiomyocytes, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(cardiomyocytes), 10)

plot1 <- VariableFeaturePlot(cardiomyocytes)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
cardiomyocytes@meta.data
all.genes <- rownames(cardiomyocytes)
cardiomyocytes <- ScaleData(cardiomyocytes, features = all.genes)
cardiomyocytes <- RunPCA(cardiomyocytes, features = VariableFeatures(object = cardiomyocytes))

print(cardiomyocytes[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(cardiomyocytes, reduction = "pca")
DimHeatmap(cardiomyocytes, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(cardiomyocytes, dims = 1:15, cells = 500, balanced = TRUE)
cardiomyocytes <- JackStraw(cardiomyocytes, num.replicate = 100)
cardiomyocytes <- ScoreJackStraw(cardiomyocytes, dims = 1:20)

JackStrawPlot(cardiomyocytes, dims = 1:15)
ElbowPlot(cardiomyocytes)

cardiomyocytes <- FindNeighbors(cardiomyocytes, dims = 1:18)
cardiomyocytes <- FindClusters(cardiomyocytes, resolution = 0.06)
cardiomyocytes <- RunUMAP(cardiomyocytes, dims = 1:18)

DimPlot(cardiomyocytes, label=T,label.size=5, reduction = "umap",cols=my36colors,raster=FALSE)+ NoLegend()
DimPlot(cardiomyocytes, group.by="Primary.Genetic.Diagnosis",label=F,reduction = "umap",cols=my36colors)

cardiomyocytes$cell_type.Primary.Genetic.Diagnosis <- paste(Idents(cardiomyocytes), cardiomyocytes$Primary.Genetic.Diagnosis, sep = "_")
cardiomyocytes$cell_type <- Idents(cardiomyocytes)
Idents(cardiomyocytes) <- "cell_type.Primary.Genetic.Diagnosis"

b.interferon.response <- FindMarkers(Cardiomyocytes, ident.1 = "Cardio_PKP2", ident.2 = "Cardio_Normal", verbose = FALSE)
write.table(b.interferon.response,sep="\t", row.names = T, file="pkp2vsnormal.xlsX")

##########################################################vanco plot######################
library(ggplot2)
library(ggrepel)
library(dplyr)

res<- read.table(file="pkp2vsnormal.xlsX",header=T)
head(res,n=10)
fix(res)    ## add "gene" group 
## 1st. add labels 
res$labels <- ifelse(abs(res$avg_log2FC) > 0.25 & res$p_val_adj < 0.05, "Both", 
                     ifelse(res$p_val_adj < 0.05, "Significant",
                            ifelse(abs(res$avg_log2FC) > 1, "log2FC > 1", "None")))

table(res$labels)

p <- ggplot(res, aes(avg_log2FC, -1*log10(p_val_adj)))
p <- p + geom_point(aes(color=labels), alpha=0.6,size=1.2) + xlim(c(-2.5, 2.5)) + ylim(c(0, 50)) # draw the basic structure of volcano plot
plot_grid(p)

### change color 
p <- p + scale_color_manual(values =c("red", "green", "grey")) 

## stage 1: #7CAE00  #00BFC4
## stage 2; #C59D07
## stage 3; #F8776F
## stage 4; #FB61D7
p + geom_text_repel(data = geneset1_left, aes(label = gene), size = 4,
                    min.segment.length = unit(0, 'lines'), nudge_y = 1,  ## label line length
                    direction= "both")  + 
  geom_text_repel(data = geneset1_right, aes(label = gene), size = 4, 
                  colour="#FB61D7", min.segment.length = unit(0, 'lines'), nudge_y = 1,  ## label line length
                  direction= "both") 

## All DEGS
mygenes<-subset(res, labels == "Both")
write.csv(mygenes, file="DEGs.CSV")
p + geom_text_repel(data = mygenes, aes(label = gene), size = 3)    ## visulalize with all DEGs

# set a and y axis limits
p <- p + xlim(c(-100, 100)) + ylim(c(0, 100))

# add title and axis titles
p <- p + labs(title="Volcano plot", x=expression(log[2](Fold_Change)), y=expression(-log[10](p_adjusted_value))) 

# set background
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 



##########################################################Mouse data
library(SeuratDisk)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)
memory.limit(10000000)
setwd("...\\heart\\cellreport")

##############################################################load data and name cluster######################
pbmc2.data <- Read10X(data.dir = "D:\\Research field\\zebrafish\\heart\\cellreport\\GSE153481_RAW\\p1_1sham")
pbmc2 <- CreateSeuratObject(counts = pbmc2.data, project = "p1_1sham", min.cells = 3, min.features = 200)
pbmc2 <- NormalizeData(pbmc2, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc2 <- NormalizeData(pbmc2)
pbmc2 <- FindVariableFeatures(pbmc2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc2)
pbmc2 <- ScaleData(pbmc2, features = all.genes)
pbmc2$type<-"sham"
pbmc2$sample<-"P1"

pbmc4.data <- Read10X(data.dir = "D:\\Research field\\zebrafish\\heart\\cellreport\\GSE153481_RAW\\p1_3sham")
pbmc4 <- CreateSeuratObject(counts = pbmc4.data, project = "p1_3sham", min.cells = 3, min.features = 200)
pbmc4 <- NormalizeData(pbmc4, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc4 <- NormalizeData(pbmc4)
pbmc4 <- FindVariableFeatures(pbmc4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc4)
pbmc4 <- ScaleData(pbmc4, features = all.genes)
pbmc4$type<-"sham"
pbmc4$sample<-"P1"

pbmc6.data <- Read10X(data.dir = "D:\\Research field\\zebrafish\\heart\\cellreport\\GSE153481_RAW\\p8_1sham")
pbmc6 <- CreateSeuratObject(counts = pbmc6.data, project = "p8_1sham", min.cells = 3, min.features = 200)
pbmc6 <- NormalizeData(pbmc6, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc6 <- NormalizeData(pbmc6)
pbmc6 <- FindVariableFeatures(pbmc6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc6)
pbmc6 <- ScaleData(pbmc6, features = all.genes)
pbmc6$type<-"sham"
pbmc6$sample<-"P8"

pbmc8.data <- Read10X(data.dir = "D:\\Research field\\zebrafish\\heart\\cellreport\\GSE153481_RAW\\p8_3sham")
pbmc8 <- CreateSeuratObject(counts = pbmc8.data, project = "p8_3sham", min.cells = 3, min.features = 200)
pbmc8 <- NormalizeData(pbmc8, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc8 <- NormalizeData(pbmc8)
pbmc8 <- FindVariableFeatures(pbmc8, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc8)
pbmc8 <- ScaleData(pbmc8, features = all.genes)
pbmc8$type<-"sham"
pbmc8$sample<-"P8"

pbmc9.data <- Read10X(data.dir = "D:\\Research field\\zebrafish\\heart\\cellreport\\GSE153481_RAW\\E18")
pbmc9 <- CreateSeuratObject(counts = pbmc9.data, project = "E18", min.cells = 3, min.features = 200)
pbmc9 <- NormalizeData(pbmc9, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc9 <- NormalizeData(pbmc9)
pbmc9 <- FindVariableFeatures(pbmc9, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc9)
pbmc9 <- ScaleData(pbmc9, features = all.genes)
pbmc9$type<-"sham"
pbmc9$sample<-"E18"

pbmc10.data <- Read10X(data.dir = "D:\\Research field\\zebrafish\\heart\\cellreport\\GSE153481_RAW\\ADULT")
pbmc10 <- CreateSeuratObject(counts = pbmc10.data, project = "Adult", min.cells = 3, min.features = 200)
pbmc10 <- NormalizeData(pbmc10, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc10 <- NormalizeData(pbmc10)
pbmc10 <- FindVariableFeatures(pbmc10, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc10)
pbmc10 <- ScaleData(pbmc10, features = all.genes)
pbmc10$type<-"sham"
pbmc10$sample<-"Adult"

pbmc <- merge(x = pbmc2, y = list(pbmc4, pbmc6, pbmc8, pbmc9, pbmc10))

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 7000  & percent.mt < 20)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


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

pbmc <- FindNeighbors(pbmc, dims = 1:35)
pbmc <- FindClusters(pbmc, resolution = 0.1)
pbmc <- RunUMAP(pbmc, dims = 1:35)

DimPlot(pbmc, label=T,label.size=4,reduction = "umap")+ NoLegend()
DimPlot(pbmc, group.by="sample",label=T,reduction = "umap")
DimPlot(pbmc, group.by="type",label=T,reduction = "umap")

library(DropletUtils)
write10xCounts(x = pbmc@assays$RNA@counts, path = '10x', version="3")
active.ident <- pbmc@active.ident
alldata <- AddMetaData(pbmc, active.ident)
Idents(object = pbmc)
pbmc[["celltype"]] <- Idents(object = pbmc)
write.csv(pbmc@meta.data,'main.csv')

FeaturePlot(pbmc, features = c("Ccl2"))

new.cluster.ids <- c("Endothelial", "Fibroblast","Pericyte","Macrophage","Endothelial","Fibroblast2","Cardiomyocyte","Epicardial cell","Granulocyte","Immune cell","Endothelial")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
saveRDS(pbmc, file = "pbmc.rds")

###############################################################################zebrafish data
library(SeuratDisk)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)

memory.limit(10000000)
setwd("...\\heart\\NG")
#########################################################################load data and name cluster
pbmc1.data <- Read10X_h5("E:\\Research field\\zebrafish\\heart\\NG\\ctr1.h5")
pbmc1 <- CreateSeuratObject(counts = pbmc1.data, project = "ctr1", min.cells = 3, min.features = 200)
pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc1 <- NormalizeData(pbmc1)
pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc1)
pbmc1 <- ScaleData(pbmc1, features = all.genes)

pbmc2.data <- Read10X_h5("E:\\Research field\\zebrafish\\heart\\NG\\ctr2.h5")
pbmc2 <- CreateSeuratObject(counts = pbmc2.data, project = "ctr2", min.cells = 3, min.features = 200)
pbmc2 <- NormalizeData(pbmc2, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc2 <- NormalizeData(pbmc2)
pbmc2 <- FindVariableFeatures(pbmc2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc2)
pbmc2 <- ScaleData(pbmc2, features = all.genes)

pbmc3.data <- Read10X_h5("E:\\Research field\\zebrafish\\heart\\NG\\ctr3.h5")
pbmc3 <- CreateSeuratObject(counts = pbmc3.data, project = "ctr3", min.cells = 3, min.features = 200)
pbmc3 <- NormalizeData(pbmc3, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc3 <- NormalizeData(pbmc3)
pbmc3 <- FindVariableFeatures(pbmc3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc3)
pbmc3 <- ScaleData(pbmc3, features = all.genes)

pbmc4.data <- Read10X_h5("D:\\Research field\\zebrafish\\heart\\NG\\3n1.h5")
pbmc4 <- CreateSeuratObject(counts = pbmc4.data, project = "3n1", min.cells = 3, min.features = 200)
pbmc4 <- NormalizeData(pbmc4, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc4 <- NormalizeData(pbmc4)
pbmc4 <- FindVariableFeatures(pbmc4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc4)
pbmc4 <- ScaleData(pbmc4, features = all.genes)

pbmc5.data <- Read10X_h5("D:\\Research field\\zebrafish\\heart\\NG\\3n2.h5")
pbmc5 <- CreateSeuratObject(counts = pbmc5.data, project = "3n2", min.cells = 3, min.features = 200)
pbmc5 <- NormalizeData(pbmc5, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc5 <- NormalizeData(pbmc5)
pbmc5 <- FindVariableFeatures(pbmc5, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc5)
pbmc5 <- ScaleData(pbmc5, features = all.genes)

pbmc6.data <- Read10X_h5("D:\\Research field\\zebrafish\\heart\\NG\\3n3.h5")
pbmc6 <- CreateSeuratObject(counts = pbmc6.data, project = "3n3", min.cells = 3, min.features = 200)
pbmc6 <- NormalizeData(pbmc6, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc6 <- NormalizeData(pbmc6)
pbmc6 <- FindVariableFeatures(pbmc6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc6)
pbmc6 <- ScaleData(pbmc6, features = all.genes)

pbmc7.data <- Read10X_h5("D:\\Research field\\zebrafish\\heart\\NG\\7n1.h5")
pbmc7 <- CreateSeuratObject(counts = pbmc7.data, project = "7n1", min.cells = 3, min.features = 200)
pbmc7 <- NormalizeData(pbmc7, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc7 <- NormalizeData(pbmc7)
pbmc7 <- FindVariableFeatures(pbmc7, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc7)
pbmc7 <- ScaleData(pbmc7, features = all.genes)

pbmc8.data <- Read10X_h5("D:\\Research field\\zebrafish\\heart\\NG\\7n2.h5")
pbmc8 <- CreateSeuratObject(counts = pbmc8.data, project = "7n2", min.cells = 3, min.features = 200)
pbmc8 <- NormalizeData(pbmc8, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc8 <- NormalizeData(pbmc8)
pbmc8 <- FindVariableFeatures(pbmc8, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc8)
pbmc8 <- ScaleData(pbmc8, features = all.genes)
#################################################################################################################
pbmc <- merge(x = pbmc1, y = list(pbmc2, pbmc3,pbmc4, pbmc5,pbmc6, pbmc7,pbmc8))
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset =nFeature_RNA < 5000  & nCount_RNA < 50000)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

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
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:20)


new.cluster.ids <- c("Cardiomyocyte", "Erythrocyte","Cardiomyocyte","Endocardium","Endocardium","Endocardium",
                     "Fibroblasts","Smootn muscle","immune cell","immune cell","Epicardium","Endothelial","Epicardium","immune cell","immune cell","Endothelial",
                     "Neutrophils","Erythrocyte")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
############################################################################
my36colors <- c('#FF3366', '#33FF33', '#FF9900','#339933', '#FF33FF', '#3399FF', '#FF6666', '#009999',
                '#33FFFF', '#E59CC4', '#AB3282', '#CC9090', '#FF0000', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
#########################################################################
DimPlot(pbmc, label=T,label.size=4,reduction = "umap")+ NoLegend()
DimPlot(pbmc, label=T,label.size=6,reduction = "umap",cols=my36colors)+ NoLegend()
DimPlot(pbmc, group.by="sample",label=T,reduction = "umap")
FeaturePlot(pbmc, features = c("pkp2"))

saveRDS(pbmc, file = "zebrafishheart.rds")
####################################################################################
