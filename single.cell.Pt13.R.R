rm(list = ls())
library(dbplyr)
library(Seurat)
library(patchwork)
getwd()
Pt13_a<-Matrix::readMM("/home/shpc_101103/liver/GSE112271/Pt13.a/GSM3064818_Pt13.a_matrix.mtx.gz")
Pt13_a[1:4,1:4]
Pt13_a1<-read.table("/home/shpc_101103/liver/GSE112271/Pt13.a/GSM3064818_Pt13.a_barcodes.tsv.gz")
Pt13_a2<-read.table("/home/shpc_101103/liver/GSE112271/Pt13.a/GSM3064818_Pt13.a_genes.tsv.gz")
colnames(Pt13_a)=Pt13_a1$V1
rownames(Pt13_a)=Pt13_a2$V2
Pt13_a[1:4,1:4]
scRNA_Pt13_a= CreateSeuratObject(counts = Pt13_a,project="Pt13_a",
                           min.cell=3,
                           min.features=200)
Pt13_b<-Matrix::readMM("/home/shpc_101103/liver/GSE112271/Pt13.b/GSM3064819_Pt13.b_matrix.mtx.gz")
Pt13_b[1:4,1:4]
Pt13_b1<-read.table("/home/shpc_101103/liver/GSE112271/Pt13.b/GSM3064819_Pt13.b_barcodes.tsv.gz")
Pt13_b2<-read.table("/home/shpc_101103/liver/GSE112271/Pt13.b/GSM3064819_Pt13.b_genes.tsv.gz")
colnames(Pt13_b)=Pt13_b1$V1
rownames(Pt13_b)=Pt13_b2$V2
Pt13_b[1:4,1:4]
scRNA_Pt13_b= CreateSeuratObject(counts = Pt13_b,project="Pt13_b",
                                 min.cell=3,
                                 min.features=200)

Pt13_c<-Matrix::readMM("/home/shpc_101103/liver/GSE112271/Pt13.c/GSM3064820_Pt13.c_matrix.mtx.gz")
Pt13_c[1:4,1:4]
Pt13_c1<-read.table("/home/shpc_101103/liver/GSE112271/Pt13.c/GSM3064820_Pt13.c_barcodes.tsv.gz")
Pt13_c2<-read.table("/home/shpc_101103/liver/GSE112271/Pt13.c/GSM3064820_Pt13.c_genes.tsv.gz")
colnames(Pt13_c)=Pt13_c1$V1
rownames(Pt13_c)=Pt13_c2$V2
Pt13_c[1:4,1:4]
scRNA_Pt13_c= CreateSeuratObject(counts = Pt13_c,project="Pt13_c",
                                 min.cell=3,
                                 min.features=200)
scRNA=merge(scRNA_Pt13_a,y=c(scRNA_Pt13_b,scRNA_Pt13_c))
head(scRNA@meta.data)
summary(scRNA@meta.data)
scRNA@assays$RNA@counts[1:4,1:4]
dim(scRNA)

table(grepl("^MT-",rownames(scRNA)))
scRNA[["percent.mt"]]<-PercentageFeatureSet(scRNA,pattern="^MT-")
head(scRNA@meta.data)
summary(scRNA@meta.data)
pctMT=10 
scRNA<-subset(scRNA,subset=percent.mt < pctMT)
dim(scRNA)
table(grepl("^ERCC-",rownames(scRNA)))

scRNA[["percent.ERCC"]]<-PercentageFeatureSet(scRNA,pattern="^ERCC-")
head(scRNA@meta.data)
summary(scRNA@meta.data)
rownames(scRNA)[grep("^ERCC-",rownames(scRNA))]
#sum(scRNA$percent.ERCC<40)

sum(scRNA$percent.ERCC<10)
pctERCC=10
scRNA<-subset(scRNA,subset=percent.ERCC<pctERCC)
dim(scRNA)

FeatureScatter(scRNA,feature1 = "nCount_RNA", feature2 = "percent.mt")
scRNA<-FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 1500)
top50<-head(VariableFeatures(scRNA),50)
top50
print(VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 3))

print(FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt"))

print(FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))

scRNA<- subset(scRNA, subset = nFeature_RNA > 200 & percent.mt < 10)
print(dim(scRNA))

scRNA = NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA<-ScaleData(scRNA,features = (rownames(scRNA)))
GetAssayData(scRNA,slot = "scale.data",assay = "RNA")[1:8,1:4]
GetAssayData(scRNA,slot = "counts",assay = "RNA")[1:8,1:4]
scRNA = RunPCA(scRNA,features = VariableFeatures(scRNA))
scRNA = RunTSNE(scRNA,features = VariableFeatures(scRNA),check_duplicates = FALSE)
scRNA = RunUMAP(scRNA,features = VariableFeatures(scRNA))
DimPlot(scRNA,reduction='pca')
DimPlot(scRNA,reduction='umap',split.by = "orig.ident")
DimPlot(scRNA,reduction='tsne',split.by = "orig.ident")
scRNA=JackStraw(scRNA,reduction='pca',dims = 20)
scRNA=ScoreJackStraw(scRNA,dims = 1:20)
JackStrawPlot(scRNA,dims = 1:20,reduction = "pca")
ElbowPlot(scRNA,ndims = 20,reduction = "pca")

scRNA<-FindNeighbors(scRNA,dims = 1:20)

scRNA<-FindClusters(scRNA,resolution =c(0.5))
table(scRNA@meta.data$seurat_clusters)


DimPlot(scRNA, reduction = "tsne", label = T)
DimPlot(scRNA, reduction = "umap", label = T,pt.size = 1)
FeaturePlot(scRNA,reduction = "tsne",features = c("ZNRF2"))
FeaturePlot(scRNA,reduction = "umap",features = c("ZNRF2"))
col.num<-length(unique(scRNA@meta.data$seurat_clusters))
VlnPlot(scRNA,
        features=c("ALB"),
        group.by="seurat_clusters",
        cols=rainbow(col.num))
VlnPlot(scRNA,
        features=c("ZNRF2"),
        group.by="orig.ident",
        cols=rainbow(col.num))
scRNA<-NormalizeData(scRNA, normalization.method = "LogNormalize")

GetAssayData(scRNA,slot = "data",assay = "RNA")[1:8,1:4]
#if test.use id negbinom, possion or deseq2, slot will be set to counts
diff.willcox= FindAllMarkers(scRNA)
head(diff.willcox)
dim(diff.willcox)
library(tidyverse)
all.makers=diff.willcox %>% select(gene, everything()) %>%
  subset(p_val<0.05 & abs(diff.willcox$avg_log2FC)>0.5)
#an adjusted pvalue <0.05 and |log 2[fold change(FC)]|0.5
#were considered the 2 cutoff criteria for identifying maker genes
dim(all.makers)
summary(all.makers)

top10= all.makers%>% group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
top10
top10=CaseMatch(search=as.vector(top10$gene),match = rownames(scRNA))
top10
length(top10)
length(unique(sort(top10)))
DoHeatmap(scRNA,features = top10,group.by="seurat_clusters")
FeaturePlot(scRNA,features = c("ALB","FGG"),reduction = 'tsne')
FeaturePlot(scRNA,features = c("ACTA2","TAGL"),reduction = 'tsne')
FeaturePlot(scRNA,features = c("KDR","VWF"),reduction = 'tsne')
FeaturePlot(scRNA,features = c("HLA-DQB1","CD68"),reduction = 'tsne')
FeaturePlot(scRNA,features = c("IGJ","CD79A"),reduction = 'tsne')
FeaturePlot(scRNA,features = c("ZNRF2"),reduction = 'tsne')

scRNA$orig.ident<-scRNA$seurat_clusters
table(Idents(scRNA))
scRNA <- RenameIdents(scRNA, `0` ="HCC", `1` = "HCC", `2` = "HCC", `3` = "HCC",
                      `4` = "Myeloid-derived", `5` = "HCC", `6` = "Endothelial cells",
                      `7` = "HCC", `8` = "HCC", `9` = "CAF", `10` = "B cells",
                      `11` = "CAF")

DimPlot(scRNA,reduction = 'tsne')
scRNA <- StashIdent(scRNA, save.name = 'orig.ident')
VlnPlot(scRNA,
        features=c("ZNRF2"),
        group.by="orig.ident",
        cols=rainbow(col.num))
