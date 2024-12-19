#下载和载入包
install.packages("Seurat")
library(Seurat)
library(Matrix)
library(dplyr)

# 一、导入数据类型与要求
  # 设置当前路径
setwd("D:\\大五上\\R包基础\\R lesson\\2010305208_刘霄(Week5)")
  # 1.分别查看三个数据
mtx <- readMM("singlecell_data\\matrix.mtx")
genes=read.table("singlecell_data\\genes.tsv", header=F, sep="\t")
barcodes=read.table("singlecell_data\\barcodes.tsv", header=F, sep="\t")
  # 2.数据读取
data=Read10X(data.dir ="singlecell_data",gene.column = 1)
  # 3.创建Seurat对象
bcdata=CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)

# 二、Seurat数据的质控
  # 2.1 三个关键参数
head(bcdata@meta.data)
VlnPlot(bcdata, features = c("nFeature_RNA", "nCount_RNA"), group.by="orig.ident",ncol = 2)
FeatureScatter(bcdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  # 2.2 阈值设置
bcdata <- subset(bcdata, subset = nFeature_RNA > 2000 & nFeature_RNA < 2500)
VlnPlot(bcdata, features = c("nFeature_RNA", "nCount_RNA"), group.by="orig.ident",ncol = 2)
  
  # 2.3手动添加注释
table(bcdata@meta.data$orig.ident)
number=as.character(bcdata@meta.data$orig.ident)

type1=c("AU565","BT20","BT474","BT483","BT549","CAL51","CAL851","CAMA1","DU4475","EFM19","EVSAT","HCC1143")
type2=c("HCC1187","HCC1500","HCC1937","HCC1954","HCC38","HCC70")
type3=c("HDQP1","HS578T","JIMT1","KPL1", "MCF12A","MCF7","MDAMB361","MDAMB415","MDAMB436","MDAMB453","MDAMB468","MX1","T47D","ZR751")
my_index=ifelse(number %in% type1,"type1",
                ifelse(number %in% type2,"type2",
                       ifelse(number %in% type3,"type3",NA)))
# summary(factor(my_index))
bcdata=AddMetaData(bcdata,metadata = my_index,col.name = "type")

#三、Seurat数据的处理与分析
  # 3.1 数据标准化
bcdata <- NormalizeData(bcdata, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # 3.2鉴定高变基因
bcdata <- FindVariableFeatures(bcdata, selection.method = "vst", nfeatures = 2000)
top_10=head(VariableFeatures(bcdata), 10)
p1=VariableFeaturePlot(bcdata)
LabelPoints(plot = p1, points = top_10, repel = TRUE)
  
  # 3.3数据放缩
all.genes <- rownames(bcdata)
bcdata <- ScaleData(bcdata, features = all.genes)
  
  # 3.4PCA降维
bcdata <- RunPCA(bcdata, features = VariableFeatures(object = bcdata))
VizDimLoadings(bcdata, dims = 1:2, reduction = "pca")
DimPlot(bcdata, reduction = "pca",group.by = "orig.ident")
DimHeatmap(bcdata, dims = 1:6, cells = 500, balanced = TRUE,nfeatures = 8)
  
  # 3.5维度确定
ElbowPlot(bcdata)
  
  # 3.6细胞聚类
bcdata <- FindNeighbors(bcdata, dims = 1:10)
bcdata <- FindClusters(bcdata, resolution = 0.5)
  
  # 3.7非线性降维
bcdata <- RunUMAP(bcdata, dims = 1:10)
bcdata <- RunTSNE(bcdata, dims = 1:10)
# 基于orig.ident区分的很好，而基于type不太好
DimPlot(bcdata, reduction = "umap",group.by = "orig.ident")
DimPlot(bcdata, reduction = "tsne",group.by = "orig.ident")
# DimPlot(bcdata, reduction = "umap",group.by = "type")
# DimPlot(bcdata, reduction = "tsne",group.by = "type")

  # 3.8寻找差异表达基因
positive_markers <- FindAllMarkers(bcdata, only.pos = TRUE)

  # 3.9细胞簇命名
levels(bcdata)#查看原始簇名
new_cluster_names=paste0("cluster",c(1:27))#新的簇名
names(new_cluster_names)=levels(bcdata)#赋值新的簇名
attributes(bcdata)#查看bcdata的属性
bcdata=RenameIdents(bcdata,new_cluster_names)#重命名簇名
DimPlot(bcdata, reduction = "umap", label = TRUE)#用新的簇名绘制UMAP图
  
# 四、个性化数据分析
    # 目标基因：ENSG00000138413
  # 4.1特定基因表达情况
    # 小提琴图
VlnPlot(bcdata, features = "ENSG00000138413",pt.size = 1, group.by = "type")
    # 山脊图
RidgePlot(bcdata, features ="ENSG00000138413",group.by = "type")
    # 结合降维结果与基因表达水平
FeaturePlot(bcdata, features = c("ENSG00000138413"), cols = c("lightgrey", "red"), reduction = "umap")
  # 4.2 多基因表达情况：以top_10为例
    #feature plot
FeaturePlot(bcdata, features = top_10)
    #dot plot
DotPlot(bcdata, features = top_10,group.by="orig.ident")
    # heatmap
DoHeatmap(bcdata, features = top_10, size = 3,group.by="type") 
  # 4.3数据提取
      # 4.3.1提取counts表达矩阵
scRNA_counts <- LayerData(bcdata,assay = 'RNA',layer = 'counts')%>%as.data.frame()
      # 4.3.2提取伪bulk数据
average_exp=AggregateExpression(bcdata,assays = 'RNA',group.by = "orig.ident")%>%as.data.frame()

#完结撒花#
