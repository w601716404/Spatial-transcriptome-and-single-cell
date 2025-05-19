rm(list=ls())
setwd("D:\\NEW-pancancer\\Singlecell2")
#加载包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
# install_github("immunogenomics/harmony")
library(harmony)
# remove.packages("Seurat")

## 将文件夹中的文件名存到"dir_name"中
#dir_name <- list.files("数据呢")
#dir_name
####Read10X(data.dir绝对路径读取####
dir_name=list.files("D:\\NEW-pancancer\\Singlecell2\\GSE200972RAW\\")
scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste("D:\\NEW-pancancer\\Singlecell2\\GSE200972RAW\\", dir_name[i], sep = ""))
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = dir_name[i],min.cells = 3, min.features = 300)
}
# a=Read10X("数据呢/T2//")



# for (i in seq_along(dir_name)){
#   assign(paste0("scs_data", i), Read10X(data.dir = paste0("E:/兼职/02_生信代分析/单细胞/1114/GSE145154_RAW/", dir_name[i], sep = "")))
# }
# 
# ### create seurat objects from cellranger files
# for (i in seq_along(dir_name)){
#   assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = dir_name[i], min.cells = 3))
# }
# # for 循环，它会遍历 samples 列表中的每一个样本





####批量计算线粒体和红细胞比例####
for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]#获取scRNAlist中的第i个Seurat对象
  # 计算线粒体比例
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")#计算以"MT-"开头的基因比例
  # 计算红细胞比例
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") #定义红细胞基因列表
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))#在Seurat对象的RNA数据中查找红细胞基因的索引位置
  HB_genes <- rownames(sc@assays$RNA)[HB_m]  # 获取匹配到的红细胞基因的行名
  HB_genes <- HB_genes[!is.na(HB_genes)]  # 删除NA值（未匹配到的基因）
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features=HB_genes)   #计算红细胞基因比例，将结果存储在名为"HB_percent"的新列中
  # 将sc赋值给scRNAlist[[i]]
  scRNAlist[[i]] <- sc
  # 删除sc
  rm(sc)
}
# 注释：以上代码通过循环遍历scRNAlist中的每个Seurat对象，计算线粒体和红细胞比例，并将结果存储在相应的列中。最后将更新后的Seurat对象重新赋值给scRNAlist。


####批量绘制质控前小提琴图####
violin_before <- list()
for(i in 1:length(scRNAlist)){
  violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                                pt.size = 0.01, 
                                ncol = 4) 
}
# 合并图片
violin_before_merge <- CombinePlots(plots = violin_before,nrow=length(scRNAlist),legend='none')
# 将图片输出到画板上
violin_before_merge
# 保存图片
ggsave("violin_before_merge.pdf", plot = violin_before_merge, width = 15, height =7)

violin_before

####批量过滤细胞、MT、HB基因####
scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                mt_percent < 10 & 
                HB_percent < 3 & 
                nCount_RNA < quantile(nCount_RNA,0.97) & 
                nCount_RNA > 1000)})
# 一般默认线粒体含量至少要小于20%，红细胞的数目要至少小于5%；
# 在这里我们将过滤严格一点，调整为：
# nFeature_RNA：每个细胞检测表达的基因数目大于300，小于7000；
# nCount_RNA:每个细胞测序的UMI count含量大于1000，且剔除最大的前3%的细胞；
# mt_percent:每个细胞的线粒体基因表达量占总体基因的比例小于10%；
# HB_percent：每个细胞红细胞基因表达量占总体基因的比例小于3%。
View(scRNAlist[[1]]@meta.data)
# ps:没有固定的阈值标准，要根据自己的数据调整参数不断尝试，才能找到最佳结果。


####merge合并样本####
scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])
## 统计细胞数
table(scRNAlist[[]]$orig.ident)
# HC1  HC2   T1   T2 
# 5872 7143 1681 1904 
##过滤后可视化
# 绘图
violin_after <- VlnPlot(scRNAlist,
                        features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                        pt.size = 0.01,
                        ncol = 4)
# 将图片输出到画板上 
violin_after
# 保存图片
ggsave("vlnplot_after_qc.pdf", plot = violin_after, width = 15, height =7) 


####数据归一化、筛选高变基因与PCA降维####
# harmony整合是基于PCA降维结果进行的。
scRNAlist <- NormalizeData(scRNAlist) %>% # 数据归一化处理
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% #筛选高变基因
  ScaleData() %>% #数据标准化
  RunPCA(npcs = 30, verbose = T)#npcs：计算和存储的PC数（默认为 50）
a=DimPlot(scRNAlist,reduction = "pca",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
a

##看下高变基因有哪些可视化
# 提取前15个高变基因ID
top15 <- head(VariableFeatures(scRNAlist), 15) 
plot1 <- VariableFeaturePlot(scRNAlist) 
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
# 合并图片
feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
feat_15
# 保存图片
ggsave(file = "feat_15.pdf",plot = feat_15,he = 10,wi = 15 )


# Seurat是5.0版在要运行下这个代码才能进行差异分析(Seurat V4可忽略此代码)：
scRNAlist <- JoinLayers(scRNAlist)#连接layers层的count数据



####细胞周期评分####
# 提取g2m特征向量
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNAlist))
# 提取s期特征向量
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNAlist))
# 对细胞周期阶段进行评分
scRNAlist <- CellCycleScoring(object=scRNAlist,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAlist=CellCycleScoring(object = scRNAlist, 
                           s.features = s_genes, 
                           g2m.features = g2m_genes, 
                           set.ident = TRUE)#set.ident 是否给每个细胞标注一个细胞周期标签
scRNAlist <- CellCycleScoring(object=scRNAlist,  g2m.features=g2m_genes,  s.features=s_genes)
# p4=VlnPlot(scRNAlist, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
#            ncol = 2, pt.size = 0.1)
# p4
scRNAlist@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()



####RunHarmony去批次####
# 整合需要指定Seurat对象和metadata中需要整合的变量名。
scRNA_harmony <- RunHarmony(scRNAlist, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
b
# 合并图片
pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated
#去批次前后对比


# 后续都是基于Harmony矫正之后的数据，不是基因表达数据和直接的PCA降维数据。
# 设置reduction = 'harmony'，后续分析是基于Harmony矫正之后的数据。
####聚类、umap/tsne降维降维####
ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:19) %>% FindClusters(resolution = 1.2)#分辨率可以自己调
##umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:19)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:19)

# 绘图
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
# 合并图片
umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
# 将图片输出到画板
umap_tsne_integrated
# 保存图片
ggsave("umap_tsne_integrated.pdf",umap_tsne_integrated,wi=25,he=15)
#保存数据
save(scRNA_harmony,scRNAlist,file = "scdata2.Rdata")
load("scdata2.Rdata")

table(scRNA_harmony@meta.data$seurat_clusters)

####细胞注释####


# 差异分析
markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)


#2.对计算好的每cluster的marker基因进行筛选
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
#筛选出P<0.05的marker基因
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
# top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
View(top10)
write.csv(top10,"12cluster_top10.csv",row.names = T)
# write.csv(top15,"2.基因注释/16cluster_top15.csv",row.names = T)
#其中①avg_lg2FC代表0cluster和除0之外的cluster进行差异分析后得到的变化倍数的lg2FC值。
#②P.val_adj:即P值，小于一定的值在认为有统计学意义。
#③以cluster0为例，pct.1代表在0这个cluster中这个基因表达的细胞比例，pct.2代表在非0cluster中，这个基因表达的细胞比例
#找marker基因时也要看pct.1和pct.2，因为我们要找的marker基因是在特定cluster中特异性高表达的基因，不能再很多cluster中都高表达。
#3.找到cellMarker网站，选择人的物种对每个cluster进行亚细胞注释(也可以通过查找文献获得marker基因对应的celltype)
#这里，我们根据ellMarker得各个细胞类型得marker基因进行标注

# 可视化maker
DoHeatmap(scRNA_harmony, features = top10$gene, slot="data") + NoLegend()#slot默认使用scaledata里边只有2k个高变gene表达 这里要使用data数据，不然有些gene会找不到表达量
VlnPlot(scRNA_harmony,features = top10$gene[1:10])#选前10个makergene看看,这里选的都是cluster0的差异基因

# 或者自己选gene看看
VlnPlot(scRNA_harmony,features = c("PTPRC","MS4A1"))


# ggsave(filename="plot.pdf",width = 210,height = 297,units = "mm")

#3.找到每一个cluster的celltype后，把每一个cluster命名成对应的celltype
View(scRNA_harmony@meta.data)#seurat_clusters这一列存放了每个细胞对应的cluster
table(scRNA_harmony@meta.data$seurat_clusters)#可以查看每个cluster有多少个细胞

# http://bio-bigdata.hrbmu.edu.cn/CellMarker/
####第二部分找细胞marker####
library(ggplot2) 
#先看下cluster1的点图，这里可以任意改top10$gene范围
p <- DotPlot(scRNA_harmony, features = top10$gene[1:10],
             assay='RNA' ,group.by = 'seurat_clusters' ) + coord_flip()+ggtitle("")
#整理的一些          'CSPG4',#Pericyte
                   'ENPP3','FCER1A','FCGR1A','FCGR2A', 'FCGR2B',#Mast cell
                   'AQP3', 'MUC5AC', 'MUC5B', 'PIGR', 'SCGB1A1',#Airway secretory cell
                   'THBD', 'CD1C','CLEC4C','CD83','HLA.DQA2','HLA.DQA1',#Dendritic cell
                   'KRT19','ITGAE', 'VCAM1', 'IL6R', 'ANPEP', 'CD24', 'CDH1',#Epithelial cell
                   'MUC1','ABCA3', 'LPCAT1', 'NAPSA', 'SFTPB', 'SFTPC', 'SLC34A2',#lung Epithelial cell
                   'ACTA2', 'PDGFRA', 'PDGFRB','THY1',#Fibroblast
                   'CD19', 'CD79A', 'MS4A1', # B cells
                   'TAGLN2','CD5',#Activated B cell
                   'CD27', 'CD38','LY9', 'LAIR1', 'ICAM1','KIT',  # plasma 
                   'NKG7','GNLY',#NK
                   'CD8B', 'ZNF683','FCGR3A', 'FCGR3B', 'NCAM1', 'KLRB1',#NKT
                   'CXCR5','CCR7',#memory T cell
                   'CD6','IL7R', 'IL2RA', 'IKZF2',#Treg
                   'GP1BA','SELL','IFNG','CXCR3','IL17A','IL4','GATA3',#T helper cell
                   'CD33', 'ENTPD1',#MDSC
                   'S100A8','S100A9','S100A12',#Neutrophil
                   'CD68',  'CD163','MRC1','MSR1','CXCL10','CCL18', ## Macrophage (belong to monocyte)
                   'PECAM1','VWF','MCAM','CD34','ESM1', ## Endothelial
                   'ALDH1A1','KRT18', 'PROM1',## Cancer stem cell
                   'HHIP','SFTPC','SFTPA','SFTPC','LAMP3',
                   'MDK','SFTPB'
                   
)
genes_to_check = c('PTPRC','CD3D','CD3E','CD4','CD8A')#Tcell marker





####12、提取上皮细胞亚聚类####
# 加载第一次聚类的数据
#load("sce_dim20_reso0.2_zhushi.Rdata")
#别忘了默认分组为celltype
Idents(scRNA) = "celltype2"
Idents(scRNA)

#提取上皮细胞
MAC=subset(scRNA,idents=c("Mon/Macr"))
table(MAC@meta.data$celltype2)
#亚群聚类常规方法是提取counts来走作标准流程，这样做肯定没问题；
#其他方法也有很多，但具有争议这里不讲。

# seuratV5要这样提取counts
GetAssayData(MAC, assay="RNA", layer='counts')
# V5这样提取count是没有行名的
MAC@assays$RNA@layers$counts#这样没有行名和列名，V5的锅

# 创建SeuratObject对象
MAC_sce = CreateSeuratObject(counts = GetAssayData(MAC, assay="RNA", layer='counts'),  # 使用提取的细胞构建新的Seurat对象
                             meta.data = MAC@meta.data)  # 保留meta.data所有信息


#标注化、归一化、高变基因、pca
MAC_sce <- NormalizeData(MAC_sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)#对数据进行归一化、找高变基因、均一化、进行PCA降维


#重点！！！！！！！！！
#这里要注意是否选择去批次，可选可不选都能解释通（建议都跑一遍）
# 比如我这里是上皮细胞里边有癌细胞，本身存在很强的样本异质性，去了批次反而去掉了样本间差异.
# 免疫细胞也可以不去批次，也存在较强异质性
# 我这里就不运行
# Epi_sce <- RunHarmony(Epi_sce, group.by.var = "orig.ident")

#降维聚类
MAC_sce
ElbowPlot(MAC_sce, ndims=50, reduction="pca")
MAC_sce <- FindNeighbors(MAC_sce, reduction = "pca", dims = 1:20)#reduction="harmony"
MAC_sce = FindClusters(MAC_sce,resolution = 0.8)
table(MAC_sce@meta.data$seurat_clusters)
# 0    1    2    3    4    5    6    7    8 
# 1134  705  150  133  106   94   75   21   18
MAC_sce <- RunUMAP(MAC_sce, reduction = "pca", dims = 1:20)##reduction="harmony"


#EPI_umap图_未注释
MAC_sce$group
plot1 =DimPlot(MAC_sce, reduction = "umap",label = T,raster=FALSE) 
plot1
plot2 = DimPlot(MAC_sce, reduction = "umap", group.by='orig.ident',raster=FALSE) 
plot2
plot3 = DimPlot(MAC_sce, reduction = "umap",split.by = "group",label = T,raster=FALSE) 
plot3
plot4 = DimPlot(MAC_sce, reduction = "umap",group.by = "group",shuffle = T,raster=FALSE) 
#自己手动保存图片吧
dev.off()



#保存未注释的上皮细胞数据#
save(MAC_sce,file = "mac_sce.Rdata")
#load("Epi_sce_9clu.Rdata")




####13、提取上皮细胞亚聚类####
# 加载第一次聚类的数据
#load("sce_dim20_reso0.2_zhushi.Rdata")
#别忘了默认分组为celltype
Idents(scRNA) = "celltype"
Idents(scRNA)

#提取上皮细胞
Epi=subset(scRNA,idents=c("Epithelial_cells"))
table(Epi@meta.data$celltype)
#亚群聚类常规方法是提取counts来走作标准流程，这样做肯定没问题；
#其他方法也有很多，但具有争议这里不讲。

# seuratV5要这样提取counts
GetAssayData(Epi, assay="RNA", layer='counts')
# V5这样提取count是没有行名的
Epi@assays$RNA@layers$counts#这样没有行名和列名，V5的锅

# 创建SeuratObject对象
Epi_sce = CreateSeuratObject(counts = GetAssayData(Epi, assay="RNA", layer='counts'),  # 使用提取的细胞构建新的Seurat对象
                             meta.data = Epi@meta.data)  # 保留meta.data所有信息


#标注化、归一化、高变基因、pca
Epi_sce <- NormalizeData(Epi_sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)#对数据进行归一化、找高变基因、均一化、进行PCA降维


#重点！！！！！！！！！
#这里要注意是否选择去批次，可选可不选都能解释通（建议都跑一遍）
# 比如我这里是上皮细胞里边有癌细胞，本身存在很强的样本异质性，去了批次反而去掉了样本间差异.
# 免疫细胞也可以不去批次，也存在较强异质性
# 我这里就不运行
# Epi_sce <- RunHarmony(Epi_sce, group.by.var = "orig.ident")

#降维聚类
Epi_sce
ElbowPlot(Epi_sce, ndims=50, reduction="pca")
Epi_sce <- FindNeighbors(Epi_sce, reduction = "pca", dims = 1:20)#reduction="harmony"
Epi_sce = FindClusters(Epi_sce,resolution = 0.8)
table(Epi_sce@meta.data$seurat_clusters)
# 0    1    2    3    4    5    6    7    8 
# 1134  705  150  133  106   94   75   21   18
Epi_sce <- RunUMAP(Epi_sce, reduction = "pca", dims = 1:20)##reduction="harmony"


#EPI_umap图_未注释
Epi_sce$group
plot1 =DimPlot(Epi_sce, reduction = "umap",label = T,raster=FALSE) 
plot1
plot2 = DimPlot(Epi_sce, reduction = "umap", group.by='orig.ident',raster=FALSE) 
plot2
plot3 = DimPlot(Epi_sce, reduction = "umap",split.by = "group",label = T,raster=FALSE) 
plot3
plot4 = DimPlot(Epi_sce, reduction = "umap",group.by = "group",shuffle = T,raster=FALSE) 
#自己手动保存图片吧
dev.off()



#保存未注释的上皮细胞数据#
save(Epi_sce,file = "Epi_sce.Rdata")
#load("Epi_sce_9clu.Rdata")










