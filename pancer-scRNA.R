
rm(list=ls())
setwd("D:\\NEW-pancancer\\Singlecell4")
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
# install_github("immunogenomics/harmony")
library(harmony)
library(data.table)
# remove.packages("Seurat")
# 查看当前工作目录
getwd()
#######################GSE207422_luad###########################
tpm<- fread("GSE207422_luad\\raw\\all_count_clean_new.txt",sep="\t",header=T,data.table = F)
head(tpm)[1:6]

rownames(tpm) <- tpm[,1]
tpm[,1] <- NULL
head(tpm)[1:6]
#test.seu1=CreateSeuratObject(counts=log(tpm+1))
GSE207422_LUAD=CreateSeuratObject(counts=tpm,
                                  min.features = 200,
                                  min.cells = 3, 
                                  project = "GSE207422_LUAD")
## 3.在GitHub储存库中r包
library(devtools)
###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(limma)
library(patchwork)
metadata<-GSE207422_LUAD@meta.data
# 使用sub函数替换掉最后一个下划线及其后面的所有内容
metadata$orig.ident<-row.names(metadata)
metadata$orig.ident <- sub("_[^_]+$", "", metadata$orig.ident)
#patients_metadata <- data.table::fread("GSE207422_luad\\raw\\xuanze.txt", header = TRUE)
#metadata$cellid<-rownames(metadata)
#metadata <- left_join(x = metadata, y = patients_metadata, by = "sample")
#row.names(metadata)<-metadata$cellid
#metadata2<-merge(metadata,patients_metadata,by = "sample")
### 通过AddMetaData直接补入4列
GSE207422_LUAD<- AddMetaData(GSE207422_LUAD, metadata = metadata)
table(GSE207422_LUAD@meta.data$orig.ident)
GSE207422_LUADsubseu<-subset(GSE207422_LUAD, subset = orig.ident %in% c('BD_immune05','BD_immune03','BD_immune13','BD_immune14',
                                                                    'BD_immune10','BD_immune04','BD_immune01','BD_immune12',
                                                                    'BD_immune06'))
rm(GSE207422_LUAD)
############################GSE235672_GBM################################
#########################h
dir='GSE235672_GBM\\h5\\'
samples=list.files(dir)
samples
scRNAlist1=lapply(samples,function(pro){
  print(pro)
  sce=CreateSeuratObject(counts=Read10X_h5(file.path(dir,pro)),
                         #project=gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',pro)),
                         project=gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_ST_','',pro)),
                         min.cells=5,
                         min.features=300)
  return(sce)
})

scRNAlist1
samples
#########10X
dir_name <- list.files('GSE235672_GBM\\raw/') 
dir_name
# [1] "P32N" "P32T" "P33N" "P33T"
# ??ʼ??һ?????б?��?洢Seurat????
scRNAlist2 <- list()
# ??ʼ????Ŀ¼???б?
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = 
                      paste("D:\\NEW-pancancer\\Singlecell4\\GSE235672_GBM\\raw\\", dir_name[i], sep = ""))
  scRNAlist2[[i]] <- CreateSeuratObject(counts, project = dir_name[i],min.cells = 3, min.features = 300)
}

scRNAlist <- append(scRNAlist1, scRNAlist2)
rm(scRNAlist1, scRNAlist2)

scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])

scRNAlist<-subset(scRNAlist, subset = orig.ident %in% c('Pt15','Pt2','Pt16'))
metadata<-scRNAlist@meta.data
metadata$sample<-metadata$orig.ident
metadata$orig.ident<-"GSE235672"
metadata$cellid<-rownames(metadata)
patients_metadata <- data.table::fread("GSE235672_GBM\\xuanze.txt", header = TRUE)
metadata <- left_join(x = metadata, y = patients_metadata, by = "sample")
row.names(metadata)<-metadata$cellid
### 通过AddMetaData直接补入4列
scRNAlist<- AddMetaData(scRNAlist, metadata = metadata)
table(scRNAlist@meta.data$sample)
GSE235672_GBMsubseu<-scRNAlist
rm(scRNAlist)
#######################GSE229353_LC###############################
dir_name=list.files("D:\\NEW-pancancer\\Singlecell4\\GSE229353_LC\\raw\\")
scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste("D:\\NEW-pancancer\\Singlecell4\\GSE229353_LC\\raw\\", dir_name[i], sep = ""))
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = dir_name[i],min.cells = 3, min.features = 300)
}
scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])
scRNAlist<-subset(scRNAlist, subset = orig.ident %in% c('P05','P07','P04'))
metadata<-scRNAlist@meta.data
metadata$sample<-metadata$orig.ident
metadata$orig.ident<-"GSE229353"
metadata$cellid<-rownames(metadata)
patients_metadata <- data.table::fread("GSE229353_LC\\xuanze.txt", header = TRUE)
metadata <- left_join(x = metadata, y = patients_metadata, by = "sample")
row.names(metadata)<-metadata$cellid
scRNAlist<- AddMetaData(scRNAlist, metadata = metadata)
table(scRNAlist@meta.data$sample)
GSE229353_LCsubseu<-scRNAlist
rm(scRNAlist)
##############GSE212217_EC###################
dir='D:\\NEW-pancancer\\Singlecell4\\GSE212217_EC\\GSE212217_RAW\\'
samples=list.files(dir)
samples
scRNAlist1=lapply(samples,function(pro){
  print(pro)
  sce=CreateSeuratObject(counts=Read10X_h5(file.path(dir,pro)),
                         #project=gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',pro)),
                         project=gsub('_scRNA_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',pro)),
                         min.cells=5,
                         min.features=300)
  return(sce)
})
scRNAlist1 <- merge(x=scRNAlist1[[1]],y=scRNAlist1[-1])
metadata<-scRNAlist1@meta.data
metadata$sample<-metadata$orig.ident
metadata$orig.ident<-"GSE212217"
metadata$cellid<-rownames(metadata)
patients_metadata <- data.table::fread("GSE212217_EC\\xuanze.txt", header = TRUE)
metadata <- left_join(x = metadata, y = patients_metadata, by = "sample")
row.names(metadata)<-metadata$cellid
scRNAlist1<- AddMetaData(scRNAlist1, metadata = metadata)
table(scRNAlist1@meta.data$sample)
GSE212217_ECsubseu<-scRNAlist1
rm(scRNAlist1)

#####################GSE169246_TNBC##################################
count <-Read10X(data.dir = 'GSE169246_TNBC\\data\\', gene.column=1)
seu<-CreateSeuratObject(counts = count, project = "TNBC_GSE169246",
                        min.features = 300,
                        min.cells = 3)
meta<-seu@meta.data
group=sapply(strsplit(rownames(meta),"\\."),"[",2)
seu$orig.ident=group

#sample<-as.data.frame(table(seu$sampleid))
#colnames(sample)<-c('sample_id','number')
#write.table(sample,'sample.csv',sep = ',',row.names = F)
GSE169246_TNBCsubseu<-subset(seu, subset = orig.ident %in% c('Pre_P001_b','Pre_P005_b',
                                                         'Pre_P007_b','Pre_P012_b',
                                                         'Pre_P017_b','Pre_P019_b'))
rm(seu,count)

table(GSE169246_TNBCsubseu@meta.data$orig.ident)#查看默认分组
#patients_metadata <- data.table::fread("GSE169246_TNBC/sample.txt", header = TRUE)
#metadata <- FetchData(GSE169246_TNBCsubseu, "sample")
#metadata$cellid <- rownames(metadata)
### dplyr中的left_join合并
#metadata <- merge(x = metadata, y = patients_metadata, by.x="sample",by.y ="sample")
### 重新加行名
#rownames(metadata) <- metadata$cellid
### 通过AddMetaData直接补入4列
#GSE169246_TNBCsubseu <- AddMetaData(GSE169246_TNBCsubseu, metadata = metadata)
#table(GSE169246_TNBCsubseu@meta.data$response)#看一下分组

###############GSE145281_BLCA######################

tpm=read.table(file = "GSE145281_BLCA\\All_matrix.txt", header=T, sep="\t", check.names=F)
scRNA=CreateSeuratObject(counts=tpm)
metadata<-scRNA@meta.data
#metadata$sample<-metadata$orig.ident
#metadata$orig.ident<-"GSE145281"
#metadata$cellid<-rownames(metadata)
# 使用sub函数替换掉最后一个下划线及其后面的所有内容
##metadata$sample <- sub("_[^_]+$", "", metadata$sample)
#patients_metadata <- data.table::fread("GSE145281_BLCA\\xuanze.txt", header = TRUE)
#metadata <- left_join(x = metadata, y = patients_metadata, by = "sample")
#row.names(metadata)<-metadata$cellid
### 通过AddMetaData直接补入4列
#scRNA<- AddMetaData(scRNA, metadata = metadata)
#table(scRNA@meta.data$sample)
GSE145281_BLCAsubseu<-scRNA
rm(scRNA)
###############GSE120575_SKCM######################
headers <- fread("GSE120575_SKCM\\GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",nrows=2,header=F,data.table = F)
headers[,1] <- NULL
tpm <- fread("GSE120575_SKCM\\GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",skip=2,header=F,data.table = F)
rownames(tpm) <- tpm[,1]
tpm[,1] <- NULL
tpm[,16292] <- NULL
colnames(tpm) <- headers[1,]
meta.sample <- headers[2,]

seurat_data3<-tpm
SKCM_GSE120575 <- CreateSeuratObject(counts = log(seurat_data3+1),
                                     min.features = 200,
                                     min.cells = 3, 
                                     project = "SKCM_GSE120575")
dt<-SKCM_GSE120575@meta.data
dt$cell<-rownames(dt)
header<-t(headers)
colnames(header)<-c("cell","orig.ident")
header<-as.data.frame(header)
dt2<-merge(dt,header,by.x="cell",by.y="cell")
rownames(dt2) <- dt2$cell
dt2<-dt2[,-2]
colnames(dt2)[4]<-c("orig.ident")
SKCM_GSE120575<- AddMetaData(SKCM_GSE120575 , metadata = dt2)

#table(SKCM_GSE120575@meta.data$orig.ident)#查看默认分组
#patients_metadata <- data.table::fread("GSE120575_SKCM/sample2.txt", header = TRUE)
#dt3<-merge(x = dt2, y = patients_metadata, by.x="sample",by.y ="sample")
####分组方法1—join函数
### 重新加行名
#rownames(dt3) <- dt3$cell
#rownames(metadata) <- metadata$cell_id
### 通过AddMetaData直接补入4列
#SKCM_GSE120575 <- AddMetaData(SKCM_GSE120575, metadata = dt3)
table(SKCM_GSE120575@meta.data$orig.ident)#看一下分组

GSE120575_SKCMsubseu<-subset(SKCM_GSE120575, subset = orig.ident %in% c('Pre_P29','Post_P8_T_enriched','Pre_P1','Pre_P3',
                                                                        'Pre_P25','Post_P5'))

########合并###########################################
#scRNAlist=list(GSE207422_LUADsubseu,scRNAlist)
scRNAlist=list(GSE207422_LUADsubseu,GSE169246_TNBCsubseu,GSE145281_BLCAsubseu,GSE120575_SKCMsubseu)

scRNAlist <- append(scRNAlist, scRNAlist3)

#scRNAlist=list(GSM5509264_diseased,GSM5509265_diseased,GSM4365601_healthy,GSM4365609_healthy)
save(scRNAlist,file = "scRNAlist.Rdata")

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
#查看第一个样本的小提琴图
violin_before[[1]] 

####批量过滤细胞、MT、HB基因####
scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                mt_percent < 15 & 
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
#rm(seu,seurat_data1,seurat_data2,seurat_data3,SKCM_GSE115978,SKCM_GSE115978subseu,SKCM_GSE120575,SKCM_GSE120575subseu
#   ,BCC_GSE123813,BCC_GSE123813subseu,header,headers,TNBC_GSE169246subseu,tpm,dt,dt2,count,meta,meta.sample,metadata,patients_metadata,dt3)
rm(GSE207422_LUADsubseu,GSE235672_GBMsubseu,GSE229353_LCsubseu,GSE212217_ECsubseu,
   GSE169246_TNBCsubseu,GSE145281_BLCAsubseu,GSE120575_SKCMsubseu,count,dt,dt2,dt3,header,headers,meta,meta.sample,metadata,patients_metadata,
   seurat_data3,SKCM_GSE120575,tpm,counts, scRNAlist3)
####merge合并样本####
scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])
#rm## 统计细胞数
table(scRNAlist[[]]$orig.ident)
# GSM4365601_healthy  GSM4365609_healthy GSM5509264_diseased GSM5509265_diseased 
# 438                8191                7777                9070 
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
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:26) %>% FindClusters(resolution = 1.5)#分辨率可以自己调
##umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:26)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:26)

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
save(scRNA_harmony,file = "scdata2.Rdata")
load("scdata2.Rdata")

table(scRNA_harmony@meta.data$seurat_clusters)

####细胞注释####
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
####细胞注释####
library(ggplot2) 
#先看下cluster1的点图，这里可以任意改top10$gene范围
p <- DotPlot(scRNA_harmony, features = top10$gene[1:10],
             assay='RNA' ,group.by = 'seurat_clusters' ) + coord_flip()+ggtitle("")
#整理的一些marker
genes_to_check = c('PTPRC','CD3D','CD3E','CD4','CD8A',## Tcells
                   'CSPG4',#Pericyte
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
##Ñ¡????????Í¼Õ¹Ê¾
# T Cells
select_genes <- c("CD3D", "CD3E","CD8A","CD3A","IL7R")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#B cells 
select_genes <- c("CD19", "CD79A","MS4A1","CD20")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#Plasma cells
select_genes <- c("IGHG1","MZB1","SDC1","CD79A")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#Monocytes and macrophages
select_genes <- c("CD68", "CD163", "CD14")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#NK Cells 
select_genes <- c("FGFBP2","FCG3RA","CX3CR1","NKG7")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#Photoreceptor cells
select_genes <- c("RCVRN")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
# Fibroblasts
select_genes <- c("FGF7","MME","ACTA2","COL1A1")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
# Endothelial cells
select_genes <- c("PECAM1", "VWF")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
# epi or tumor
select_genes <- c("EPCAM","KRT19","PROM1", "ALDH1A1", "CD24")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
# immune
select_genes <- c("CD45","PTPRC")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#epithelial/cancer
select_genes <- c("EpCAM")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1

#myoepithelial cells
select_genes <- c("ACTA2","TAGLN")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#basal epithelial
select_genes <- c("KRT14","ITGA6","KRT5","TP63","KRT17","MME")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1

#luminal epithelial  KRT8 KRT18 KRT19 FOXA1 GATA3  MUC1 CD24 
select_genes <- c("KRT8","KRT18","KRT19","FOXA1","GATA3","MUC1","CD24")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#luminal progenitor KIT GABRP 
select_genes <- c("KIT","GABRP")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#stroma  FAP  COL1A1 COL3A1 COL5A1 ACTA2 TAGLN LUM FBLN1 COL6A3 COL1A2 COL6A1 COL6A2 
select_genes <- c("FAP","COL1A1","COL3A1","ACTA2","TAGLN","LUM","FBLN1","COL6A3","COL1A2","COL6A1","COL6A2")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#endothelial
select_genes <- c("PECAM1","VWF","CDH5","SELE")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#T cell
select_genes <- c("CD2","CD3D","CD3E","CD3G","CD8A","CD8B")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#basal 
select_genes <- c("KRT5","ACTA2","MYLK","SNAI2")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#luminal progenitor
select_genes <- c("TNFRSF11A","KIT")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#mature luminal cells
select_genes <- c("ESR1","PGR","FOXA1")
p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1

#?Ú¶??Ö·?????SingleR????Ï¸??????
library(SingleR)
##??Ê¦???????Ä°Ù¶??Æ´??? ???????Ðµ??Ëµ????Ý¿â£¬??Îª????Ã»??vpn??????singler?????Ý¿?Ã»??????
###???Øº????Ý¿??ó£¬°?ref_Human_all.Rdata???Øµ??????Ð£????????Ç¶????Ý¿??Ä¼??Ø£??Í¿??Ô°???singler???ã·¨À´??Ï¸????Èº???Ð¶????Ë¡?
load("ref_Human_all.RData")
###???Ç¿??Ô¿????Ú»????Ð¶???Ò»????ref_Human_all???Ä¼? ??Ð¡Îª113mb  ???????????Ý¿?
####È»?????Ç°Ñ»????Ðµ?ref_Human_all??Öµ??refdata
refdata <- ref
###??rna??×ªÂ¼??????????È¡
?GetAssayData
testdata <- GetAssayData(scRNA_harmony, slot="data")
###??scRNA?????Ðµ?seurat_clusters??È¡??À´??×¢???????????????Íµ?
clusters <- scRNA_harmony@meta.data$seurat_clusters
DefaultAssay(scRNA_harmony)
###??Ê¼??singler????
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
###????Ï¸?????Íµ?×¢???Ä¼?
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
###????Ò»??
write.csv(celltype,"celltype_singleR.csv",row.names = FALSE)
celltype[22,2]<-c("T_cells")
celltype[30,2]<-c("Monocyte")
celltype[32,2]<-c("Fibroblasts")
celltype[33,2]<-c("T_cells")
celltype[15,2]<-c("Epithelial_cells")

##??singler??×¢??Ð´??metadata?? ??Á½?Ö·???
###????Ò»
scRNA_harmony@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
###??Îª?Ò°?singler??×¢?Í¼??Øµ?metadata??Ê±?????????????Ö½?celltype?????Ô»?Í¼Ê±????group.by="celltype"
#seurat_clusters
DimPlot(scRNA_harmony, group.by="seurat_clusters", label=T, label.size=5, reduction='tsne',raster=FALSE)
DimPlot(scRNA_harmony, group.by="celltype", label=T, label.size=5, reduction='tsne',raster=FALSE)
DimPlot(scRNA_harmony, group.by="celltype", label=T, label.size=5, reduction='umap',raster=FALSE)
DimPlot(scRNA_harmony, group.by="response", label=T, label.size=5, reduction='tsne',raster=FALSE)
DimPlot(scRNA_harmony, group.by="response", label=T, label.size=5, reduction='umap',raster=FALSE)
scRNA<-scRNA_harmony
rm(scRNA_harmony)
save(scRNA,file="scRNA.rdata")
############增加分组
rm(list=ls())
load("scRNA.rdata")
metadata<-scRNA@meta.data
colnames(meta)
patients_metadata <- data.table::fread("clinical_zong.txt", header = TRUE)
metadata$cellid<-rownames(metadata)
metadata <- left_join(x = metadata, y = patients_metadata, by = "orig.ident")
row.names(metadata)<-metadata$cellid
#metadata<-metadata[,-6]
#metadata2<-merge(metadata,patients_metadata,by = "sample")
### 通过AddMetaData直接补入4列
scRNA<- AddMetaData(scRNA, metadata = metadata)
table(scRNA@meta.data$response)
save(scRNA,file="scRNA.rdata")
#"GSM4909299","GSM4909305","GSM4909313","GSM4909308","GSM4909306","GSM4909290","GSM4909287","GSM4909281",
#"GSM4909285","GSM4909282"
scRNA<-scRNA_harmony
rm(scRNA_harmony)
###????????
#celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F) 
#scRNA1@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
###??Îª?Ò°?singler??×¢?Í¼??Øµ?metadata??Ê±?????????????Ö½?singleR?????Ô»?Í¼Ê±????group.by="singleR"
#DimPlot(scRNA1, group.by="singleR", label=T, label.size=5, reduction='tsne')
###???Ç¿??Ô¿???  Á½?Ö·????Ãµ??Ä½???????Ò»???Ä£??????Ò±È½?Ï²???Ú¶??Ö·???
save(scRNA,file="scRNA.rdata")

DimPlot(scRNA, group.by="celltype", label=T, label.size=5, reduction='umap',raster=FALSE)
DimPlot(scRNA, group.by="response", label=T, label.size=5, reduction='tsne',raster=FALSE)
DimPlot(scRNA, group.by="GEO", label=T, label.size=5, reduction='tsne',raster=FALSE)
select_genes <- c("MDH1")
p1 <- VlnPlot(scRNA, features = select_genes, pt.size=0, group.by="response", ncol=2)
p1

DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters")
?VlnPlot
VlnPlot(scRNA, features = c("MDH1"), group.by="response")
FeaturePlot(scRNA, features = c("MDH1"), split.by="response")
?FeaturePlot

#########################################提取巨噬细胞##############
rm(list=ls())
load("scRNA.rdata")
####??È¡Ï¸???Ó¼?--??È¡scRNA1?Ðµ?????Ï¸??
##?Ò¹??Úº?Ò²???? ??È¡Ï¸???Ä·?????Á½??
####????Ò»????which????
phe=scRNA@meta.data
###?é¿´Ò»??immune_annotation?Ðµ???Ï¢????Òª??Îª?Ë·??ã¸´??Õ³??
table(phe$celltype)
###??whichÆ¥??immune??È»???Ãµ?Ï¸????ID ????metadata??????
#cells.use <- row.names(scRNA@meta.data)[which(phe$celltype=='Macrophage')]
scRNA@meta.data[which(scRNA@meta.data$celltype =="Macrophage"),"celltype2"]<-"Mon/Macr"
scRNA@meta.data[which(scRNA@meta.data$celltype =="Monocyte"),"celltype2"]<-"Mon/Macr"
###??scRNA1?Ðµ?????Ï¸????È¡??À´
#scRNAsub <-subset(scRNA, cells=cells.use)
##?é¿´Ò»???Ð¶?????????Ï¸??
#scRNAsub
##??????????Á½??subset????
cells.use1 <- subset(scRNA@meta.data, celltype2=="Mon/Macr")
###??scRNA1?Ðµ?????Ï¸????È¡??À´
scRNAsub <- subset(scRNA, cells=row.names(cells.use1))
##?é¿´Ò»???Ð¶?????????Ï¸??
#scRNAsub
###???Ç·???scRNAsub??scRNA.imm??Ò»????Ï¸????Ä¿ Ëµ????È¡Ã»?Ð´?

##???Â½?Î¬???????Â½?Î¬???à£º??Îª?Ù¾?????Ï¸??Ö®???????È½?Ð¡?????Ô¾??àº¯??FindClusters()???Æ·Ö±??ÊµÄ²??????????ßµ?resolution = 0.9
##PCA??Î¬????Ö®Ç°??Ò»?? ???ï²»?Ø¸?×¸????
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub, ndims=20, reduction="pca")
pc.num=1:17
##Ï¸?????à£¬??Ö®Ç°??Ò»?? ???ï²»?Ø¸?×¸????
scRNAsub <- FindNeighbors(scRNAsub, dims = 1:17) 
scRNAsub <- FindClusters(scRNAsub, resolution = 0.8)
table(scRNAsub@meta.data$seurat_clusters)
metadata <- scRNAsub@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
#write.csv(cell_cluster,'cell_cluster.csv',row.names = F)
##?????Ô½?Î¬ ????Ö®Ç°??Ò»?? ???ï²»?Ø¸?×¸????
#tSNE
scRNAsub = RunTSNE(scRNAsub, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub, 'tsne')
#write.csv(embed_tsne,'embed_tsne.csv')
plot1 = DimPlot(scRNAsub, reduction = "tsne") 
plot1
ggsave("tSNE.pdf", plot = plot1, width = 8, height = 7)
#UMAP
scRNAsub <- RunUMAP(scRNAsub, dims = pc.num)
embed_umap <- Embeddings(scRNAsub, 'umap')
#write.csv(embed_umap,'embed_umap.csv') 
plot2 = DimPlot(scRNAsub, reduction = "umap") 
plot2
ggsave("UMAP.pdf", plot = plot2, width = 8, height = 7)

#?Ï²?tSNE??UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
plotc
ggsave("tSNE_UMAPmac.pdf", plot = plotc, width = 10, height = 5)

#?ßµÍ·???????Î¢????????
DimPlot(scRNAsub, reduction = "umap",group.by = "celltype")#?ßµÍ·???????????Î¢????????
DimPlot(scRNA, reduction = "umap",split.by = "risk",label = T) + NoLegend()
?DimPlot





####????À´????×¢??  ×¢?ÍµÄ·??????????? Ç°?????????Ð½??? ?????????Í²???×¸????
##Cluster????????---??Êµ??????marker??????????marker????È¥??Õ¾  ???Ö¶?????????×¢??
#diff.wilcox = FindAllMarkers(scRNAsub)
#all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
#top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#write.csv(all.markers, "subcluster/diff_genes_wilcox.csv", row.names = F)
#write.csv(top10, "subcluster/top10_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub,file="scRNAmac.rdata")
DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='tsne',raster=FALSE)
DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='umap',raster=FALSE)

#
select_genes <- c("MDH1")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#
select_genes <- c("CTSD")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
# 
select_genes <- c("FOLR2")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#
select_genes <- c("HSPA6")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#
select_genes <- c("IFIT1")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
select_genes <- c("MKI67")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
select_genes <- c("S100A9")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
select_genes <- c("SPP1")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
select_genes <- c("TFPI")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
##Ï¸?????Í¼???
library(SingleR)
library(celldex)
##SingleRÏ¸???????????????????Çµ??Ãµ??????ßµ????Ý¿â£¬???Ï´Î²?Ò»??Ó´  ?Ï´????Ãµ??Ëµ????Ý¿?À´×¢?Íµ?
###????Í¬???Ä·??? ??Ê¦???????Ä°Ù¶????Ä¼????Øµ??????Ô¼??Ñ¾????Ãº?Â·?????Ä¼????Ð£???Òª?Â´??Ø·???  ????ref_Mouse_all.Rdata???Ä¼??????Â´??Ë£?????Ö®??Òª???Øµ???????
###??Ö®Ç°Ò²????  ?????Ä¼?????????Á½?Ö·?Ê½ Ò»?????Ö¶?  Ò»?????Ãº?????È¡  ?????????Ãº?????È¡
#load("ref_Monaco_114s.RData")
ref_Monaco <- MonacoImmuneData()
####???Ç½????Øµ????Ý¸?Öµ??refdata
refdata <- ref_Monaco
###????À´????????Ö®Ç°??singler×¢????Ê½??Ò»Ä£Ò»???? ???ï²»??×¸??
testdata <- GetAssayData(scRNAsub, slot="data")
clusters <- scRNAsub@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"celltype_singleR.csv",row.names = F)
###?????Ç½?singler×¢?Íµ???Ï¢???Øµ?metadata??

celltype[1,2]<-c("low_MDH1")
celltype[2,2]<-c("low_MDH1")
celltype[3,2]<-c("low_MDH1")
celltype[4,2]<-c("high_MDH1")
celltype[5,2]<-c("low_MDH1")
celltype[6,2]<-c("high_MDH1")
celltype[7,2]<-c("high_MDH1")
celltype[8,2]<-c("low_MDH1")
celltype[9,2]<-c("low_MDH1")
celltype[10,2]<-c("high_MDH1")
celltype[11,2]<-c("low_MDH1")
celltype[12,2]<-c("low_MDH1")
celltype[13,2]<-c("high_MDH1")
celltype[14,2]<-c("low_MDH1")
celltype[15,2]<-c("low_MDH1")
celltype[16,2]<-c("high_MDH1")
celltype[17,2]<-c("high_MDH1")
celltype[18,2]<-c("high_MDH1")
celltype[19,2]<-c("high_MDH1")
celltype[20,2]<-c("high_MDH1")
celltype[21,2]<-c("low_MDH1")
celltype[22,2]<-c("low_MDH1")
celltype[23,2]<-c("high_MDH1")
celltype[24,2]<-c("low_MDH1")
celltype[25,2]<-c("high_MDH1")
celltype[26,2]<-c("low_MDH1")
celltype[27,2]<-c("high_MDH1")

scRNAsub@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNAsub@meta.data[which(scRNAsub@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='tsne')
p1
p2 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='umap')
p2
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
p3
p4 = DimPlot(scRNAsub, group.by="seurat_clusters", label=T, label.size=5, reduction='tsne')
p4
p5 = DimPlot(scRNAsub, group.by="seurat_clusters", label=T, label.size=5, reduction='umap')
p5
ggsave("tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("UMAP_celltype.pdf", p2, width=7 ,height=6)
ggsave("celltype.png", p3, width=10 ,height=5)

save(scRNAsub,file="scRNAmac.rdata")
####??????marker??????singler×¢?Í£????Ç»???Ò»?Ö·??????Ç¸??????Ðµ?????Ñ§ÖªÊ¶???????×£???dotplotÀ´?Ö¶?×¢??
####????????Ò²??Ö®Ç°??Ò»?? ????Ò²????×¸???Ë£??????Ðµ??? ?????Ç±È½Ï¾?×¼??Ò»??×¢?Í·???

#####multi-group
################################################
### ?????Ä±?
## ???????????????ÂµÄ»????????????ã£©
data <- as.data.frame(table(scRNA$risk,scRNA$celltype))
colnames(data) <- c("risk","celltype","Freq")
library(dplyr)
df <- data %>% 
  group_by(risk) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = Freq/Total)

df$celltype  <- factor(df$celltype,levels = unique(df$celltype))

#write.csv(df,file = "output/cell_percent.csv",row.names = F,quote = F)

library(ggplot2)
p <- ggplot(df, aes(x = risk, y = Percent, fill = celltype)) +
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 1, width = 0.95) +
  #scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
p

### ??????Ê½
### facet_wrap ???????????Ð²???
ggplot(df,aes(x = risk, y = Percent,fill=risk))+
  geom_bar(stat="identity")+
  facet_wrap(~ celltype, nrow = 2)+
  theme_bw()+
  NoLegend()

### É¸Ñ¡???? ??Ç°????seurat??subset?È½Ï£????â·º?Íº???
data <- subset(df, subset = !CellType %in% c("Mono/Mk Doublets", "pDC", "Mk", "Eryth"))
ggplot(data,aes(x = group, y = Percent,fill=group))+
  geom_bar(stat="identity")+
  facet_wrap(~ CellType, nrow = 2)+
  theme_bw()+
  NoLegend()

ggplot(data,aes(x = CellType, y = Percent,fill=group))+
  geom_bar(stat="identity")+
  facet_wrap(~ group, nrow = 2)+
  theme_bw()+
  NoLegend()


################################################
### ???Ó»?Õ¹Ê¾DotPlot ?? Clustered_DotPlot
all_markers <- readRDS(file = "output/Seurat_stim_all_markers.rds")
library(dplyr)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  slice(1:5) %>%#Ñ¡Ç°5??
  ungroup() %>%
  pull(gene) %>%
  unique()

DotPlot(scobj, features = top_markers)
DotPlot(scobj, features = top_markers) + RotatedAxis()
DotPlot(scobj, features = top_markers) + coord_flip()+ RotatedAxis()
DotPlot(scobj, features = top_markers,dot.scale = 3) + coord_flip()+ RotatedAxis()+
  theme(axis.text.y = element_text(size = 8))

scCustomize::Clustered_DotPlot(scobj, features = top_markers)

### ????dotplot
### ×¢??????, ??Òª??????Ê±????+?Å£??áµ¼?Â»?Í¼Ê§?Ü£?????CD8+ T
Idents(scRNA) <- 'celltype'
levels(Idents(scRNA))


dput(levels(Idents(scRNA)))
c("Epithelial_cells", "T_cells", "Endothelial_cells", "Tissue_stem_cells", "Macrophage", 
  "CMP", "B_cell")
### ????Ë³??
Idents(scRNA) <- factor(Idents(scRNA), 
                        levels = c("Epithelial_cells", "T_cells", "Endothelial_cells", "Fibroblasts", "Macrophage", 
                                   "Neutrophils", "B_cell","NK_cell"))
### Ñ¡????Í¼????
markers.to.plot <- c("MDH1")
### ??Í¼,??ÒªÊ±??cols??????????Á½????É«??split.by???Ö¶???
DotPlot(scRNA, features = markers.to.plot, 
        cols = c("blue", "red"), 
        dot.scale = 6, 
        split.by = "risk") +
  RotatedAxis()

DotPlot(scRNA, features = markers.to.plot, 
        cols = c("blue", "red"), 
        dot.scale = 6, 
) +
  RotatedAxis()

rm(list=ls())
###############################提取巨噬细胞#########################
load("scRNA.rdata")
####??È¡Ï¸???Ó¼?--??È¡scRNA1?Ðµ?????Ï¸??
##?Ò¹??Úº?Ò²???? ??È¡Ï¸???Ä·?????Á½??
####????Ò»????which????
phe=scRNA@meta.data
###?é¿´Ò»??immune_annotation?Ðµ???Ï¢????Òª??Îª?Ë·??ã¸´??Õ³??
table(phe$celltype)
###??whichÆ¥??immune??È»???Ãµ?Ï¸????ID ????metadata??????
#cells.use <- row.names(scRNA@meta.data)[which(phe$celltype=='Macrophage')]
#scRNA@meta.data[which(scRNA@meta.data$celltype =="Macrophage"),"celltype2"]<-"Mon/Macr"
#scRNA@meta.data[which(scRNA@meta.data$celltype =="Monocyte"),"celltype2"]<-"Mon/Macr"
###??scRNA1?Ðµ?????Ï¸????È¡??À´
#scRNAsub <-subset(scRNA, cells=cells.use)
##?é¿´Ò»???Ð¶?????????Ï¸??
#scRNAsub
##??????????Á½??subset????
cells.use1 <- subset(scRNA@meta.data, celltype=="Monocyte")
###??scRNA1?Ðµ?????Ï¸????È¡??À´
scRNAsub <- subset(scRNA, cells=row.names(cells.use1))
##?é¿´Ò»???Ð¶?????????Ï¸??
#scRNAsub
###???Ç·???scRNAsub??scRNA.imm??Ò»????Ï¸????Ä¿ Ëµ????È¡Ã»?Ð´?

##???Â½?Î¬???????Â½?Î¬???à£º??Îª?Ù¾?????Ï¸??Ö®???????È½?Ð¡?????Ô¾??àº¯??FindClusters()???Æ·Ö±??ÊµÄ²??????????ßµ?resolution = 0.9
##PCA??Î¬????Ö®Ç°??Ò»?? ???ï²»?Ø¸?×¸????
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub, ndims=20, reduction="pca")
pc.num=1:17
##Ï¸?????à£¬??Ö®Ç°??Ò»?? ???ï²»?Ø¸?×¸????
scRNAsub <- FindNeighbors(scRNAsub, dims = 1:17) 
scRNAsub <- FindClusters(scRNAsub, resolution = 0.8)
table(scRNAsub@meta.data$seurat_clusters)
metadata <- scRNAsub@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
#write.csv(cell_cluster,'cell_cluster.csv',row.names = F)
##?????Ô½?Î¬ ????Ö®Ç°??Ò»?? ???ï²»?Ø¸?×¸????
#tSNE
scRNAsub = RunTSNE(scRNAsub, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub, 'tsne')
#write.csv(embed_tsne,'embed_tsne.csv')
plot1 = DimPlot(scRNAsub, reduction = "tsne") 
plot1
ggsave("tSNE.pdf", plot = plot1, width = 8, height = 7)
#UMAP
scRNAsub <- RunUMAP(scRNAsub, dims = pc.num)
embed_umap <- Embeddings(scRNAsub, 'umap')
#write.csv(embed_umap,'embed_umap.csv') 
plot2 = DimPlot(scRNAsub, reduction = "umap") 
plot2
ggsave("UMAP.pdf", plot = plot2, width = 8, height = 7)

#?Ï²?tSNE??UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
plotc
ggsave("tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)

#?ßµÍ·???????Î¢????????
DimPlot(scRNAsub, reduction = "umap",group.by = "celltype")#?ßµÍ·???????????Î¢????????
DimPlot(scRNA, reduction = "umap",split.by = "risk",label = T) + NoLegend()
?DimPlot





####????À´????×¢??  ×¢?ÍµÄ·??????????? Ç°?????????Ð½??? ?????????Í²???×¸????
##Cluster????????---??Êµ??????marker??????????marker????È¥??Õ¾  ???Ö¶?????????×¢??
#diff.wilcox = FindAllMarkers(scRNAsub)
#all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
#top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#write.csv(all.markers, "subcluster/diff_genes_wilcox.csv", row.names = F)
#write.csv(top10, "subcluster/top10_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub,file="scRNAmac.rdata")
DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='tsne',raster=FALSE)
DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='umap',raster=FALSE)

#
select_genes <- c("MDH1")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#
select_genes <- c("CTSD")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
# 
select_genes <- c("FOLR2")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#
select_genes <- c("HSPA6")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
#
select_genes <- c("IFIT1")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
select_genes <- c("MKI67")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
select_genes <- c("S100A9")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
select_genes <- c("SPP1")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
select_genes <- c("TFPI")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
##Ï¸?????Í¼???
library(SingleR)
library(celldex)
##SingleRÏ¸???????????????????Çµ??Ãµ??????ßµ????Ý¿â£¬???Ï´Î²?Ò»??Ó´  ?Ï´????Ãµ??Ëµ????Ý¿?À´×¢?Íµ?
###????Í¬???Ä·??? ??Ê¦???????Ä°Ù¶????Ä¼????Øµ??????Ô¼??Ñ¾????Ãº?Â·?????Ä¼????Ð£???Òª?Â´??Ø·???  ????ref_Mouse_all.Rdata???Ä¼??????Â´??Ë£?????Ö®??Òª???Øµ???????
###??Ö®Ç°Ò²????  ?????Ä¼?????????Á½?Ö·?Ê½ Ò»?????Ö¶?  Ò»?????Ãº?????È¡  ?????????Ãº?????È¡
#load("ref_Monaco_114s.RData")
ref_Monaco <- MonacoImmuneData()
####???Ç½????Øµ????Ý¸?Öµ??refdata
refdata <- ref_Monaco
###????À´????????Ö®Ç°??singler×¢????Ê½??Ò»Ä£Ò»???? ???ï²»??×¸??
testdata <- GetAssayData(scRNAsub, slot="data")
clusters <- scRNAsub@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"celltype_singleR.csv",row.names = F)
###?????Ç½?singler×¢?Íµ???Ï¢???Øµ?metadata??

celltype[1,2]<-c("low-MDH1")
celltype[2,2]<-c("low-MDH1")
celltype[3,2]<-c("low-MDH1")
celltype[4,2]<-c("high-MDH1")
celltype[5,2]<-c("high-MDH1")
celltype[6,2]<-c("low-MDH1")
celltype[7,2]<-c("low-MDH1")
celltype[8,2]<-c("high-MDH1")
celltype[9,2]<-c("high-MDH1")
celltype[10,2]<-c("low-MDH1")
celltype[11,2]<-c("high-MDH1")
celltype[12,2]<-c("low-MDH1")
celltype[13,2]<-c("high-MDH1")
celltype[14,2]<-c("high-MDH1")
celltype[15,2]<-c("high-MDH1")
celltype[16,2]<-c("high-MDH1")
celltype[17,2]<-c("high-MDH1")

scRNAsub@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNAsub@meta.data[which(scRNAsub@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='tsne')
p1
p2 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='umap')
p2
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
p3
p4 = DimPlot(scRNAsub, group.by="seurat_clusters", label=T, label.size=5, reduction='tsne')
p4
p5 = DimPlot(scRNAsub, group.by="seurat_clusters", label=T, label.size=5, reduction='umap')
p5
ggsave("tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("UMAP_celltype.pdf", p2, width=7 ,height=6)
ggsave("celltype.png", p3, width=10 ,height=5)

save(scRNAsub,file="scRNAmac.rdata")
####??????marker??????singler×¢?Í£????Ç»???Ò»?Ö·??????Ç¸??????Ðµ?????Ñ§ÖªÊ¶???????×£???dotplotÀ´?Ö¶?×¢??
####????????Ò²??Ö®Ç°??Ò»?? ????Ò²????×¸???Ë£??????Ðµ??? ?????Ç±È½Ï¾?×¼??Ò»??×¢?Í·???

#####multi-group
data <- as.data.frame(table(scRNA$risk,scRNA$celltype))
colnames(data) <- c("risk","celltype","Freq")
library(dplyr)
df <- data %>% 
  group_by(risk) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = Freq/Total)

df$celltype  <- factor(df$celltype,levels = unique(df$celltype))

#write.csv(df,file = "output/cell_percent.csv",row.names = F,quote = F)

library(ggplot2)
p <- ggplot(df, aes(x = risk, y = Percent, fill = celltype)) +
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 1, width = 0.95) +
  #scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
p

### ??????Ê½
### facet_wrap ???????????Ð²???
ggplot(df,aes(x = risk, y = Percent,fill=risk))+
  geom_bar(stat="identity")+
  facet_wrap(~ celltype, nrow = 2)+
  theme_bw()+
  NoLegend()

### É¸Ñ¡???? ??Ç°????seurat??subset?È½Ï£????â·º?Íº???
data <- subset(df, subset = !CellType %in% c("Mono/Mk Doublets", "pDC", "Mk", "Eryth"))
ggplot(data,aes(x = group, y = Percent,fill=group))+
  geom_bar(stat="identity")+
  facet_wrap(~ CellType, nrow = 2)+
  theme_bw()+
  NoLegend()

ggplot(data,aes(x = CellType, y = Percent,fill=group))+
  geom_bar(stat="identity")+
  facet_wrap(~ group, nrow = 2)+
  theme_bw()+
  NoLegend()


### ???Ó»?Õ¹Ê¾DotPlot ?? Clustered_DotPlot
all_markers <- readRDS(file = "output/Seurat_stim_all_markers.rds")
library(dplyr)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  slice(1:5) %>%#Ñ¡Ç°5??
  ungroup() %>%
  pull(gene) %>%
  unique()

DotPlot(scobj, features = top_markers)
DotPlot(scobj, features = top_markers) + RotatedAxis()
DotPlot(scobj, features = top_markers) + coord_flip()+ RotatedAxis()
DotPlot(scobj, features = top_markers,dot.scale = 3) + coord_flip()+ RotatedAxis()+
  theme(axis.text.y = element_text(size = 8))

scCustomize::Clustered_DotPlot(scobj, features = top_markers)

### ????dotplot
### ×¢??????, ??Òª??????Ê±????+?Å£??áµ¼?Â»?Í¼Ê§?Ü£?????CD8+ T
Idents(scRNA) <- 'celltype'
levels(Idents(scRNA))


dput(levels(Idents(scRNA)))
c("Epithelial_cells", "T_cells", "Endothelial_cells", "Tissue_stem_cells", "Macrophage", 
  "CMP", "B_cell")
### ????Ë³??
Idents(scRNA) <- factor(Idents(scRNA), 
                        levels = c("Epithelial_cells", "T_cells", "Endothelial_cells", "Fibroblasts", "Macrophage", 
                                   "Neutrophils", "B_cell","NK_cell"))
### Ñ¡????Í¼????
markers.to.plot <- c("MDH1")
### ??Í¼,??ÒªÊ±??cols??????????Á½????É«??split.by???Ö¶???
DotPlot(scRNA, features = markers.to.plot, 
        cols = c("blue", "red"), 
        dot.scale = 6, 
        split.by = "risk") +
  RotatedAxis()

DotPlot(scRNA, features = markers.to.plot, 
        cols = c("blue", "red"), 
        dot.scale = 6, 
) +
  RotatedAxis()
###############################提取上皮细胞#########################
rm(list=ls())
load("scRNA.rdata")
####??È¡Ï¸???Ó¼?--??È¡scRNA1?Ðµ?????Ï¸??
##?Ò¹??Úº?Ò²???? ??È¡Ï¸???Ä·?????Á½??
####????Ò»????which????
phe=scRNA@meta.data
###?é¿´Ò»??immune_annotation?Ðµ???Ï¢????Òª??Îª?Ë·??ã¸´??Õ³??
table(phe$celltype)
###??whichÆ¥??immune??È»???Ãµ?Ï¸????ID ????metadata??????
#cells.use <- row.names(scRNA@meta.data)[which(phe$celltype=='Macrophage')]
#scRNA@meta.data[which(scRNA@meta.data$celltype =="Macrophage"),"celltype2"]<-"Mon/Macr"
#scRNA@meta.data[which(scRNA@meta.data$celltype =="Monocyte"),"celltype2"]<-"Mon/Macr"
###??scRNA1?Ðµ?????Ï¸????È¡??À´
#scRNAsub <-subset(scRNA, cells=cells.use)
##?é¿´Ò»???Ð¶?????????Ï¸??
#scRNAsub
##??????????Á½??subset????
cells.use1 <- subset(scRNA@meta.data, celltype=="Epithelial_cells")
###??scRNA1?Ðµ?????Ï¸????È¡??À´
scRNAsub <- subset(scRNA, cells=row.names(cells.use1))
##?é¿´Ò»???Ð¶?????????Ï¸??
#scRNAsub
###???Ç·???scRNAsub??scRNA.imm??Ò»????Ï¸????Ä¿ Ëµ????È¡Ã»?Ð´?

##???Â½?Î¬???????Â½?Î¬???à£º??Îª?Ù¾?????Ï¸??Ö®???????È½?Ð¡?????Ô¾??àº¯??FindClusters()???Æ·Ö±??ÊµÄ²??????????ßµ?resolution = 0.9
##PCA??Î¬????Ö®Ç°??Ò»?? ???ï²»?Ø¸?×¸????
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub, ndims=20, reduction="pca")
pc.num=1:17
##Ï¸?????à£¬??Ö®Ç°??Ò»?? ???ï²»?Ø¸?×¸????
scRNAsub <- FindNeighbors(scRNAsub, dims = 1:17) 
scRNAsub <- FindClusters(scRNAsub, resolution = 0.8)
table(scRNAsub@meta.data$seurat_clusters)
metadata <- scRNAsub@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
#write.csv(cell_cluster,'cell_cluster.csv',row.names = F)
##?????Ô½?Î¬ ????Ö®Ç°??Ò»?? ???ï²»?Ø¸?×¸????
#tSNE
scRNAsub = RunTSNE(scRNAsub, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub, 'tsne')
#write.csv(embed_tsne,'embed_tsne.csv')
plot1 = DimPlot(scRNAsub, reduction = "tsne") 
plot1
ggsave("tSNEepi.pdf", plot = plot1, width = 8, height = 7)
#UMAP
scRNAsub <- RunUMAP(scRNAsub, dims = pc.num)
embed_umap <- Embeddings(scRNAsub, 'umap')
#write.csv(embed_umap,'embed_umap.csv') 
plot2 = DimPlot(scRNAsub, reduction = "umap") 
plot2
ggsave("UMAPepi.pdf", plot = plot2, width = 8, height = 7)

#?Ï²?tSNE??UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
plotc
ggsave("tSNE_UMAPepi.pdf", plot = plotc, width = 10, height = 5)

#?ßµÍ·???????Î¢????????
DimPlot(scRNAsub, reduction = "umap",group.by = "celltype")#?ßµÍ·???????????Î¢????????
DimPlot(scRNA, reduction = "umap",split.by = "risk",label = T) + NoLegend()
?DimPlot





####????À´????×¢??  ×¢?ÍµÄ·??????????? Ç°?????????Ð½??? ?????????Í²???×¸????
##Cluster????????---??Êµ??????marker??????????marker????È¥??Õ¾  ???Ö¶?????????×¢??
#diff.wilcox = FindAllMarkers(scRNAsub)
#all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
#top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#write.csv(all.markers, "subcluster/diff_genes_wilcox.csv", row.names = F)
#write.csv(top10, "subcluster/top10_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub,file="scRNAepi.rdata")
DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='tsne',raster=FALSE)
DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='umap',raster=FALSE)

#
select_genes <- c("MDH1")
p1 <- VlnPlot(scRNAsub, features = select_genes, pt.size=0, group.by="seurat_clusters", ncol=2)
p1
##Ï¸?????Í¼???
library(SingleR)
library(celldex)
##SingleRÏ¸???????????????????Çµ??Ãµ??????ßµ????Ý¿â£¬???Ï´Î²?Ò»??Ó´  ?Ï´????Ãµ??Ëµ????Ý¿?À´×¢?Íµ?
###????Í¬???Ä·??? ??Ê¦???????Ä°Ù¶????Ä¼????Øµ??????Ô¼??Ñ¾????Ãº?Â·?????Ä¼????Ð£???Òª?Â´??Ø·???  ????ref_Mouse_all.Rdata???Ä¼??????Â´??Ë£?????Ö®??Òª???Øµ???????
###??Ö®Ç°Ò²????  ?????Ä¼?????????Á½?Ö·?Ê½ Ò»?????Ö¶?  Ò»?????Ãº?????È¡  ?????????Ãº?????È¡
#load("ref_Monaco_114s.RData")
ref_Monaco <- MonacoImmuneData()
####???Ç½????Øµ????Ý¸?Öµ??refdata
refdata <- ref_Monaco
###????À´????????Ö®Ç°??singler×¢????Ê½??Ò»Ä£Ò»???? ???ï²»??×¸??
testdata <- GetAssayData(scRNAsub, slot="data")
clusters <- scRNAsub@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"celltype_singleR.csv",row.names = F)
###?????Ç½?singler×¢?Íµ???Ï¢???Øµ?metadata??


celltype[1,2]<-c("low_MDH1")
celltype[2,2]<-c("low_MDH1")
celltype[3,2]<-c("high_MDH1")
celltype[4,2]<-c("high_MDH1")
celltype[5,2]<-c("low_MDH1")
celltype[6,2]<-c("low_MDH1")
celltype[7,2]<-c("low_MDH1")
celltype[8,2]<-c("high_MDH1")
celltype[9,2]<-c("high_MDH1")
celltype[10,2]<-c("high_MDH1")
celltype[11,2]<-c("low_MDH1")
celltype[12,2]<-c("high_MDH1")
celltype[13,2]<-c("low_MDH1")
celltype[14,2]<-c("high_MDH1")
celltype[15,2]<-c("low_MDH1")
celltype[16,2]<-c("high_MDH1")



scRNAsub@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNAsub@meta.data[which(scRNAsub@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='tsne')
p1
p2 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='umap')
p2
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
p3
p4 = DimPlot(scRNAsub, group.by="seurat_clusters", label=T, label.size=5, reduction='tsne')
p4
p5 = DimPlot(scRNAsub, group.by="seurat_clusters", label=T, label.size=5, reduction='umap')
p5
ggsave("tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("UMAP_celltype.pdf", p2, width=7 ,height=6)
ggsave("celltype.png", p3, width=10 ,height=5)

save(scRNAsub,file="scRNAepi.rdata")
####??????marker??????singler×¢?Í£????Ç»???Ò»?Ö·??????Ç¸??????Ðµ?????Ñ§ÖªÊ¶???????×£???dotplotÀ´?Ö¶?×¢??
####????????Ò²??Ö®Ç°??Ò»?? ????Ò²????×¸???Ë£??????Ðµ??? ?????Ç±È½Ï¾?×¼??Ò»??×¢?Í·???

#####multi-group
data <- as.data.frame(table(scRNA$risk,scRNA$celltype))
colnames(data) <- c("risk","celltype","Freq")
library(dplyr)
df <- data %>% 
  group_by(risk) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = Freq/Total)

df$celltype  <- factor(df$celltype,levels = unique(df$celltype))

#write.csv(df,file = "output/cell_percent.csv",row.names = F,quote = F)

library(ggplot2)
p <- ggplot(df, aes(x = risk, y = Percent, fill = celltype)) +
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 1, width = 0.95) +
  #scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
p

### ??????Ê½
### facet_wrap ???????????Ð²???
ggplot(df,aes(x = risk, y = Percent,fill=risk))+
  geom_bar(stat="identity")+
  facet_wrap(~ celltype, nrow = 2)+
  theme_bw()+
  NoLegend()

### É¸Ñ¡???? ??Ç°????seurat??subset?È½Ï£????â·º?Íº???
data <- subset(df, subset = !CellType %in% c("Mono/Mk Doublets", "pDC", "Mk", "Eryth"))
ggplot(data,aes(x = group, y = Percent,fill=group))+
  geom_bar(stat="identity")+
  facet_wrap(~ CellType, nrow = 2)+
  theme_bw()+
  NoLegend()

ggplot(data,aes(x = CellType, y = Percent,fill=group))+
  geom_bar(stat="identity")+
  facet_wrap(~ group, nrow = 2)+
  theme_bw()+
  NoLegend()


### ???Ó»?Õ¹Ê¾DotPlot ?? Clustered_DotPlot
all_markers <- readRDS(file = "output/Seurat_stim_all_markers.rds")
library(dplyr)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  slice(1:5) %>%#Ñ¡Ç°5??
  ungroup() %>%
  pull(gene) %>%
  unique()

DotPlot(scobj, features = top_markers)
DotPlot(scobj, features = top_markers) + RotatedAxis()
DotPlot(scobj, features = top_markers) + coord_flip()+ RotatedAxis()
DotPlot(scobj, features = top_markers,dot.scale = 3) + coord_flip()+ RotatedAxis()+
  theme(axis.text.y = element_text(size = 8))

scCustomize::Clustered_DotPlot(scobj, features = top_markers)

### ????dotplot
### ×¢??????, ??Òª??????Ê±????+?Å£??áµ¼?Â»?Í¼Ê§?Ü£?????CD8+ T
Idents(scRNA) <- 'celltype'
levels(Idents(scRNA))


dput(levels(Idents(scRNA)))
c("Epithelial_cells", "T_cells", "Endothelial_cells", "Tissue_stem_cells", "Macrophage", 
  "CMP", "B_cell")
### ????Ë³??
Idents(scRNA) <- factor(Idents(scRNA), 
                        levels = c("Epithelial_cells", "T_cells", "Endothelial_cells", "Fibroblasts", "Macrophage", 
                                   "Neutrophils", "B_cell","NK_cell"))
### Ñ¡????Í¼????
markers.to.plot <- c("MDH1")
### ??Í¼,??ÒªÊ±??cols??????????Á½????É«??split.by???Ö¶???
DotPlot(scRNA, features = markers.to.plot, 
        cols = c("blue", "red"), 
        dot.scale = 6, 
        split.by = "risk") +
  RotatedAxis()

DotPlot(scRNA, features = markers.to.plot, 
        cols = c("blue", "red"), 
        dot.scale = 6, 
) +
  RotatedAxis()
