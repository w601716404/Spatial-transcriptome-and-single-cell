# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 15:05:55 2024

@author: box
"""
import backports
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sc.logging.print_versions()
sc.set_figure_params(facecolor="white",figsize=(8,8))
sc.settings.verbosity = 3

#file_path = 'E:\\NEW-pancancer_2\\spatial\\GSE277206_RAW\\V0801\\filtered_feature_bc_matrix.h5'

# 读取数据
#adata = sc.read_10x_h5(file_path)


print(os.getcwd())#查看当前路径
os.chdir("E:\\NEW_pancancer_2\\spatial\\GSE203612\\BRCA1")
print(os.getcwd())#检查路径是否成功
#该文件下的h5文件
adata = sc.read_visium('E:\\NEW_pancancer_2\\spatial\\GSE203612\\BRCA1\\')

adata
adata.var_names_make_unique()#把基因名变唯一的名字，去除重复基因
#去除线粒体，赋值与VAR中新的值
adata.var["mt"] = adata.var_names.str.startswith("MT-")
#线粒体百分比的计算
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata
#两种复制方法，第一种的话adata1的变化随着adata变化而变化，第二种就是两个数据
#adata1=adata
#adata1=adata.copy()
#质控
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False, bins=40,ax=axs[1],)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False, bins=60, ax=axs[3],)
#过滤
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)
#归一化
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
#降维，聚类
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters" )

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)

sc.pl.umap(adata, color=["total_counts" ], wspace=0.4)
sc.pl.umap(adata, color=[ "n_genes_by_counts"], wspace=0.4)
sc.pl.umap(adata, color=[ "clusters"], wspace=0.4)

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])

sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)
#筛选，1和3cluster 的分布
sc.pl.spatial(adata,
    img_key="hires",color="clusters",
    groups=["1", "3"],crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,size=1.3,)
#marker基因
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")

#两个基因表达量的分布
sc.pl.spatial(adata, img_key="hires", color=["clusters" ])

#sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)
#sc.pl.spatial(adata, img_key="hires", color=["Mdh1", "Cd68"], alpha=0.7)
sc.pl.spatial(adata, img_key="hires", color=["MDH1", "CD68"], alpha=0.7)
#保存数据
os.chdir('E:\\SpatialD\\Spatial\\k3.python\\k3.data')
adata.write("GSE194329.data.anndata.h5ad")


