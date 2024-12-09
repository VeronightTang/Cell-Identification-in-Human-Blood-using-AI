import pandas as pd
import scanpy as sc
import numpy as np
import anndata

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")

# # 指定结果文件
results_file = "Marker_genes.h5ad"  # 用于存储分析结果的文件

# 读取10X Genomics数据并创建AnnData对象
# adata = sc.read_10x_mtx(
#     "68kPBMC/hg19/",  # 包含`.mtx`文件的目录
#     var_names="gene_symbols",  # 使用基因符号作为变量名称
#     cache=True,  # 写入缓存文件以便后续快速读取
# )

# 读取所有 .mtx 文件并合并成一个 AnnData 对象
data_dirs = [
    "gene_bc_matrices/CD4-T cell sample1/",
    "gene_bc_matrices/CD4-T cell sample2/",
    "gene_bc_matrices/CD4-T cell sample3/",
    "gene_bc_matrices/CD8-T cell sample1/",
    "gene_bc_matrices/NKT cells sample1/",
    "gene_bc_matrices/PBMC sample1/",
    "gene_bc_matrices/PBMC sample2/",
    "gene_bc_matrices/PBMC sample3/",
    "gene_bc_matrices/PBMC sample4/",
]

adatas = []
for data_dir in data_dirs:
    adatar = sc.read_10x_mtx(
        data_dir,  # 包含 `.mtx` 文件的目录
        var_names="gene_symbols",  # 使用基因符号作为变量名称
        cache=True,  # 写入缓存文件以便后续快速读取
    )
    adatas.append(adatar)

# 合并所有 AnnData 对象
adata = anndata.concat(adatas, label="batch", keys=data_dirs)

# 确保基因名称唯一
adata.var_names_make_unique()

# 显示每个细胞中表达最高的基因
sc.pl.highest_expr_genes(adata, n_top=20)

# 基本过滤细胞
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# 计算质量度量的小提琴图
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

# 去除线粒体基因表达过多或总数过多的细胞
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

# 归一化和对数化数据
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 寻找高变异基因
sc.pp.highly_variable_genes(adata, min_mean=0.05, max_mean=5, min_disp=0.35)
sc.pl.highly_variable_genes(adata)

adata.raw = adata

# 回归出总计数和线粒体基因百分比的影响
adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

# 运行PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color="CST3")
sc.pl.pca_variance_ratio(adata, log=True)
adata.write(results_file)

# 计算邻近图
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# 计算UMAP
sc.tl.umap(adata)
sc.pl.umap(adata, color=["CST3", "NKG7", "PPBP"])

# 聚类分析
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=["leiden", "CST3", "NKG7"])
adata.write(results_file)

# 计算每个簇中高度差异基因的排名
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
adata.write(results_file)

adata = sc.read(results_file)

# 在数据帧中显示每个簇 0、1、...、7 排名靠前的 10 个基因
pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(5)

# 获取包含分数和组的表格
result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
pd.DataFrame(
    {
        group + "_" + key[:1]: result[key][group]
        for group in groups
        for key in ["names", "pvals"]
    }
).head(5)

# 与单个集群比较
sc.tl.rank_genes_groups(adata, "leiden", groups=["0"], reference="1", method="wilcoxon")
sc.pl.rank_genes_groups(adata, groups=["0"], n_genes=20)
sc.pl.rank_genes_groups_violin(adata, groups="0", n_genes=8)

# 跨组比较某个基因
sc.pl.violin(adata, ["CST3", "NKG7", "PPBP"], groupby="leiden")

# 实际标记细胞类型
new_cluster_names = [
    "CD4 T",
    "Plasma",
    "FCGR3A+ Monocytes",
    "CD56^dimCD16+ NK",
    "CD56^brightCD16- NK",
    "CD8 T",
    "CD14+ Monocytes",
    "CD16+ Monocytes",
    "Conventional Dendritic",
    "Plasmacytoid Dendritic",
    "Megakaryocytes",
    "Mature B",
    "Cytotoxic T lymphocytes",
    "Macrophages",
    "Neutrophils",
]
adata.rename_categories("leiden", new_cluster_names)

sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="right margin",  # 将图例放置在右边缘
    legend_fontsize=10,  # 调整图例字体大小
    title="",
    frameon=False,
    save=".pdf"
)

# 可视化标记基因
# 定义一个标记基因列表
marker_genes = [
    *["IL7R", "CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "CD14"],
    *["LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1"],
    *["FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"],
]

sc.pl.dotplot(adata, marker_genes, groupby="leiden")
sc.pl.stacked_violin(adata, marker_genes, groupby="leiden")

# `compression='gzip'` saves disk space, and slightly slows down writing and subsequent reading
adata.write(results_file, compression="gzip")