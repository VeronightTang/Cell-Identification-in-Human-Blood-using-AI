import anndata
import numpy as np
import pandas as pd
import autogenes as ag
import scanpy as sc
import matplotlib.pyplot as plt

bulk_data = pd.read_csv('Simulated_bulk_data.csv', index_col=0)
adata = sc.read('Marker_genes_part1.h5ad', cache=True)

# # 标准化单细胞数据
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
#
# # 选择最变异的基因
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# adata = adata[:, adata.var.highly_variable]

# 确保基因名称一致
common_genes = adata.var_names.intersection(bulk_data.index)
adata_sc = adata[:, common_genes]
bulk_data = bulk_data.loc[common_genes]

ag.init(adata,use_highly_variable=True,celltype_key='leiden')
ag.optimize(ngen=5000,nfeatures=400,seed=0,mode='fixed')
ag.plot(weights=(-1,0))
# 保存图像为PDF文件
plt.savefig('correlation_distance_plot.pdf', format='pdf')
index = ag.select(index=0)
results = ag.deconvolve(bulk_data.T, model='nnls')
# nusvr linear nnls
# 将结果转换为DataFrame并保存到文件
results_df = pd.DataFrame(results, index=bulk_data.columns, columns=adata_sc.obs['leiden'].cat.categories)
# results_df.to_csv('autogenes_results_nusvr.csv')
# results_df.to_csv('autogenes_results_linear.csv')
results_df.to_csv('autogenes_results_nnls.csv')

print(results_df)
