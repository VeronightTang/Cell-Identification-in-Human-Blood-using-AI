import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 读取数据
# adata = sc.read_10x_mtx(
#     "6kpbmc/hg19/",  # 替换成你的数据目录
#     var_names="gene_symbols",
#     cache=True
# )
adata = sc.read("Marker_genes.h5ad")

# 打乱并分割数据
indices = np.random.permutation(adata.shape[0])
split_index = len(indices) // 2
adata_1 = adata[indices[:split_index]]
adata_2 = adata[indices[split_index:]]

# 保存分割后的数据
adata_1.write("Marker_genes_part1.h5ad")
adata_2.write("Marker_genes_part2.h5ad")

# # 可视化分割后的数据
# plt.figure(figsize=(6, 4))
# plt.bar(['Part 1', 'Part 2'], [adata_1.shape[0], adata_2.shape[0]])
# plt.xlabel('Dataset Part')
# plt.ylabel('Number of Cells')
# plt.title('Number of Cells in Each Part')
# plt.show()