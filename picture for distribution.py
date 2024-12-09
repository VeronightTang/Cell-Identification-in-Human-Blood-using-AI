import scanpy as sc
import matplotlib.pyplot as plt

# 读取生成的.h5ad文件
adata = sc.read_h5ad('Marker_genes.h5ad')

# 计算每个细胞的基因数量
gene_counts_per_cell = (adata.X > 0).sum(axis=1)

# 将结果转换为一维数组
gene_counts_per_cell = gene_counts_per_cell.A1 if hasattr(gene_counts_per_cell, 'A1') else gene_counts_per_cell

# 计算基因数量的最大值，用于设置直方图的范围
max_genes = gene_counts_per_cell.max()

# 设置每个bin的宽度为4
bin_width = 4
bins = range(0, max_genes + bin_width, bin_width)


# 绘制直方图
plt.figure(figsize=(10, 6))
plt.hist(gene_counts_per_cell, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Number of Genes')
plt.ylabel('Number of Cells')
# plt.title('Distribution of Gene Counts per Cell')
plt.grid(True)
plt.savefig('Distribution of Gene Counts per Cell.pdf', format='pdf')
plt.show()
