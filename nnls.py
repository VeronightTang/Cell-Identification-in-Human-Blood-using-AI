import pandas as pd
import numpy as np
import scanpy as sc
from scipy.optimize import nnls
import matplotlib.pyplot as plt

# 加载并标准化单细胞数据
adata = sc.read_h5ad('Marker_genes_part1.h5ad')

# 加载bulk数据
bulk_data = pd.read_csv('Simulated_bulk_data.csv', index_col=0)

# 确保基因名称一致
common_genes = adata.var_names.intersection(bulk_data.index)
adata = adata[:, common_genes]
bulk_data = bulk_data.loc[common_genes]

# 创建伪散装数据
pseudo_bulk_data = adata.to_df().T.groupby(adata.obs['leiden']).sum().T

# 确保基因顺序一致
pseudo_bulk_data = pseudo_bulk_data.loc[common_genes]

# 准备目标矩阵和权重矩阵
X = pseudo_bulk_data.T
y = bulk_data.values

# 进行NNLS反卷积
nnls_results = []
for i in range(y.shape[1]):
    coef, _ = nnls(X, y[:, i])
    nnls_results.append(coef)
nnls_results = np.array(nnls_results)

# 将结果转换为DataFrame并保存到文件
results_df = pd.DataFrame(nnls_results, columns=adata.obs['leiden'].cat.categories, index=bulk_data.columns)
results_df.to_csv('nnls_results.csv')

# 打印结果
print(results_df)

# 计算真实细胞类型比例
true_proportions = adata.obs['leiden'].value_counts(normalize=True)
true_proportions_df = true_proportions * 100  # 转换为百分比
true_proportions_df = true_proportions_df.reset_index()
true_proportions_df.columns = ['Cell Type', 'True Proportion']

# 计算NNLS估计比例
estimated_proportions = results_df.sum(axis=0)
total_estimated_expression = estimated_proportions.sum()
estimated_proportions_percentage = (estimated_proportions / total_estimated_expression) * 100
estimated_proportions_df = estimated_proportions_percentage.reset_index()
estimated_proportions_df.columns = ['Cell Type', 'Estimated Proportion']

# 合并真实比例和估计比例
comparison_df = pd.merge(true_proportions_df, estimated_proportions_df, on='Cell Type', how='outer').fillna(0)

# 绘制比较图
fig, ax = plt.subplots(figsize=(12, 8))
width = 0.35  # the width of the bars
x = np.arange(len(comparison_df['Cell Type']))

# 绘制柱状图
ax.bar(x - width/2, comparison_df['True Proportion'], width, label='True Proportion')
ax.bar(x + width/2, comparison_df['Estimated Proportion'], width, label='Estimated Proportion')

ax.set_xlabel('Cell Type')
ax.set_ylabel('Proportion (%)')
ax.set_title('Comparison of True and Estimated Cell Type Proportions using NNLS')
ax.set_xticks(x)
ax.set_xticklabels(comparison_df['Cell Type'], rotation=90)
ax.legend()

plt.tight_layout()
# 保存图像为PDF文件
plt.savefig('comparison_proportions_nnls.pdf', format='pdf')
plt.show()
