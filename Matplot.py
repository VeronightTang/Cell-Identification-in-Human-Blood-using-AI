import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# 加载真实数据
adata_real = sc.read_h5ad('Marker_genes_part1.h5ad')

# 提取细胞类型注释
cell_types = adata_real.obs['leiden']

# 计算每种细胞类型的比例
true_proportions = cell_types.value_counts(normalize=True)
true_cell_type_proportions = true_proportions * 100  # 转换为百分比

# 将结果转换为 DataFrame
true_proportions_df = true_cell_type_proportions.to_frame().reset_index()
true_proportions_df.columns = ['Cell Type', 'True Proportion']

# 打印真实比例
print(true_proportions_df)

# 加载反卷积结果
autogenes_results = pd.read_csv('autogenes_results_nnls.csv', index_col=0)

# 计算反卷积结果中每种细胞类型的总估计表达值
estimated_proportions = autogenes_results.sum(axis=0)

# 转化为百分比
total_estimated_expression = estimated_proportions.sum()
estimated_proportions_percentage = (estimated_proportions / total_estimated_expression) * 100

# 将结果转换为 DataFrame
estimated_proportions_df = estimated_proportions_percentage.to_frame().reset_index()
estimated_proportions_df.columns = ['Cell Type', 'Estimated Proportion']

# 打印估计比例
print(estimated_proportions_df)

# 合并真实比例和估计比例
comparison_df = pd.merge(true_proportions_df, estimated_proportions_df, on='Cell Type', how='outer').fillna(0)

# 绘制比较图
fig, ax = plt.subplots(figsize=(12, 8))
width = 0.35  # the width of the bars

# 设置位置
x = np.arange(len(comparison_df['Cell Type']))

# 绘制柱状图
ax.bar(x - width/2, comparison_df['True Proportion'], width, label='True Proportion')
ax.bar(x + width/2, comparison_df['Estimated Proportion'], width, label='Estimated Proportion')

ax.set_xlabel('Cell Type')
ax.set_ylabel('Proportion (%)')
# ax.set_title('Comparison of True and Estimated Cell Type Proportions using nusvr')
# ax.set_title('Comparison of True and Estimated Cell Type Proportions using Linear')
# ax.set_title('Comparison of True and Estimated Cell Type Proportions using nnls')
ax.set_xticks(x)
ax.set_xticklabels(comparison_df['Cell Type'], rotation=90)
ax.legend()

plt.tight_layout()
# 保存图像为PDF文件
plt.savefig('comparison_proportions_nnls.pdf', format='pdf')
plt.show()
