import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# 读取已经生成的 h5ad 文件
adata = sc.read("Marker_genes.h5ad")

# 确保细胞类型信息存在于adata.obs中（列名为'leiden'）
if 'leiden' not in adata.obs:
    raise ValueError("细胞类型信息不存在于adata.obs中，请确保已进行细胞类型注释")

# 获取细胞类型信息
cell_types = adata.obs['leiden']

# 计算每种细胞类型的频率
cell_type_counts = cell_types.value_counts()

# 将数据转换为DataFrame以便绘图
cell_type_df = pd.DataFrame({'Cell Type': cell_type_counts.index, 'Count': cell_type_counts.values})

# 绘制直方图
plt.figure(figsize=(10, 6))
sns.barplot(x='Cell Type', y='Count', data=cell_type_df)
plt.xticks(rotation=90)
plt.xlabel('Cell Type')
plt.ylabel('Count')

# 调整底部空间以显示完整的X轴标签
plt.subplots_adjust(bottom=0.35)

# 保存为 PDF 文件
plt.savefig('Cell Type Distribution.pdf', format='pdf')
plt.show()
