import scanpy as sc
import numpy as np
import pandas as pd

# 加载单细胞数据，假设保存为AnnData对象'sc_data.h5ad'
adata = sc.read_h5ad('Marker_genes.h5ad')

# 通过汇总表达值模拟散装RNA测序数据
bulk_data = adata.X.sum(axis=0)

# 将结果转换为pandas DataFrame
bulk_df = pd.DataFrame(bulk_data.T, index=adata.var_names, columns=['bulk_expression'])

# 将模拟的散装数据保存为CSV文件（可选）
bulk_df.to_csv('Simulated_bulk_data.csv')

# 打印出模拟的散装数据的前几行
print(bulk_df.head())
