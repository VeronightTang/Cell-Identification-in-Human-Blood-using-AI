# Cell-Identification-in-Human-Blood-using-AI
Leiden clustering and deconvolution with AutoGeneS optimized NNLS

Project Overview - Single Cell Analysis
This project involves single-cell RNA sequencing (scRNA-seq) data analysis, specifically focusing on the deconvolution of bulk RNA sequencing (bulk-seq) data using various Python scripts and tools. Below is a detailed description of the purpose and function of each file in the project directory.

File Descriptions
1. Bulk_data_for_generate.py
Purpose: This script is used to generate bulk RNA-seq data from single-cell RNA-seq data. The generated bulk data is simulated and is later used in the deconvolution process to estimate cell type proportions.
2. Divide_data.py
Purpose: This script first splits the Marker_genes.h5ad file into two datasets for subsequent analysis, specifically for clustering and deconvolution. Additionally, it further divides the generated datasets into training and validation sets for the deconvolution analysis, ensuring that the data is appropriately partitioned to support model training and validation.
3. picture_for_distribution.py
Purpose: This script generates visualizations, specifically for analyzing the distribution of different cell types across the dataset. The visualizations help in understanding the data structure before proceeding with the deconvolution analysis.
4. proportion_for_cells.py
Purpose: This script calculates the proportions of each cell type within the single-cell RNA-seq data. These proportions are then used as reference values during the deconvolution process to compare against estimated values from bulk RNA-seq data.
5. AutogeneS.py
Purpose: This script contains the implementation of the AutoGeneS-optimized NNLS algorithm. AutoGeneS is used to optimize the selection of marker genes, enhancing the accuracy of deconvolution and improving the performance of the NNLS algorithm.
6. Matplot.py
Purpose: This script is used to compare the deconvolution results with the validation set and visualize them in the form of bar charts. It is a key tool for analyzing and presenting the results, primarily using the Matplotlib library for generating plots.
7. nnls.py
Purpose: This script implements the traditional Non-Negative Least Squares (NNLS) algorithm for deconvolution. It estimates the proportions of different cell types within the bulk RNA-seq data based on the single-cell reference data.
8. Leiden for scanpy.py
Purpose: This script applies the Leiden algorithm for clustering the single-cell RNA-seq data. It is used in conjunction with the Scanpy library to identify distinct cell types within the dataset, which is essential for accurate deconvolution.
9. Simulated_bulk_data.csv
Purpose: This file contains the simulated bulk RNA-seq data generated from the single-cell RNA-seq data. It serves as the input data for the deconvolution process.
10. Marker_genes.h5ad
Purpose: This file is an .h5ad file that contains single-cell RNA-seq data, including the expression levels of various marker genes. It is used as a reference dataset for clustering, marker gene selection, and deconvolution.
