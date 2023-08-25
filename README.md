## MUDI - Multiomic Data Integration for single-cell sequencing

----------------------------------
MUDI is a python package developed for integration of multiomic data for single cell such as single-cell Hi-C, single-cell RNA-seq, single-cell ATAC-seq, single-cell ChIP-seq, single-cell DNA methylation and so on.

Current version supports integration of scHi-C and scRNA-seq.

### Citation:

Please cite as follows when you use this tool:

MUDI: Multiomic Data Integration for single-cell sequencing. (Yufan Zhou. 2023) Github. https://github.com/yufanzhouonline/MUDI

Thank you.

### Installation Requirements:

Python: 3.7.16 or higher

Numpy: 1.21.6 or higher

Pandas: 1.2.3 or higher

Matplotlib: 3.5.3 or higher

Scikit-learn: 1.0.2 or higher

### Installation

Directly download these package to your local and run as API. Please follow Usage section for the detail.

### Input:

1. config.json: setup the configure file

2. refseq.bed: reference genome bed file

3. cell_1000000_abs.bed: bin bed files of genome, four columns: chrosome, bin start, bin end, bin number

4. sample.txt: files of contact pairs of invididual cells of scHiC, include three columns: cell id, sample type or cell type, path stored contact pairs

5. cluster_degs.txt: files of DEGs of scRNAseq clusters, include p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene

6. cluster_number.txt: cluster number of scRNA-seq, include cluster id and cell number of each cluster

7. var_gene.rnk: Top 2000 variable genes of scRNA-seq

Please refer to examples of these input files in the folder data.

Please refer to codes in seruat.r to get scRNA-seq input files.

### Output:

pdcluster.txt: output of clusters of scHi-C data and eigenvectors, columns include cell_id, sample_type, cluster_id, PC1, PC2, PC3, and path stored contact pairs

### Usage:

Please refer to tutorial:

https://github.com/yufanzhouonline/MUDI/blob/master/tutorial.py


