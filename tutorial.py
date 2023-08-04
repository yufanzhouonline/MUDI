###############################################
### schic cluster and plot

from cluster import *
from plot import *

### read configuration file
mudi = MUDI('config.json')

### instance of class Cluster
cluster = Cluster(mudi.config)

### down sampling cells
cluster.down_sampling(frac = 0.5)
cluster.down_sampling(frac = 0.1)
cluster.down_sampling(frac = 0.75)
cluster.down_sampling(frac = 0.25)

### read chromatin contacts from files
cluster.read_contact()

### preprocessing of contacts
cluster.preprocessing()

### call schicluster to get PCA
cluster.call_pca()

### plot 2D for all cells
Plot.all2D(data=cluster.RAW_PCA, cell=cluster.SAMPLE_COUNT)

### plot 3D for all cells
Plot.all3D(data=cluster.RAW_PCA, cell=cluster.SAMPLE_COUNT)

### determine the optimal cluster number by Silhouette Coefficient
Plot.cluster_number(X = cluster.RAW_PCA[1][:,0:3])

### define clusters by GaussianMixture, n_cluster is the number of clusters
cluster.call_cluster(X = cluster.RAW_PCA[1][:,0:3], n_cluster = 14)
cluster.call_cluster(X = cluster.RAW_PCA[1][:,0:3], n_cluster = 9)

### if pca has been saved, read from file raw_pca_1.txt, raw_pca_2.txt
cluster.read_pca()

### if cluster has been saved, read from file pdcluster.txt
cluster.read_cluster()

### plot 2D cluster
Plot.cluster2D(X = cluster.RAW_PCA[1][:,0:3], yhat=cluster.CLUSTER_LABEL, pc1=0, pc2=1)

### plot 3D cluster
Plot.cluster3D(X = cluster.RAW_PCA[1][:,0:3], yhat=cluster.CLUSTER_LABEL)
###################################
### scRNAseq cluster and plot
### please read R source code in seurat.r as reference

#########################################################
from scintegrator import *

### instance of class Integrator
cell = Integrator(mudi.config)

### load file of cell labeling
cell.cell_label()

### show cell_id, sample_type, scHi-C data file path
cell.CELL

### if down sampling has performed
cell.CELL = cluster.CELL.copy()

### read cell cluster file
cell.cell_cluster()

### show cell cluster: sample, cell, cluster, PC1, PC2, PC3, path, label
cell.PDCLUSTER

### list of sample name in each cluster
cell.CLUSTER_LIST

### number of scHi-C clusters
cell.CLUSTER_NUM

### return the number of CADs in cluster indicated
cell.common_binnum(sample_name = cell.CELL, cluster_list = cell.CLUSTER_LIST[0], binvalue_cutoff = 10)

### Return CADs number per cell of all clusters
### cell.common_bin_per_cell()

### Return the set of CADs, each set is aaaabbbb, aaaa is bin1, bbbb is bin2 
cell.common_binset(sample_name = cell.CELL, cluster_list = cell.CLUSTER_LIST[0], binvalue_cutoff = 0)

### Calculate the score of CADs
cell.cads_score(sample_name = cell.CELL, cluster_list = cell.CLUSTER_LIST[0], binvalue_cutoff = 0)

### list of interaction frequencies of each cell in a cluster
cell.BIN_LIST

### data frame of interaction frequencies, common bins (# of cells) of a cluster
cell.BIN_SCORE

### Return the list of genes in CADs of each cluster to cell.GENE_LIST
#cell.cluster_gene_list(cutoff_list = [200] * 14)
cell.cluster_gene_list(cutoff_list = [70, 50, 120, 60, 80, 60, 80, 70, 50, 50, 60, 50, 100, 70])
#cell.cluster_gene_list(cutoff_list = [70, 50, 120, 60, 80, 60, 80, 70, 50, 50, 60, 50, 100, 70])
cell.cluster_gene_list(cutoff_list = [6, 0, 6, 0, 6, 0, 6, 0, 0])
cell.cluster_gene_list(cutoff_list = [0] * cell.CLUSTER_NUM)
