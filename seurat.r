#################################################################################################################
#################################################################################################################
#################################################################################################################
### Reference R coding for treating scRNA-seq with Seurat R package
### Please find more tutorial at Seurat official web site:
### https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(Seurat)
library(dplyr)
library(Matrix)
library(hdf5r)

### function for reading single-cell expression matrix
readdata <- function(readsdir, projectname, mincells, minfeatures, percentmt){
    print(paste("Reading", readsdir))
    #pbmc.data <- Read10X(data.dir = readsdir)
    pbmc.data <- read.table(file=readsdir,sep="\t")
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = projectname, min.cells = mincells, min.features = minfeatures)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    pbmc <- subset(pbmc, subset = percent.mt < percentmt)
    return(pbmc)
}

### function for reading general 10X Genomics sequencing data from folder 
readfolder <- function(readsdir, projectname, mincells, minfeatures, percentmt){
    print(paste("Reading", readsdir))
    pbmc.data <- Read10X(data.dir = readsdir)
    #pbmc.data <- read.table(file=readsdir,sep="\t")
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = projectname, min.cells = mincells, min.features = minfeatures)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    pbmc <- subset(pbmc, subset = percent.mt < percentmt)
    return(pbmc)
}

### function for reading single-cell expression data stored in HDF5 format
readhdf5 <- function(readsdir, projectname, mincells, minfeatures, percentmt){
    print(paste("Reading", readsdir))
    pbmc.data <- Read10X_h5(readsdir, use.names = TRUE, unique.features = TRUE)
    #pbmc.data <- read.table(file=readsdir,sep="\t")
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = projectname, min.cells = mincells, min.features = minfeatures)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    pbmc <- subset(pbmc, subset = percent.mt < percentmt)
    return(pbmc)
}


###MCF7R1
mcf7r1path <- "/scRNAseq_data/CellRanger_Counts/MCF7R1_Count_outs/filtered_feature_bc_matrix"
###MCF7R2
mcf7r2path <- "/scRNAseq_data/CellRanger_Counts/MCF7R2_Count_outs/filtered_feature_bc_matrix"
###M1R1
m1r1path <- "/scRNAseq_data/CellRanger_Counts/M1R1_Count_outs/filtered_feature_bc_matrix"
###M1R2
m1r2path <- "/scRNAseq_data/CellRanger_Counts/M1R2_Count_outs/filtered_feature_bc_matrix"
###TRR1
trr1path <- "/scRNAseq_data/CellRanger_Counts/TRR1_Count_outs/filtered_feature_bc_matrix"
###TRR2
trr2path <- "/scRNAseq_data/CellRanger_Counts/TRR2_Count_outs/filtered_feature_bc_matrix"

mcf7r1 <- readfolder(mcf7r1path, 'MCF7_R1', 3, 200, 30)
mcf7r2 <- readfolder(mcf7r2path, 'MCF7_R2', 3, 200, 30)
m1r1 <- readfolder(m1r1path, 'MCF7M1_R1', 3, 200, 30)
m1r2 <- readfolder(m1r2path, 'MCF7M1_R2', 3, 200, 30)
trr1 <- readfolder(trr1path, 'MCF7TR_R1', 3, 200, 30)
trr2 <- readfolder(trr2path, 'MCF7TR_R2', 3, 200, 30)

mcf7all <- merge(mcf7r1, y = c(mcf7r2, m1r1, m1r2, trr1, trr2), add.cell.ids = c("MCF7_R1", "MCF7_R2", "MCF7M1_R1", "MCF7M1_R2", "MCF7TR_R1", "MCF7TR_R2"), project = "MCF7 Single cells")
pbmc <- mcf7all

### total gene expression per cell before normalization
hist(colSums(pbmc$RNA@data),
     breaks = 100,
     main = "Total expression before normalization",
     xlab = "Sum of expression")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

### total gene expression per cell after normalization
hist(colSums(pbmc$RNA@data),
     breaks = 100,
     main = "Total expression after normalization",
     xlab = "Sum of expression")

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

### Top 10 of variable genes
top10 <- head(VariableFeatures(pbmc), 10) #head(pbmc$RNA@var.features,10)
# "PPBP"   "LYZ"    "S100A9" "IGLL5"  "GNLY"   "FTL"    "PF4"    "FTH1"   "GNG11"  "S100A8"

###varible genes with vst.variance
head(pbmc@assays$RNA@meta.features)
vargene <- pbmc@assays$RNA@meta.features
vargene[top10,]
vargene <- vargene[order(-vargene$vst.variance.standardized),]
vargene <- vargene[which(vargene$vst.variable==TRUE),]
vargene$gene <- rownames(vargene)

write.table(vargene[,c("gene", "vst.variance.standardized")],file="vargene.rnk",row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

### plot for variable genes
plot1 <- VariableFeaturePlot(pbmc)
### plot for variable genes with top10 labelled
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

### scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

### remove low quality data identified by scale data
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

### PCA analysis
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

### print major eigenvectors
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

### visualization of eigenvectors

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

### deomonstrate two eigenvectors

DimPlot(pbmc, reduction = "pca",split.by = 'ident')

### show heatmap of eigenvectors

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

### show heatmap of more eigenvectors
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

### dimension of data
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

### cluster of single cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 1)

### look for clusters of cells
head(Idents(pbmc), 5)

### Save old identity classes (the cluster labels) for reference.
pbmc[["old.ident"]] <- Idents(object = pbmc)

# Rename clusters
pbmc <- RenameIdents(object = pbmc, `0` = "D1", `1` = "D2", `2` = "D3", `3` = "D4", `4` = "D5", `5` = "D6", `6` = "D7", `7` = "D8", `8` = "D9", `9` = "D10", `10` = "D11", `11` = "D12", `12` = "D13")

### PCA visualization
pbmc <- RunPCA(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "pca")

### PCA with cell cluster label
DimPlot(pbmc, reduction = "pca",label = TRUE)
LabelClusters(DimPlot(pbmc, reduction = "pca"),id = 'ident')

### UMAP visualization
### install UMAP： reticulate::py_install(packages ='umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

### UMAP with cell cluster label
DimPlot(pbmc, reduction = "umap",label = TRUE)
LabelClusters(DimPlot(pbmc, reduction = "umap"),id = 'ident')

### TSNE visualization
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")

### TSNE with cell cluster label
DimPlot(pbmc, reduction = "tsne",label = TRUE)
LabelClusters(DimPlot(pbmc, reduction = "tsne"),id = 'ident')

###Get the tsne coordinates
tsne1 <- pbmc@reductions$tsne@cell.embeddings[,1]
tsne2 <- pbmc@reductions$tsne@cell.embeddings[,2]

###Pull number of cells in cluster from seurat object
cluster <- pbmc@meta.data

table(cluster$seurat_clusters)

#cluster$clusterlabel <- paste0("N", as.numeric(cluster$seurat_clusters))
cluster$clusterlabel <- paste0("FF", as.numeric(cluster$seurat_clusters))
table(cluster$clusterlabel)
write.table(cluster,file="scRNAseqlusters.txt", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(table(cluster$seurat_clusters),file="clusterNumber.txt",row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

celltype <- table(cluster$orig.ident, cluster$seurat_clusters)
print(celltype)
write.table(celltype,file="clusterCellType.txt",row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

### find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
pbmc.markers$gene <- rownames(pbmc.markers)

write.table(pbmc.markers, file="clusterDEGs.txt",row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

### visualization of interested genes
VlnPlot(pbmc, features = c("ESR1", "EGFR"))
DotPlot(pbmc, features = c("ESR1", "EGFR"))
FeaturePlot(pbmc, features = c("ESR1", "EGFR"))
