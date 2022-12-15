library(Seurat)
library(dplyr)
library(stringr)

EA001_expr <- Seurat::Read10X('~/Passegue_bioinformatic_data/EA001/analysis/200522_EMMANUELLE_AMELIE_10_MOUSE_10X-EA001-cellranger-count-default/EA001_cellranger_count_outs/filtered_feature_bc_matrix/')
EA002_expr <- Seurat::Read10X('~/Passegue_bioinformatic_data/EA002/analysis/200522_EMMANUELLE_AMELIE_10_MOUSE_10X-EA002-cellranger-count-default/EA002_cellranger_count_outs/filtered_feature_bc_matrix/')
AA003_expr <- Seurat::Read10X('~/Passegue_bioinformatic_data/AA003/analysis/210222_AMELIE_AMELIE_3_MOUSE_10X-AA003-cellranger-count-default/AA003_cellranger_count_outs/filtered_feature_bc_matrix/')
EJ001_expr <- Seurat::Read10X('~/Passegue_bioinformatic_data/EJ001/EJ001_cellranger_count_outs/filtered_feature_bc_matrix/')
EO0035_expr <- Seurat::Read10X('~/Passegue_bioinformatic_data/EO0035/EO035_cellranger_count_outs/filtered_feature_bc_matrix/')
EO0036_expr <- Seurat::Read10X('~/Passegue_bioinformatic_data/EO0036/EO036_cellranger_count_outs/filtered_feature_bc_matrix/')

colnames(EA001_expr) <- paste('EA001_', colnames(EA001_expr), sep='')
EA001 <- CreateSeuratObject(counts = EA001_expr, min.cells = 3, min.features = 200)

colnames(EA002_expr) <- paste('EA002_', colnames(EA002_expr), sep='')
EA002 <- CreateSeuratObject(counts = EA002_expr, min.cells = 3, min.features = 200)

colnames(AA003_expr) <- paste('AA003_', colnames(AA003_expr), sep='')
AA003 <- CreateSeuratObject(counts = AA003_expr, min.cells = 3, min.features = 200)


colnames(EJ001_expr) <- paste('EJ001_', colnames(EJ001_expr), sep='')
EJ001 <- CreateSeuratObject(counts = EJ001_expr, min.cells = 3, min.features = 200)

colnames(EO0035_expr) <- paste('EO0035_', colnames(EO0035_expr), sep='')
EO0035 <- CreateSeuratObject(counts = EO0035_expr, min.cells = 3, min.features = 200)

colnames(EO0036_expr) <- paste('EO0036_', colnames(EO0036_expr), sep='')
EO0036 <- CreateSeuratObject(counts = EO0036_expr, min.cells = 3, min.features = 200)

############# EA001

# QC and selecting cells for further analysis
EA001[["percent.mt"]] <- PercentageFeatureSet(EA001, pattern = "^MT-")
VlnPlot(EA001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(EA001, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(EA001, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
EA001 <- subset(EA001, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA <90000)

# Remove Rpl|Rps|Mrpl|Mrps|mt- genes
Remin_Gene <- setdiff(rownames(EA001), rownames(EA001)[str_detect(rownames(EA001), 'RPL|RPS|MRPL|MRPS|MT-')])
EA001 <- EA001[Remin_Gene,]

# Normalize data
EA001 <- NormalizeData(EA001, normalization.method = "LogNormalize", scale.factor = 10000)
EA001 <- FindVariableFeatures(EA001, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(EA001), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(EA001)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data
EA001 <- ScaleData(EA001, features = rownames(EA001))

# Run PCA
EA001 <- RunPCA(EA001, features = VariableFeatures(object = EA001))
print(EA001[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(EA001, dims = 1:2, reduction = "pca")
DimPlot(EA001, reduction = "pca")
ElbowPlot(EA001,ndims = 50)

EA001 <- FindNeighbors(EA001, dims = 1:30)
EA001 <- FindClusters(EA001, resolution = 0.5)
head(Idents(EA001), 5)

# Run UMAP
EA001 <- RunUMAP(EA001, dims = 1:10)
DimPlot(EA001, reduction = "umap")

saveRDS(EA001, file = "EA001.rds")

# Finding differentially expressed features (cluster biomarkers)
EA001.markers <- FindAllMarkers(EA001, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
EA001.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

VlnPlot(EA001, features = c("Hlf", "Ctla2a"), slot = "counts", log = TRUE)
EA001.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(EA001, features = top10$gene) + NoLegend()


######### EA002

EA002[["percent.mt"]] <- PercentageFeatureSet(EA002, pattern = "^MT-")
VlnPlot(EA002, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(EA002, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(EA002, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
EA002 <- subset(EA002, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA <90000)

# Remove Rpl|Rps|Mrpl|Mrps|mt- genes
Remin_Gene <- setdiff(rownames(EA002), rownames(EA002)[str_detect(rownames(EA002), 'RPL|RPS|MRPL|MRPS|MT-')])
EA002 <- EA002[Remin_Gene,]

# Normalize data
EA002 <- NormalizeData(EA002, normalization.method = "LogNormalize", scale.factor = 10000)
EA002 <- FindVariableFeatures(EA002, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(EA002), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(EA002)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data
EA002 <- ScaleData(EA002, features = rownames(EA002))

# Run PCA
EA002 <- RunPCA(EA002, features = VariableFeatures(object = EA002))
print(EA002[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(EA002, dims = 1:2, reduction = "pca")
DimPlot(EA002, reduction = "pca")
ElbowPlot(EA002,ndims = 50)

EA002 <- FindNeighbors(EA002, dims = 1:30)
EA002 <- FindClusters(EA002, resolution = 0.5)
head(Idents(EA002), 5)

# Run UMAP
EA002 <- RunUMAP(EA002, dims = 1:10)
DimPlot(EA002, reduction = "umap")

saveRDS(EA002, file = "EA002.rds")

# Finding differentially expressed features (cluster biomarkers)
EA002.markers <- FindAllMarkers(EA002, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
EA002.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

VlnPlot(EA002, features = c("S100a8", "Elane"), slot = "counts", log = TRUE)
EA002.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(EA002, features = top10$gene) + NoLegend()

######### AA003

# QC and selecting cells for further analysis
AA003[["percent.mt"]] <- PercentageFeatureSet(AA003, pattern = "^MT-")
VlnPlot(AA003, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(AA003, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AA003, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
AA003 <- subset(AA003, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA < 90000)

# Remove Rpl|Rps|Mrpl|Mrps|mt- genes
Remin_Gene <- setdiff(rownames(AA003), rownames(AA003)[str_detect(rownames(AA003), 'RPL|RPS|MRPL|MRPS|MT-')])
AA003 <- AA003[Remin_Gene,]

# Normalize data
AA003 <- NormalizeData(AA003, normalization.method = "LogNormalize", scale.factor = 10000)
AA003 <- FindVariableFeatures(AA003, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AA003), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AA003)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data
AA003 <- ScaleData(AA003, features = rownames(AA003))

# Run PCA
AA003 <- RunPCA(AA003, features = VariableFeatures(object = AA003))
print(AA003[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(AA003, dims = 1:2, reduction = "pca")
DimPlot(AA003, reduction = "pca")
ElbowPlot(AA003,ndims = 50)

AA003 <- FindNeighbors(AA003, dims = 1:30)
AA003 <- FindClusters(AA003, resolution = 0.5)
head(Idents(AA003), 5)

# Run UMAP
AA003 <- RunUMAP(AA003, dims = 1:10)
DimPlot(AA003, reduction = "umap")

saveRDS(AA003, file = "AA003.rds")

# Finding differentially expressed features (cluster bio-markers)
AA003.markers <- FindAllMarkers(AA003, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AA003.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

VlnPlot(AA003, features = c("Ifitm1", "Wfdc17"), slot = "counts", log = TRUE)
AA003.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(AA003, features = top10$gene) + NoLegend()

####### EJ001

# QC and selecting cells for further analysis
EJ001[["percent.mt"]] <- PercentageFeatureSet(EJ001, pattern = "^MT-")
VlnPlot(EJ001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(EJ001, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(EJ001, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
EJ001 <- subset(EJ001, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA < 90000)

# Remove Rpl|Rps|Mrpl|Mrps|mt- genes
Remin_Gene <- setdiff(rownames(EJ001), rownames(EJ001)[str_detect(rownames(EJ001), 'RPL|RPS|MRPL|MRPS|MT-')])
EJ001 <- EJ001[Remin_Gene,]

# Normalize data
EJ001 <- NormalizeData(EJ001, normalization.method = "LogNormalize", scale.factor = 10000)
EJ001 <- FindVariableFeatures(EJ001, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(EJ001), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(EJ001)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data
EJ001 <- ScaleData(EJ001, features = rownames(EJ001))

# Run PCA
EJ001 <- RunPCA(EJ001, features = VariableFeatures(object = EJ001))
print(EJ001[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(EJ001, dims = 1:2, reduction = "pca")
DimPlot(EJ001, reduction = "pca")
ElbowPlot(EJ001,ndims = 50)

EJ001 <- FindNeighbors(EJ001, dims = 1:30)
EJ001 <- FindClusters(EJ001, resolution = 0.5)
head(Idents(EJ001), 5)

# Run UMAP
EJ001 <- RunUMAP(EJ001, dims = 1:10)
DimPlot(EJ001, reduction = "umap")

saveRDS(EJ001, file = "EJ001.rds")

# Finding differentially expressed features (cluster biomarkers)
EJ001.markers <- FindAllMarkers(EJ001, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
EJ001.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

VlnPlot(EJ001, features = c("Ung", "Mcm3"), slot = "counts", log = TRUE)
EJ001.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(EJ001, features = top10$gene) + NoLegend()


###### EO0035

# QC and selecting cells for further analysis
EO0035[["percent.mt"]] <- PercentageFeatureSet(EO0035, pattern = "^MT-")
VlnPlot(EO0035, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(EO0035, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(EO0035, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
EO0035 <- subset(EO0035, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA < 90000)

# Remove Rpl|Rps|Mrpl|Mrps|mt- genes
Remin_Gene <- setdiff(rownames(EO0035), rownames(EO0035)[str_detect(rownames(EO0035), 'RPL|RPS|MRPL|MRPS|MT-')])
EO0035 <- EO0035[Remin_Gene,]

# Normalize data
EO0035 <- NormalizeData(EO0035, normalization.method = "LogNormalize", scale.factor = 10000)
EO0035 <- FindVariableFeatures(EO0035, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(EO0035), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(EO0035)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data
EO0035 <- ScaleData(EO0035, features = rownames(EO0035))

# Run PCA
EO0035 <- RunPCA(EO0035, features = VariableFeatures(object = EO0035))
print(EO0035[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(EO0035, dims = 1:2, reduction = "pca")
DimPlot(EO0035, reduction = "pca")
ElbowPlot(EO0035,ndims = 50)

EO0035 <- FindNeighbors(EO0035, dims = 1:30)
EO0035 <- FindClusters(EO0035, resolution = 0.5)
head(Idents(EO0035), 5)

# Run UMAP
EO0035 <- RunUMAP(EO0035, dims = 1:10)
DimPlot(EO0035, reduction = "umap")

saveRDS(EO0035, file = "EO0035.rds")

# Finding differentially expressed features (cluster biomarkers)
EO0035.markers <- FindAllMarkers(EO0035, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
EO0035.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

VlnPlot(EO0035, features = c("Wfdc17", "Dntt"), slot = "counts", log = TRUE)
EO0035.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(EO0035, features = top10$gene) + NoLegend()


####### EO0036

# QC and selecting cells for further analysis
EO0036[["percent.mt"]] <- PercentageFeatureSet(EO0036, pattern = "^MT-")
VlnPlot(EO0036, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(EO0036, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(EO0036, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
EO0036 <- subset(EO0036, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000)

# Remove Rpl|Rps|Mrpl|Mrps|mt- genes
Remin_Gene <- setdiff(rownames(EO0036), rownames(EO0036)[str_detect(rownames(EO0036), 'RPL|RPS|MRPL|MRPS|MT-')])
EO0036 <- EO0036[Remin_Gene,]

# Normalize data
EO0036 <- NormalizeData(EO0036, normalization.method = "LogNormalize", scale.factor = 10000)
EO0036 <- FindVariableFeatures(EO0036, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(EO0036), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(EO0036)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data
EO0036 <- ScaleData(EO0036, features = rownames(EO0036))

# Run PCA
EO0036 <- RunPCA(EO0036, features = VariableFeatures(object = EO0036))
print(EO0036[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(EO0036, dims = 1:2, reduction = "pca")
DimPlot(EO0036, reduction = "pca")
ElbowPlot(EO0036,ndims = 50)

EO0036 <- FindNeighbors(EO0036, dims = 1:30)
EO0036 <- FindClusters(EO0036, resolution = 0.5)
head(Idents(EO0036), 5)

# Run UMAP
EO0036 <- RunUMAP(EO0036, dims = 1:10)
DimPlot(EO0036, reduction = "umap")

saveRDS(EO0036, file = "EO0036.rds")

# Finding differentially expressed features (cluster biomarkers)
EO0036.markers <- FindAllMarkers(EO0036, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
EO0036.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

VlnPlot(EO0036, features = c("Wfdc17", "Ccl9"), slot = "counts", log = TRUE)
EO0036.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(EO0036, features = top10$gene) + NoLegend()


######### integrate
features <- SelectIntegrationFeatures(object.list = list(EA001, EA002, AA003, EJ001, EO0035, EO0036))
EA.anchors <- FindIntegrationAnchors(object.list = list(EA001, EA002, AA003, EJ001, EO0035, EO0036), anchor.features = features)
EA.combined <- IntegrateData(anchorset = EA.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(EA.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
EA.combined <- ScaleData(EA.combined)
EA.combined <- RunPCA(EA.combined, npcs = 50)
ElbowPlot(EA.combined,ndims = 50)
EA.combined <- RunUMAP(EA.combined, reduction = "pca", dims = 1:30)
EA.combined <- FindNeighbors(EA.combined, reduction = "pca", dims = 1:30)
EA.combined <- FindClusters(EA.combined, resolution = 0.5)
#saveRDS(EA.combined, 'EA.combined.rds')

