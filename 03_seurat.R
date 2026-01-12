#-----------------------------------------------
# STEP 1: Load libraries and configure environment
#-----------------------------------------------
BiocManager::install("DropletUtils")
# Core single-cell analysis
library(Seurat)              # Main scRNA-seq toolkit
library(SeuratObject)        # Seurat object infrastructure

# QC-specific packages
library(DropletUtils)        # Empty droplet detection
library(scDblFinder)         # Doublet detection
library(SoupX)               # Ambient RNA correction
library(SingleCellExperiment) # SCE objects

# Visualization
library(ggplot2)             # Custom plots when needed
library(patchwork)           # Combine plots
library(dplyr)               # Data manipulatio
# Load other packages ----
library(tidyverse)
library(scater) # for quality control and visualization for scRNA-seq data
library(scran) # for low level processing of scRNA-seq data
library(Matrix)
library(scales)
library(rjson)
library(R2HTML)
library(DT)

#----------------------------------------------------

# Create Seurat object

#----------------------------------------------------
# if you're working with data from CellRanger, you will still use the Read10X function below to read in your filtered feature_bc_matrix file
# OPTION A: Load RAW matrix 
#counts <- Read10X(data.dir = file.path(cellranger_output, "raw_feature_bc_matrix"))

# OPTION B: Load FILTERED matrix 

datadir <- 'counts_filtered'
list.files(datadir)

expression_matrix <- Read10X(
  datadir,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

# creating the Seurat Object
pbmc.1k.seurat <- CreateSeuratObject(
  counts = expression_matrix,
  min.cells = 3 #filter cells with atleast 3 genes
) 
# If you already have a Seurat object, you can load it directly
#mySeuratObject <- LoadSeuratRds("mySeuratObject.rds")

# normalization Scales gene counts per cell to account for different sequencing depths
# find variable features: Identifies genes that vary most across cells Keeps top 2000 informative genes
#scalling to mean =0 std =1 Prevents high-expression genes from dominating PCA
#PCA Captures major biological variation, dimensionality reduction
#find neighbours :Determines which cells are similar
#find clusters : Groups cells into clusters using the neighbor graph. Assigns a cluster label to every cell
pbmc.1k.seurat <- pbmc.1k.seurat %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE)



# --------------------------------------------------
# Step 1: Calculate percentage of mitochondrial reads
# --------------------------------------------------
# NOTE:
# - Human genes: "^MT-"
# - Mouse genes: "^mt-"

pbmc_qc <- pbmc.1k.seurat

pbmc_qc[["percent_mito"]] <- PercentageFeatureSet(
  object = pbmc_qc,
  pattern = "^MT-"
)

# Violin plots showing distribution of:
# - nCount_RNA   : total UMI counts per cell
# - nFeature_RNA: number of genes detected per cell
# - percent_mito: fraction of mitochondrial reads

VlnPlot(
  pbmc_qc,
  features = c("nCount_RNA", "nFeature_RNA", "percent_mito"),
  pt.size = 0.1
)


# Filter cells based on QC thresholds
# These thresholds remove:
# - low-complexity cells (low UMIs / low genes)
# - potential doublets (extremely high UMIs)
# - stressed or dying cells (high mitochondrial content)

pbmc_qc <- subset(
  pbmc_qc,
  subset =
    nCount_RNA > 1000 &
    nCount_RNA < 20000 &
    nFeature_RNA > 1000 &
    percent_mito < 20
)
# NOTE:
# QC thresholds must be chosen carefully to avoid removing
# rare or biologically meaningful cell populations.

library(ggplot2)

ggplot(
  pbmc_qc@meta.data,
  aes(x = nCount_RNA, y = nFeature_RNA)
) +
  geom_point(alpha = 0.7, size = 0.5) +
  labs(
    x = "Total UMI counts per cell",
    y = "Number of genes detected per cell",
    title = "Gene complexity vs sequencing depth"
  ) +
  theme_classic()


FeatureScatter(
  pbmc_qc,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)
# Interpretation of UMI vs gene scatter plot:
# 1. Bottom-left: low UMIs and low genes → low-quality cells
# 2. Bottom-right: high UMIs but low genes →
#    possible dying cells or low-complexity cell types (e.g. RBC erythrocytes)
# 3. Top-right: high-quality cells with good complexity


# Identify top highly variable genes
top_variable_genes <- head(VariableFeatures(pbmc_qc), 10)
top_variable_genes


# List of all genes for scaling
all_genes <- rownames(pbmc_qc)

# Scaling ensures:
# 1. Mean expression per gene = 0
# 2. Variance per gene = 1
# This prevents highly expressed genes from dominating PCA

#----------------plot UMAP----------------------------
pbmc.1k.seurat <- ScaleData(pbmc.1k.seurat, features = all_genes)
pbmc.1k.seurat <- RunPCA(pbmc.1k.seurat, npcs = 40, verbose = FALSE)

## Visualize genes driving variation in the first two PCs
VizDimLoadings(pbmc.1k.seurat, dims = 1:2, reduction = "pca")

# # Elbow plot helps identify PCs capturing meaningful biological signal
# PCs before the "elbow" capture structured biological variation
# PCs after the elbow mostly represent noise
# Here, we retain the first 20 PCs

ElbowPlot(pbmc.1k.seurat) 


# UMAP visualization of single cells based on the top 20 principal components.
pbmc.1k.seurat <- RunUMAP(pbmc.1k.seurat, reduction = "pca", dims = 1:20)

# Build shared nearest neighbor (SNN) graph
pbmc.1k.seurat <- FindNeighbors(pbmc.1k.seurat, reduction = "pca", dims = 1:20)

# Identify transcriptionally distinct cell clusters
pbmc.1k.seurat <- FindClusters(pbmc.1k.seurat, resolution = 0.5)

#visualize clusters
DimPlot(pbmc.1k.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

DimPlot(pbmc.1k.seurat, reduction = "pca", split.by = "orig.ident", label = TRUE)

#------------------------------------------------------
# find cluster markers
#------------------------------------------------------

devtools::install_github("immunogenomics/presto")
library(presto) # install with devtools::install_github("immunogenomics/presto")


# There are three main approaches to identify marker genes in Seurat:
#
# 1. FindMarkers()
#    → Compare ONE cluster vs all others (used here)
#
# 2. FindAllMarkers()
#    → Identify markers for EACH cluster
#
# 3. FindConservedMarkers()
#    → Identify genes conserved across conditions or samples

# Identify marker genes for Cluster 1 compared to all other clusters
# min.pct = minimum fraction of cells expressing the gene in either group

cluster1_markers <- FindMarkers(
  object   = pbmc.1k.seurat,
  ident.1  = 1,
  min.pct  = 0.25
)


# Calculate difference in expression prevalence between groups
cluster1_markers$pct_diff <- 
  cluster1_markers$pct.1 - cluster1_markers$pct.2


cluster1_markers_tbl <- 
  as_tibble(cluster1_markers, rownames = "gene")


#Represents genes most strongly upregulated in Cluster 1 Ranked by average log2 fold change
top_markers_cluster1 <- cluster1_markers_tbl %>%
  arrange(desc(avg_log2FC)) %>%     # rank by effect size
  slice_head(n = 20)                # keep top 20 genes

# make list of genes into an interactive table
datatable(top_markers_cluster1, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: Cluster 1 genes',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

# plot genes of interest on UMAP eg GZMB, IGHM> imunoglobulins bcells
FeaturePlot(pbmc.1k.seurat, 
            reduction = "umap", 
            features = c("IGHM","GZMB"),
            pt.size = 0.4, 
            order = TRUE,
            #split.by = "orig.ident",
            min.cutoff = 'q10',
            label = FALSE)

# now let's try with FindAllMarkers
pbmc.1k.markers <- FindAllMarkers(pbmc.1k.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


#top 10 marker genes for each cluster and plot as a heatmap
top10 <- pbmc.1k.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(pbmc.1k.seurat, features = top10$gene)



