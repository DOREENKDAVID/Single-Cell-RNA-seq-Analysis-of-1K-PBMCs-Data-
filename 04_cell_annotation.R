#-----------------------------------------------------

# Assigning identity to cell clusters  ----

#------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "SingleR",
  "celldex",
  "TFBSTools",
  "BSgenome.Hsapiens.UCSC.hg38"
))

install.packages("pheatmap")
install.packages("devtools")
devtools::install_github("satijalab/azimuth")


library(SingleR)      # automated cell-type annotation
library(celldex)      # curated reference datasets
library(pheatmap)     # heatmap visualization
library(Azimuth)      # scRNA-seq atlas-based annotation
library(SingleCellExperiment)
library(Seurat)

#-----------------------------------------------------

# Convert Seurat → SingleCellExperiment

#-----------------------------------------------------
#SingleR works natively with SingleCellExperiment (SCE) objects.


# option 1 - use 'read10xCounts' function from DropletUtils package read directly from 10x output
pbmc_sce <- read10xCounts(datadir)

# option 2 - convert from Seurat
pbmc_sce <- as.SingleCellExperiment(pbmc.1k.seurat)

# the singleCellExperiment data structure is easy to work with
rownames(pbmc_sce)        # genes
colnames(pbmc_sce)        # cells
assays(pbmc_sce)          # count / logcounts matrices
reducedDims(pbmc_sce)     # PCA / UMAP embeddings

assays(pbmc_sce)
subset <- pbmc_sce[,c(1,2,8)]
rowData(pbmc_sce)$gene_symbol <- rownames(pbmc_sce)

#------------------------------------------

#Cell-type annotation using SingleR

#---------------------------------------------
#SingleR compares your cells against reference expression profiles.

#Load reference datasets
monaco_ref  <- MonacoImmuneData(ensembl = FALSE)   # immune-focused (human)
hemato_ref  <- NovershternHematopoieticData(ensembl = FALSE)


# immune reference
singleR_monaco <- SingleR(
  test  = pbmc_sce,
  ref   = monaco_ref,
  labels = monaco_ref$label.main
)
#visualize annotation
plotScoreHeatmap(singleR_monaco)

#hematopoietic reference
singleR_hemato <- SingleR(
  test  = pbmc_sce,
  ref   = hemato_ref,
  labels = hemato_ref$label.main
)


plotScoreHeatmap(singleR_hemato)


#add predicted labels to objects
pbmc_sce$SingleR_label_monaco <- singleR_monaco$labels
pbmc_sce$SingleR_label_hemato <- singleR_hemato$labels


#Visualize SingleR annotations on UMAP
plotUMAP(pbmc_sce, colour_by = "SingleR_label_monaco")
plotUMAP(pbmc_sce, colour_by = "SingleR_label_hemato")








#Cell-type annotation using Azimuth

#Azimuth performs single-cell to single-cell mapping using curated atlases.

#Run Azimuth PBMC reference

pbmc.1k.seurat <- RunAzimuth(
  pbmc.1k.seurat,
  reference = "pbmcref",
  query.modality = "RNA",
  umap.name = "azimuth.umap"
)
DimPlot(
  pbmc.1k.seurat,
  reduction = "umap",
  group.by  = "predicted.celltype.l1",
  label     = TRUE
) + NoLegend()


FeaturePlot(pbmc.1k.seurat, features = c("CD3D", "MS4A1", "NKG7"))

pbmc.1k.sce[["SingleR.labels"]] <- predictions_2$labels
plotUMAP(pbmc.1k.sce, colour_by = "SingleR.labels")
# OPTION 2: assign identity to cell clusters using Azimuth
# Azimuth using reference 'atlas' scRNA-seq datasets to ID clusters in a query dataset
# Azimuth lets you choose different levels of specificity for cell type annotation
# The RunAzimuth function can take a Seurat object as input
pbmc.1k.seurat <- RunAzimuth(pbmc.1k.seurat, reference = "pbmcref", 
                             query.modality = "RNA", 
                             umap.name = "azimuth.umap")
# visualize UMAP with Azimuth labels
DimPlot(pbmc.1k.seurat, reduction = "umap", group.by  = "predicted.celltype.l1", label = TRUE) + NoLegend()

