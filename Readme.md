# Single-Cell RNA-seq Analysis of 1K PBMCs with R (Seurat) & Python (Scanpy) 

# Introduction

Single-cell RNA sequencing (scRNA-seq) of peripheral blood mononuclear cells (PBMCs) is a powerful approach to profile immune heterogeneity across diverse biological contexts, including infection, autoimmunity, aging, and neurodegeneration. scRNA-seq enables high-resolution profiling of gene expression at the individual cell level, allowing the identification of distinct cell populations within heterogeneous tissues. In this project, we analyze a dataset consisting of 1k peripheral blood mononuclear cell (PBMCs) from a Healthy Donor and are freely available from 10x Genomics using both R (Seurat) and Python (Scanpy) to ensure reproducibility, cross-platform validation, and methodological comparison.
The workflow spans raw data preprocessing, quality control, dimensionality reduction, clustering, marker gene identification, and cell type annotation using both automated reference-based approaches and manual marker-based annotation.

# Project Objectives

The main objectives of this project are to:
- Preprocess raw scRNA-seq data using Kallisto-BUStools
- Construct and analyze a Seurat object in R
- Perform rigorous quality control (QC) using standard metrics
- Generate cell clusters and visualize them using UMAP
- Identify cluster-specific marker genes
- Annotate cell types using:
 * Reference-based methods (SingleR, CellAssign)
 *  Manual annotation using canonical marker genes in Scanpy
- Replicate and validate the analysis pipeline in Python (Scanpy)

# Data Preprocessing (Kallisto-BUStools)

Raw FASTQ files were preprocessed using kb-python (Kallisto-BUStools) to generate gene-by-cell count matrices. Pseudoalignment of reads using kallisto. UMI correction and count matrix generation using bustools

#### Reference Preparation
A human reference transcriptome was prepared using GRCh38 cDNA from Ensembl. The reference index and transcript-to-gene mapping were generated as follows:
```
kb ref \
  -d human \
  -i Homo_sapiens.GRCh38.cdna.all.index \
  -g t2g.txt

```
Where:

-d human specifies the organism
-i defines the output kallisto index
-g generates a transcript-to-gene mapping file required for gene-level quantification

#### UMI Counting and Matrix Generation
Paired-end FASTQ files were processed using the 10x Genomics v3 chemistry. The command below performs pseudoalignment, UMI correction, and count matrix generation in Cell Ranger–compatible format:
```
kb count \
  file1.fastq.gz file2.fastq.gz \
  -i Homo_sapiens.GRCh38.cdna.all.index \
  -g t2g.txt \
  -x 10XV3 \
  -t 8 \
  --cellranger
```

Where:
-x 10XV3 specifies 10x Genomics v3 chemistry
-t 8 sets the number of processing threads
--cellranger outputs matrices in Cell Ranger–compatible format (matrix.mtx, genes.tsv, barcodes.tsv)

#### Output Files
This step produces the following files for downstream analysis in Seurat (R) and Scanpy (Python):
- matrix.mtx – sparse UMI count matrix
- genes.tsv / features.tsv – gene identifiers and symbols
- barcodes.tsv – cell barcodes

# R-Based Analysis (Seurat)
Data Import: Preprocessed count matrices were imported into R and converted into a Seurat object.
Quality Control (QC): performed using DropUtils and Seurat’s built-in metrics, including:
-Number of detected genes per cell (nFeature_RNA)
-Total UMI counts per cell (nCount_RNA)
-Percentage of mitochondrial gene expression (percent.mt)
![metrics Violin plots](https://github.com/DOREENKDAVID/Single-Cell-RNA-seq-Analysis-of-1K-PBMCs-Data-/blob/main/figures/QC_violin_plot.png)

A comprehensive QC report was generated using diagnostic plots such as:
![QC report](https://github.com/DOREENKDAVID/Single-Cell-RNA-seq-Analysis-of-1K-PBMCs-Data-/blob/main/figures/barcode_rank.png)

### Normalization and Feature Selection:
Data normalization using NormalizeData, Identification of highly variable genes using FindVariableFeatures, Data scaling to ensure equal contribution of genes
![normalization](https://github.com/DOREENKDAVID/Single-Cell-RNA-seq-Analysis-of-1K-PBMCs-Data-/blob/main/figures/normalization.png)
![FindVariableFeatures](https://github.com/DOREENKDAVID/Single-Cell-RNA-seq-Analysis-of-1K-PBMCs-Data-/blob/main/figures/HVG.png)



### Dimensionality Reduction and Clustering: 
Principal Component Analysis (PCA), Nearest-neighbor graph construction, Clustering using the Louvain/Leiden algorithm
![PCA](https://github.com/DOREENKDAVID/Single-Cell-RNA-seq-Analysis-of-1K-PBMCs-Data-/blob/main/figures/pca_plots.png)
### Marker Gene Identification : 
Cluster-specific marker genes were identified using: FindAllMarkers, Statistical tests including Wilcoxon rank-sum. Top marker genes were visualized using heatmaps and feature plots.

### Automated Cell Type Annotation: Clusters were annotated using:
- SingleR with public immune reference datasets (e.g., Monaco Immune Data)
- CellAssign for probabilistic cell type assignment
- cluster annotation in Scanpy was performed manually, based on: Canonical immune marker genes (e.g., CD3D, MS4A1, NKG7, FCGR3A) 
- major cell types: CD4 T cells, CD8 T cells, B cells, NK cells, CD14+ Monocytes, FCGR3A+ Monocytes, Dendritic cells
![cluster annotation](https://github.com/DOREENKDAVID/Single-Cell-RNA-seq-Analysis-of-1K-PBMCs-Data-/blob/main/figures/cluster_annotation_plot.png)
# Results and Interpretation
Consistent immune cell populations were identified across both R and Python pipelines
Marker genes matched known PBMC cell-type signatures
major cell types: CD4 T cells, CD8 T cells, B cells, NK cells, CD14+ Monocytes, FCGR3A+ Monocytes, Dendritic cells were present

# Tools and Technologies
Kallisto-BUStools (kb-python)
R: Seurat, DropUtils, SingleR, CellAssign
Python: Scanpy, AnnData, CellTypist
Visualization: UMAP, heatmaps, violin plots, dot plots

# Conclusion
This project demonstrates a full scRNA-seq analysis pipeline from raw data preprocessing to biological interpretation, using both R and Python ecosystems. By combining automated reference-based annotation with manual marker-driven interpretation, the analysis ensures robust and biologically meaningful cell type identification.
