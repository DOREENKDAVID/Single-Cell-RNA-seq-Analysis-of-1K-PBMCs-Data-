# in this project i aim to do quality assessment (QA) and analysis of single cell RNA-seq data
#Dataset:  (~1000 cell) dataset from human peripheral blood mononuclear cells (PBMCs).from the public datasets on the 10X Genomics website: https://www.10xgenomics.com/resources/datasets

#-----------------------------------------------
# Installation: Required packages for single-sample QC
#-----------------------------------------------

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install Seurat 5
install.packages("Seurat")

# Install SeuratObject (Seurat's data structure package)
install.packages("SeuratObject")

# Install visualization packages
install.packages(c(
  "ggplot2",        # For custom plots when needed
  "patchwork",      # For combining Seurat plots
  "dplyr"           # Data manipulation
))

# Install Bioconductor manager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages for QC
BiocManager::install(c(
  "DropletUtils",           # Empty droplet detection
  "scDblFinder",            # Doublet detection (better than DoubletFinder for Seurat 5)
  "SingleCellExperiment" ,   # Data structure for some QC tools
  "scran",
  "scater"
))

# Install SoupX for ambient RNA correction
install.packages("SoupX")



#-----------------------------------------------
# STEP 1: Load libraries and configure environment
#-----------------------------------------------

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


# Set random seed for reproducibility
set.seed(42)

#setwd
setwd("D:/DIY_transcriptomics/scRNA_seq")


#-----------------------------------------------
# STEP 2: Load 10x Genomics data
#-----------------------------------------------

#load unfiltered data with all droplets including emptyDrops()
# load raw data matrix using the readMM function from the Matrix package
#matrix.mtx: Sparse UMI count matrix
raw_mtx <- readMM('counts_unfiltered/cellranger/matrix.mtx')
head(raw_mtx)

# load genes.tsv:Gene names and IDs 
genes <- read.csv('counts_unfiltered/cellranger/genes.tsv', sep = '\t', header = F)
head(genes)


# add ensemble gene_ids to the data matrix as rownames
rownames(raw_mtx) <- genes[,1] 

# add cell barcodes as column names
#load barcodes.tsv: Cell barcodes
colnames(raw_mtx) <- read.csv('counts_unfiltered/cellranger/barcodes.tsv', sep = '\t', header = F)[,1] 


#-----------------------------------------------
# STEP 3: Empty droplet detection
#-----------------------------------------------
#Three types of droplets:
#Cell-containing droplets: High UMI counts, diverse gene expression
#Empty droplets: Low UMI counts, only ambient RNA
#Damaged cells: Intermediate UMI counts (need filtering later)
# use DropletUtils package to get probability that each barcode is a cell
out <- emptyDrops(raw_mtx) 
head(out)
# set threshold probability for calling a cell if not emptydrop is less than 0.05
keep <- out$FDR <= 0.05 
# use threshold to remove empty drops
keep[is.na(keep)] <- FALSE
filtered_mtx <- raw_mtx[,keep] 

#FDR False Discovery Rate. (EmptyDrops / scRNA-seq), 
#it’s a statistical measure of confidence that a barcode (droplet) contains a real cell rather than just ambient RNA
plot_df <- data.frame(
  total_umi = Matrix::colSums(raw_mtx),
  FDR = out$FDR,
  is_cell = keep
)

# Make classification explicit (important for legend clarity)
plot_df$is_cell <- factor(
  plot_df$is_cell,
  levels = c(FALSE, TRUE),
  labels = c("Empty Droplet", "Cell")
)
p1 <- ggplot(plot_df, aes(x = log10(total_umi + 1), fill = is_cell)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("#A23B72", "#2E86AB")) +
  labs(
    title = "EmptyDrops: Cell vs Empty Droplet Detection",
    subtitle = paste(
      sum(plot_df$is_cell == "Cell"),
      "cells called from",
      ncol(raw_mtx),
      "total droplets"
    ),
    x = "log10(UMI + 1)",
    y = "Number of Droplets",
    fill = "Classification"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 14))

p1

emptydrops_summary <- function(raw_mtx, keep) {
  total <- ncol(raw_mtx)
  cells <- sum(keep)
  empty <- total - cells
  
  cat(
    "Droplets tested:", total, "\n",
    "Cells called:", cells, "(", round(cells/total*100, 1), "% )\n",
    "Empty droplets removed:", empty, "(", round(empty/total*100, 1), "% )\n",
    "After EmptyDrops:", cells, "cells retained\n"
  )
}

emptydrops_summary(raw_mtx, keep)

# Create plots directories
dir.create("plots", showWarnings = FALSE)
ggsave("plots/01_emptydrops_umi_distribution.png", p1,
       width = 10, height = 6, dpi = 300)


# write out filtered results with only validated cells
write10xCounts('counts_filtered', gene.symbol = genes[,2], filtered_mtx, overwrite=T) 


