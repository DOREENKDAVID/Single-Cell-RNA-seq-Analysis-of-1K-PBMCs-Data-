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

#-------------------------------------------------------

# Generate QC report ----

#--------------------------------------------------------

# this QC plot has two important features:
# 1. 'knee point' - is the point where the signed curvature is minimized. 
# This corresponds to a transition between a distinct subset of barcodes with large totals and the majority of barcodes with smaller totals
# 2. 'inflection point' - is the point on the curve where the first derivative is minimized. 
# This corresponds to the point past which cells cannot reliably be distinguished from background



# source the R script that contains the bc_rank_plot and print_HTML functions we'll use to produce a QC report
# this script comes from Sarah Ennis's github repo here:  https://github.com/Sarah145/scRNA_pre_process
source('./functions.R') 

# load filtered mtx
filtered_mtx <- readMM('counts_filtered/matrix.mtx') 

# load run info from JSON files produced by Kb
kb_stats <- c(fromJSON(file = 'inspect.json'), 
              fromJSON(file = 'run_info.json')) 
# -------------------------------
# 1. Determine the sequencing chemistry
# -------------------------------
# Extract 10X chemistry version from kb_stats$call
chemistry_version <- grep('10X(.*)', strsplit(kb_stats$call, '\\s')[[1]], value = TRUE)

# -------------------------------
# 2. Summarize sequencing/alignment statistics
# -------------------------------
# Create a clean table of sequencing metrics
sequencing_summary <- data.frame(
  stat = c(
    'Sequencing technology', 
    'Number of reads processed', 
    '% reads pseudoaligned', 
    '% reads on whitelist'
  ),
  value = prettyNum(
    c(
      chemistry_version, 
      kb_stats$n_processed, 
      kb_stats$p_pseudoaligned, 
      round(kb_stats$percentageReadsOnWhitelist, 2)
    ), 
    big.mark = ','
  )
)

# -------------------------------
# 3. Calculate cell-level metrics
# -------------------------------
# Percent of total counts that are in filtered cells
percent_counts_in_cells <- round((sum(filtered_mtx) / sum(raw_mtx)) * 100, 2)

# Median UMI counts per cell
median_counts_per_cell <- median(colSums(filtered_mtx))

# Median number of genes detected per cell
median_genes_per_cell <- median(apply(filtered_mtx, 2, function(x) sum(x >= 1)))

# Total number of genes detected across all cells
total_genes_detected <- sum(rowSums(filtered_mtx) >= 1)

# Combine cell-level metrics into a table
cell_summary <- data.frame(
  stat = c(
    'Estimated number of cells',
    '% counts in cells',
    'Median counts per cell',
    'Median genes per cell',
    'Total genes detected'
  ),
  value = prettyNum(
    c(
      ncol(filtered_mtx),
      percent_counts_in_cells,
      median_counts_per_cell,
      median_genes_per_cell,
      total_genes_detected
    ),
    big.mark = ','
  )
)

# -------------------------------
# 4. Generate barcode rank statistics
# -------------------------------
barcode_stats <- barcodeRanks(raw_mtx)

# -------------------------------
# 5. Load barcode lists
# -------------------------------
raw_barcodes <- read.csv(
  'counts_unfiltered/cellranger/barcodes.tsv', 
  header = FALSE, sep = '\t'
)[,1]

filtered_barcodes <- read.csv(
  'counts_filtered/barcodes.tsv', 
  header = FALSE, sep = '\t'
)[,1]

# -------------------------------
# 6. Create and save barcode rank plot
# -------------------------------
bc_rank_plot(
  stats = barcode_stats,
  raw_cells = raw_barcodes,
  filt_cells = filtered_barcodes,
  save = 'counts_filtered/barcode_rank.png'
)

# Save barcode rank statistics to text file
write.table(
  as.data.frame(barcode_stats),
  file = "counts_filtered/barcode_rank.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------
# 7. Normalize output path
# -------------------------------
output_dir <- normalizePath("counts_filtered", winslash = "/", mustWork = TRUE)

# -------------------------------
# 8. Output HTML summary of sequencing run
# -------------------------------
print_HTML(
  seq_stats = sequencing_summary,
  cell_stats = cell_summary,
  dir = output_dir,
  sample_id = NULL
)


#-----------------------------------------------------------

# Ambient RNA contamination testing with SoupX

#------------------------------------------------------------
#tod ->table of counts()
#table of droplets
library(Matrix)

# RAW
raw_mtx <- readMM("counts_unfiltered/cellranger/matrix.mtx")
raw_barcodes <- readLines("counts_unfiltered/cellranger/barcodes.tsv")
raw_genes <- readLines("counts_unfiltered/cellranger/genes.tsv")

colnames(raw_mtx) <- raw_barcodes
rownames(raw_mtx) <- sapply(strsplit(raw_genes, "\t"), `[`, 2)

# FILTERED
filtered_mtx <- readMM("counts_filtered/matrix.mtx")
filtered_barcodes <- readLines("counts_filtered/barcodes.tsv")
filtered_genes <- readLines("counts_filtered/genes.tsv")

colnames(filtered_mtx) <- filtered_barcodes
rownames(filtered_mtx) <- sapply(strsplit(filtered_genes, "\t"), `[`, 2)

tod <- raw_mtx          # raw counts (all droplets)
toc <- filtered_mtx     # filtered  counts (cells only)


# create soupX channel
sc <- SoupChannel(
  tod = tod,
  toc = toc,
  calcSoupProfile = TRUE
)

#create temp surat object for clustering
temp_seurat <- CreateSeuratObject(counts = filtered_mtx)

temp_seurat <- NormalizeData(temp_seurat, verbose = FALSE)
temp_seurat <- FindVariableFeatures(temp_seurat, nfeatures = 2000, verbose = FALSE)
temp_seurat <- ScaleData(temp_seurat, verbose = FALSE)
temp_seurat <- RunPCA(temp_seurat, npcs = 30, verbose = FALSE)
temp_seurat <- FindNeighbors(temp_seurat, dims = 1:30, verbose = FALSE)
temp_seurat <- FindClusters(temp_seurat, resolution = 0.8, verbose = FALSE)


#assign cluster to soupX
sc <- setClusters(
  sc,
  clusters = setNames(
    as.character(temp_seurat$seurat_clusters),
    colnames(temp_seurat)
  )
)
# Check barcode overlap
sum(colnames(filtered_mtx) %in% colnames(temp_seurat))
length(colnames(filtered_mtx))
length(colnames(temp_seurat))



#estimate ambient dna contamination

sc <- tryCatch({
  autoEstCont(sc, verbose = FALSE)
}, error = function(e) {
  cat("⚠ autoEstCont failed – sample likely very clean\n")
  sc$fit$rho <- NULL
  return(sc)
})

#report contamination

contamination_fraction <- sc$fit$rho

cat(
  "Ambient RNA contamination:",
  ifelse(
    is.null(contamination_fraction),
    "NULL (very clean sample)",
    paste0(round(contamination_fraction * 100, 2), "%")
  ),
  "\n"
)


