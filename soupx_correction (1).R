library(Seurat)
library(SoupX)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(R.utils)
library(SeuratData)
library(dplyr)
library(foreach)
library(doParallel)
library(foreach)
library(doParallel)



#' Run a standard clustering pipeline on single-cell data
#'
#' @param counts A matrix or data frame of raw counts for single-cell data.
#' @return A Seurat object with clustering information.
get_clusters <- function(counts){
  srat  <- CreateSeuratObject(counts = counts)
  srat    <- SCTransform(srat, verbose = F) # normalize
  srat    <- RunPCA(srat, verbose = F) # PCA
  srat    <- RunUMAP(srat, dims = 1:30, verbose = F) # UMAP
  srat    <- FindNeighbors(srat, dims = 1:30, verbose = F) # Neigbours
  srat    <- FindClusters(srat, verbose = T) # cluster cells
  return (srat)
}

#' Load STARsolo Velocyto Outputs into a Seurat Object
#'
#' This function reads and processes spliced, unspliced, and ambiguous matrices from a STARsolo output directory
#' to create a Seurat object, which can be used for single-cell RNA-seq analysis, including RNA velocity.
#'
#' @param data_dir A string specifying the directory containing STARsolo output files.
#'        The directory should contain `unspliced.mtx.gz`, `spliced.mtx.gz`, `ambiguous.mtx.gz`,
#'        `barcodes.tsv.gz`, and `features.tsv.gz` files.
#' @param return_multi_layers A logical value indicating whether to return the separate spliced, unspliced,
#'        and ambiguous layers within the Seurat object. Default is FALSE.
#'
#' @return A Seurat object containing the combined counts matrix as well as metadata columns with
#'         statistics for spliced, unspliced, and ambiguous reads. If `return_multi_layers` is TRUE,
#'         the function also includes separate layers for spliced, unspliced, and ambiguous counts.
#'
#' @details
#' - The function sums the spliced, unspliced, and ambiguous matrices to create a combined counts matrix.
#' - Barcodes and feature names are assigned as column and row names, respectively.
#' - Metadata columns are added to the Seurat object, containing the sum of spliced, unspliced, and ambiguous counts
#'   per cell, as well as the percentage of unspliced reads.
#' - Setting `return_multi_layers` to TRUE will store spliced, unspliced, and ambiguous matrices as
#'   separate layers within the Seurat object.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- read_STARsolo_velocyto("/path/to/STARsolo_output", return_multi_layers = TRUE)
#' }
#'
#' @export
read_STARsolo_velocyto <- function(data_dir, return_multi_layers = FALSE){
  unspliced <- readMM(file.path(data_dir, 'unspliced.mtx.gz'))
  spliced <- readMM(file.path(data_dir, 'spliced.mtx.gz'))
  ambiguous <- readMM(file.path(data_dir, 'ambiguous.mtx.gz'))

  # Summing matrices to create a combined counts matrix
  counts <- unspliced + spliced + ambiguous

  # Load cell barcodes and feature (gene) identifiers
  barcodes_d_path <- file.path(data_dir, "barcodes.tsv.gz")
  barcodes_d <- read.table(gzfile(barcodes_d_path), header = FALSE, stringsAsFactors = FALSE)
  features_d_path <- file.path(data_dir, "features.tsv.gz")
  features_d <- read.table(gzfile(features_d_path), header = FALSE, stringsAsFactors = FALSE)

  # Assign row and column names to the counts matrix
  rownames(counts) <- features_d$V1
  colnames(counts) <- barcodes_d$V1

  # Create a Seurat object using the combined counts matrix
  seurat_obj <- CreateSeuratObject(counts = counts)

  # Conditionally store spliced, unspliced, and ambiguous data in separate layers if return_multi_layers is TRUE
  if (return_multi_layers) {
    LayerData(seurat_obj, assay = "RNA", layer = "spliced") <- spliced
    LayerData(seurat_obj, assay = "RNA", layer = "unspliced") <- unspliced
    LayerData(seurat_obj, assay = "RNA", layer = "ambiguous") <- ambiguous
  }

  # Add metadata columns with counts and percentages for spliced, unspliced, and ambiguous data
  seurat_obj <- AddMetaData(seurat_obj, metadata = colSums(spliced), col.name = "sum_spliced")
  seurat_obj <- AddMetaData(seurat_obj, metadata = colSums(unspliced), col.name = "sum_unspliced")
  seurat_obj <- AddMetaData(seurat_obj, metadata = colSums(ambiguous), col.name = "sum_ambiguous")
  seurat_obj <- AddMetaData(seurat_obj, metadata = 100 * (colSums(unspliced) / (colSums(spliced) + colSums(unspliced))), col.name = "pct_unspliced")

  # Return the Seurat object with metadata for downstream analysis
  return(seurat_obj)
}

#' Process Single-Cell RNA-seq Data with RNA Velocity and SoupX Correction
#'
#' This function processes single-cell RNA-seq data for a given sample, including reading and combining
#' spliced, unspliced, and ambiguous matrices, performing clustering, UMAP embedding, and adjusting counts
#' using SoupX to correct for ambient RNA contamination. It also generates diagnostic plots and saves processed
#' outputs in a specified directory.
#'
#' @param sample A string specifying the sample name. This name is used to locate input data files and to
#'        define output directories and filenames.
#' @param cont_fraction A numeric value e.g. 0.2 representing the contamination fraction to use for SoupX
#'                      (default is NULL). If `NULL`, the contamination fraction will be estimated automatically.
#'
#' @return This function does not return any values but generates and saves the following outputs in a specified directory:
#' - Adjusted gene count matrices (after SoupX correction) in 10x format, compressed with gzip
#' - Various metadata files: `soupProfile.csv`, `fit_markersUsed.csv`, `metaData.csv`, and `obs_metadata.csv`
#' - A PDF with diagnostic plots for UMAP, barcode rank, and marker distributions
#'
#' @details
#' - This function first checks for existing processed output. If found, it skips processing for that sample.
#' - Velocyto matrices (spliced, unspliced, ambiguous) are read for both raw and filtered data.
#' - Clustering and dimensionality reduction are performed, and results are used in SoupX to identify and remove
#'   ambient RNA contamination.
#' - The corrected gene counts are saved in 10x format and compressed, along with various diagnostic and metadata files.
#'
#' @examples
#' \dontrun{
#' process_sample("sample_name")
#' }
#'
#' @export
process_sample <- function(sample, cont_fraction = NULL){
  # Print the sample name to track progress
  print(sample)

  # Define paths for sample input and output directories
  sample_data_dir = file.path(raw_data_location, tissue, sample)
  specific_export_velocyto_dir <- file.path(velocyto_export_dir, tissue, sample)
  if (!file.exists(specific_export_velocyto_dir)) { dir.create(specific_export_velocyto_dir, recursive = TRUE) }
  counts_dir <- file.path(specific_export_velocyto_dir, 'counts')

  # Skip processing if adjusted counts already exist for this sample
  if(file.exists(counts_dir)) { return() }

  # Load unspliced, spliced, and ambiguous matrices for both raw and filtered gene counts
  tod <- read_STARsolo_velocyto(data_dir = file.path(sample_data_dir, "Solo.out/Velocyto/raw/"))
  toc <- read_STARsolo_velocyto(data_dir = file.path(sample_data_dir, "Solo.out/Velocyto/filtered/"))

  # Perform clustering and UMAP embedding
  srat <- get_clusters(LayerData(toc, assay = "RNA", layer = "counts"))
  umap_pl <- DimPlot(srat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings

  # Set up SoupChannel for SoupX contamination correction
  sc <- SoupChannel(LayerData(tod, assay = "RNA", layer = "counts"), LayerData(toc, assay = "RNA", layer = "counts"), calcSoupProfile = FALSE)
  sc <- estimateSoup(sc)
  sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc <- setDR(sc, umap)

  # Estimate or manually set contamination fraction
  if (!is.null(cont_fraction)) {
    sc <- setContaminationFraction(sc, cont_fraction)
  } else {
    sc <- autoEstCont(sc)
  }

  # Open a PDF file to save diagnostic plots
  pdf(file.path(specific_export_velocyto_dir, "soupx_plots.pdf"))
  adj.matrix <- adjustCounts(sc, roundToInt = TRUE)

  # Save adjusted counts in 10x format and compress with gzip
  DropletUtils::write10xCounts(counts_dir, adj.matrix)
  for (file in list.files(counts_dir)) {
    R.utils::gzip(file.path(counts_dir, file))
  }

  # Copy Velocyto feature names file for consistency
  file.copy(from = file.path(sample_data_dir, "Solo.out/Velocyto/filtered/features.tsv.gz"),
            to = file.path(counts_dir, 'features.tsv.gz'))

  # Save metadata files for inspection and downstream analysis
  write.csv(sc$soupProfile, file.path(specific_export_velocyto_dir, 'soupProfile.csv'))
  write.csv(sc$fit$markersUsed, file.path(specific_export_velocyto_dir, 'fit_markersUsed.csv'))
  write.csv(sc$metaData, file.path(specific_export_velocyto_dir, 'metaData.csv'))
  write.csv(toc@meta.data, file.path(specific_export_velocyto_dir, 'obs_metadata.csv'))

  # Calculate genes per droplet and prepare data for plotting
  total_genes_per_droplet <- colSums(LayerData(tod, assay = "RNA", layer = "counts") > 0)

  # Create a rank plot for barcode quality assessment
  df <- data.frame(
    total_genes = sort(total_genes_per_droplet, decreasing = TRUE),
    index = seq_along(total_genes_per_droplet)
  )

  # Create the rank plot
  rankplot <- ggplot(df, aes(x = index, y = total_genes)) +
    geom_line() +
    geom_vline(xintercept = dim(toc)[2], color = 'red') +
    scale_y_log10() +
    scale_x_log10() +
    labs(
      title = paste0("Raw Barcode Rank Plot for ", sample, " - Gene"),
      x = "Droplets",
      y = "Number of UMIs"
    ) +
    theme_classic()

  # Plot and save UMAP, rank plot, and marker distribution plot in the PDF
  plot(umap_pl)
  plot(rankplot)
  try(plotMarkerDistribution(sc))

  # Close the PDF
  dev.off()
}

gzip <- R.utils::gzip




##### Example to run SoupX correction in parallel

# Set up parallel processing
num_cores <- parallel::detectCores() - 2  # Use all but two cores for parallel execution
cl <- parallel::makeCluster(num_cores) # Initialize a cluster with the specified number of cores
doParallel::registerDoParallel(cl) # Register the cluster for parallel processing

# Define file paths for raw data and export directories
raw_data_location <- file.path('/Users/emilioskarwan/Documents/data/macque/macaque_sc/raw_data')

# Define export directory paths for SoupX and Velocyto outputs
export_dir <- file.path('/Users/emilioskarwan/Documents/data/macque/macaque_sc/soupx_out_velocombo')
velocyto_export_dir <- file.path(export_dir, 'Velocyto')

# Create export directories if they do not exist
if (!file.exists(export_dir)) { dir.create(export_dir, recursive = TRUE)}
if (!file.exists(velocyto_export_dir)) { dir.create(velocyto_export_dir, recursive = TRUE)}

# List all tissue directories within the raw data directory
tissues <- list.files(raw_data_location)
tissues <- tissues
tissue <- tissues[[1]]

# Iterate over each tissue directory in parallel
foreach(tissue = tissues, .packages = c("Seurat", "SoupX", "DropletUtils", "Matrix", "ggplot2", "foreach", "R.utils")) %dopar% {
  # List all sample subdirectories for the current tissue
  samples <- list.files(file.path(raw_data_location, tissue))

    # Parallelize processing of each sample within the current tissue
  foreach(sample = samples) %dopar% {
    process_sample(sample) # Process each sample using the `process_sample` function
  }
}
# Stop the parallel cluster once all tasks are completed
parallel::stopCluster(cl)





##### Example using manual estimation.
gzip <- R.utils::gzip
# Set the number of cores to use
num_cores <- parallel::detectCores() - 6  # Use one less than the number of available cores
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

raw_data_location <- file.path('/Users/emilioskarwan/Documents/data/macque/macaque_sc/raw_data')

export_dir <- file.path('/Users/emilioskarwan/Documents/data/macque/macaque_sc/soupx_out_velocombo')
velocyto_export_dir <- file.path(export_dir, 'Velocyto_manual')

if (!file.exists(export_dir)) { dir.create(export_dir, recursive = TRUE)}
if (!file.exists(velocyto_export_dir)) { dir.create(velocyto_export_dir, recursive = TRUE)}


tissues <- list.files(raw_data_location)
tissues <- c("PBMC")
tissue <- "PBMC"

samples <- list.files(file.path(raw_data_location, tissue))
for (sample in samples){
  process_sample(sample, cont_fraction = 0.2)
}
#
#tissue <- 'Lymph_node'
#sample <- 'M3-lym-2_21'
#process_sample(sample, cont_fraction = 0.2)
#
#
#tissue <- 'Prostate'
#samples <- list.files(file.path(raw_data_location, tissue))
#for (sample in samples){
#  process_sample(sample, cont_fraction = 0.2)
#}
#
#
#tissue <- 'Epididymis'
#samples <- list.files(file.path(raw_data_location, tissue))
#for (sample in samples){
#  process_sample(sample, cont_fraction = 0.2)
#}








