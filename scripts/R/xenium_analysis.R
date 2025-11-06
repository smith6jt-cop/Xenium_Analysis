# Advanced Statistical Analysis with R
# This script demonstrates R integration for Xenium analysis
# Uses Seurat, ggplot2, and other R packages

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(reticulate)

# Set up Python integration
use_condaenv("xenium_analysis")

# Import Python modules
sc <- import("scanpy")
pd <- import("pandas")

#' Load Xenium Data from h5ad
#'
#' @param file_path Path to h5ad file
#' @return Seurat object
load_xenium_h5ad <- function(file_path) {
  message("Loading data from: ", file_path)
  
  # Read h5ad using scanpy
  adata <- sc$read_h5ad(file_path)
  
  # Convert to Seurat
  # Extract expression matrix
  counts <- t(as.matrix(adata$X))
  rownames(counts) <- adata$var_names$to_list()
  colnames(counts) <- adata$obs_names$to_list()
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = as.data.frame(adata$obs)
  )
  
  # Add spatial coordinates if available
  if ("spatial" %in% names(adata$obsm)) {
    spatial_coords <- as.matrix(adata$obsm[['spatial']])
    colnames(spatial_coords) <- c("x", "y")
    seurat_obj@meta.data <- cbind(seurat_obj@meta.data, spatial_coords)
  }
  
  message("Loaded Seurat object with ", ncol(seurat_obj), " cells and ", 
          nrow(seurat_obj), " genes")
  
  return(seurat_obj)
}

#' Perform Differential Expression with DESeq2-like approach
#'
#' @param seurat_obj Seurat object
#' @param group_var Variable for grouping
#' @param ident_1 First group
#' @param ident_2 Second group
#' @return Data frame with DE results
run_de_analysis <- function(seurat_obj, group_var, ident_1, ident_2) {
  message("Running DE analysis: ", ident_1, " vs ", ident_2)
  
  Idents(seurat_obj) <- group_var
  
  de_results <- FindMarkers(
    seurat_obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  
  de_results$gene <- rownames(de_results)
  de_results <- de_results %>%
    arrange(p_val_adj)
  
  return(de_results)
}

#' Create enhanced volcano plot
#'
#' @param de_results Data frame with DE results
#' @param title Plot title
#' @return ggplot object
plot_volcano <- function(de_results, title = "Volcano Plot") {
  de_results$significant <- ifelse(
    de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 0.5,
    "Significant",
    "Not Significant"
  )
  
  p <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                               color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
    theme_minimal() +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

#' Spatial feature plots
#'
#' @param seurat_obj Seurat object with spatial coordinates
#' @param features Features to plot
#' @return ggplot object
plot_spatial_features <- function(seurat_obj, features) {
  if (!all(c("x", "y") %in% colnames(seurat_obj@meta.data))) {
    stop("Spatial coordinates (x, y) not found in metadata")
  }
  
  plots <- lapply(features, function(feature) {
    if (feature %in% rownames(seurat_obj)) {
      # Gene expression
      expr <- FetchData(seurat_obj, vars = feature)
      data <- cbind(seurat_obj@meta.data[, c("x", "y")], expr)
      
      p <- ggplot(data, aes(x = x, y = y, color = .data[[feature]])) +
        geom_point(size = 0.5) +
        scale_color_viridis_c() +
        theme_void() +
        labs(title = feature) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right"
        )
    } else {
      # Metadata
      data <- seurat_obj@meta.data[, c("x", "y", feature)]
      
      p <- ggplot(data, aes(x = x, y = y, color = .data[[feature]])) +
        geom_point(size = 0.5) +
        theme_void() +
        labs(title = feature) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right"
        )
    }
    
    return(p)
  })
  
  # Combine plots
  combined <- wrap_plots(plots, ncol = 2)
  return(combined)
}

#' Cell type composition analysis
#'
#' @param seurat_obj Seurat object
#' @param celltype_col Cell type column
#' @param group_col Grouping column
#' @return ggplot object
plot_composition <- function(seurat_obj, celltype_col, group_col) {
  composition <- seurat_obj@meta.data %>%
    group_by(.data[[group_col]], .data[[celltype_col]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(.data[[group_col]]) %>%
    mutate(
      total = sum(count),
      percentage = count / total * 100
    )
  
  p <- ggplot(composition, aes(x = .data[[celltype_col]], y = percentage, 
                                fill = .data[[group_col]])) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      title = "Cell Type Composition",
      x = "Cell Type",
      y = "Percentage (%)",
      fill = group_col
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

#' Main analysis workflow
main <- function() {
  # Example usage
  message("=== Xenium R Analysis Workflow ===\n")
  
  # Define paths
  data_dir <- "../data/processed"
  output_dir <- "../figures/R_analysis"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load data
  # seurat_obj <- load_xenium_h5ad(file.path(data_dir, "xenium_sample_01_annotated.h5ad"))
  
  message("\nAnalysis functions loaded successfully!")
  message("Use the following functions for analysis:")
  message("  - load_xenium_h5ad(): Load h5ad files")
  message("  - run_de_analysis(): Differential expression")
  message("  - plot_volcano(): Volcano plots")
  message("  - plot_spatial_features(): Spatial visualization")
  message("  - plot_composition(): Cell type composition")
}

# Run main if script is executed
if (!interactive()) {
  main()
}
