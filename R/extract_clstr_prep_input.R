#' Extract Cluster from Seurat Object and Return Pseudobulk and Variance Matrices 
#'
#' This function processes a Seurat object to filter cells based on thresholds, extract normalized counts, compute pseudobulk data, 
#' and compute variance matrices, while keeping data in a sparse matrix format for as long as possible to improve performance and reduce memory usage.
#' It returns a list containing the pseudobulk normalized counts and variance matrices.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param cluster_id The cluster ID or name to subset the Seurat object.
#' @param assay The assay to use (default is "RNA").
#' @param top_n_genes The number of top genes to retain based on mean expression (default is 10000).
#' @param cells_threshold The minimum number of cells required per individual (default is 10).
#' @param clip_sd The number of standard deviations to use for clipping cells based on counts (default is 2).
#' @return A list containing pseudobulk normalized counts and variance matrices.
#'
#' @examples
#' # Example usage:
#' # pseudobulk_data <- extract_clstr_prep_input(seurat_obj, cluster_id = "2")
#' @importFrom Seurat NormalizeData GetAssayData
#' @importFrom SeuratObject Layers
#' @importFrom stats sd
#' @importFrom rlang is_empty
#' @importFrom Matrix rowMeans
#' @importFrom matrixStats rowVars
#' @export
extract_clstr_prep_input = function(seurat_obj, cluster_id, assay = "RNA", top_n_genes = 10000, cells_threshold = 10, clip_sd = 2){
    # Subset to the cluster of interest
    cluster <- subset(seurat_obj, subset = seurat_clusters == cluster_id)

    # Filter out cells with outlier total counts
    if (ncol(cluster) > 1) {
        log_counts <- log(cluster$nCount_RNA)
        max_cts_lim <- exp(mean(log_counts) + clip_sd * sd(log_counts))
        min_cts_lim <- exp(mean(log_counts) - clip_sd * sd(log_counts))
        cluster <- subset(cluster, subset = nCount_RNA > min_cts_lim & nCount_RNA < max_cts_lim)
        print(paste0("Removing cells with counts < ", round(min_cts_lim, 2), " or > ", round(max_cts_lim, 2)))
    }
    
    # Normalize data
    cluster <- NormalizeData(cluster, assay = assay, normalization.method = "RC", scale.factor = mean(cluster$nCount_RNA))
    
    # Get all data layers and combine them into a single sparse matrix
    layers <- Layers(cluster[[assay]], search = "data")
    norm_counts_list <- lapply(layers, function(layer) {
        GetAssayData(cluster[[assay]], layer = layer)
    })
    rc_data <- do.call(cbind, norm_counts_list)

    # Keep top genes by mean expression across all cells
    top_n_genes <- min(top_n_genes, nrow(rc_data))
    gene_means <- Matrix::rowMeans(rc_data)
    gene_mean_idx <- order(gene_means, decreasing = TRUE)[1:top_n_genes]
    rc_data <- rc_data[gene_mean_idx, ]

   
    # Get metadata and identify individuals to keep based on cell count threshold
    meta <- cluster@meta.data
    cell_counts <- table(meta$orig.ident)
    individuals_to_keep <- names(cell_counts[cell_counts >= cells_threshold])

    if (is_empty(individuals_to_keep)) {
        warning(paste("No individuals found with at least", cells_threshold, "cells for cluster", cluster_id))
        # Return empty list or empty data frames, consistent with the function's output format
        empty_df <- data.frame(matrix(ncol = 0, nrow = nrow(rc_data)))
        rownames(empty_df) <- rownames(rc_data)
        return(list(pseudobulk_rc_matrix = empty_df, var_df = empty_df))
    }
    
    # Pre-allocate result matrices for efficiency
    n_genes <- nrow(rc_data)
    n_individuals <- length(individuals_to_keep)
    gene_names <- rownames(rc_data)
    
    pseudobulk_rc_matrix <- matrix(NA, nrow = n_genes, ncol = n_individuals, 
                                   dimnames = list(gene_names, individuals_to_keep))
    var_df <- matrix(NA, nrow = n_genes, ncol = n_individuals, 
                     dimnames = list(gene_names, individuals_to_keep))

    # Loop through each individual, subset the sparse matrix, and calculate stats
    for (id in individuals_to_keep) {
        # Get column indices for cells belonging to the current individual
        cell_indices <- which(meta$orig.ident == id)
        
        # Subset the sparse matrix for the current individual. 
        # 'drop = FALSE' is crucial to maintain matrix structure
        individual_sparse_matrix <- rc_data[, cell_indices, drop = FALSE]
        
        # Calculate mean per gene 
        pseudobulk_rc_matrix[, id] <- Matrix::rowMeans(individual_sparse_matrix)
        
        # Calculate variance per gene 
        # For a single cell, variance is NA. rowVars handles this.
        var_df[, id] <- matrixStats::rowVars(as.matrix(individual_sparse_matrix))
    }
    
    # Final check for consistency 
    stopifnot(all(rownames(pseudobulk_rc_matrix) == rownames(var_df)))
    stopifnot(all(colnames(pseudobulk_rc_matrix) == colnames(var_df)))
    
    # Convert to data.frame 
    return(list(pseudobulk_rc_matrix = as.data.frame(pseudobulk_rc_matrix), 
                var_df = as.data.frame(var_df)))
}