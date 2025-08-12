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
#' @param normalize Logical indicating whether to normalize the data using edgeR's TMM method (default is TRUE).
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#' @return A list containing pseudobulk normalized counts and variance matrices.
#'
#' @examples
#' # Example usage:
#' # pseudobulk_data <- extract_clstr_prep_input(seurat_obj, cluster_id = "2")
#' @importFrom Seurat NormalizeData GetAssayData
#' @importFrom SeuratObject Layers
#' @importFrom stats sd
#' @importFrom rlang is_empty
#' @importFrom Matrix rowMeans colSums
#' @importFrom matrixStats rowVars
#' @importFrom edgeR calcNormFactors DGEList
#' @export
extract_clstr_prep_input = function(seurat_obj, cluster_id, assay = "RNA", top_n_genes = 10000, cells_threshold = 10, clip_sd = 2, normalize = TRUE, verbose = TRUE){
    if(verbose) {
        message("Extracting cluster ", cluster_id, " from Seurat object...")
    }
    # Subset to the cluster of interest
    cluster <- subset(seurat_obj, subset = seurat_clusters == cluster_id)

    # Filter out cells with outlier total counts
    if (ncol(cluster) > 1) {
        log_counts <- log(cluster$nCount_RNA)
        max_cts_lim <- exp(mean(log_counts) + clip_sd * sd(log_counts))
        min_cts_lim <- exp(mean(log_counts) - clip_sd * sd(log_counts))
        cluster <- subset(cluster, subset = nCount_RNA > min_cts_lim & nCount_RNA < max_cts_lim)
        if(verbose) {
            message("Removing cells with counts < ", round(min_cts_lim, 0), " or > ", round(max_cts_lim, 0))
        }
    }
    
    # Normalize data with edgeR's TMM method
    counts <- cluster@assays[[assay]]$counts
    if(normalize){
        if(verbose) {
            message("Normalizing data using edgeR's TMM method...")
            message("Creating DGEList object and calculating normalization factors...")
        }
        y <- DGEList(counts = counts)
        y <- calcNormFactors(y, method = "TMM")
        edgeR_factors <- y$samples$norm.factors

        # calc normalized counts, keep the data on the raw count scale
        effective_lib_size <- Matrix::colSums(counts) * edgeR_factors
        common_lib_size <- mean(effective_lib_size)
        norm_factors_per_sample <- effective_lib_size / common_lib_size
        norm_factors_expanded <- rep.int(norm_factors_per_sample, diff(counts@p))
        count_data_out <- counts
        count_data_out@x <- counts@x / norm_factors_expanded
    }else{
        if(verbose) {
            message("Skipping normalization, using raw counts.")
        }
        count_data_out <- counts
    }


    # Keep top genes by mean expression across all cells
    top_n_genes <- min(top_n_genes, nrow(count_data_out))
    gene_means <- Matrix::rowMeans(count_data_out)
    gene_mean_idx <- order(gene_means, decreasing = TRUE)[1:top_n_genes]
    if(verbose) {
        message("Keeping top ", top_n_genes, " genes based on mean expression.")
    }
    count_data_out <- count_data_out[gene_mean_idx, ]

   
    # Get metadata and identify individuals to keep based on cell count threshold
    meta <- cluster@meta.data
    cell_counts <- table(meta$orig.ident)
    individuals_to_keep <- names(cell_counts[cell_counts >= cells_threshold])

    if (is_empty(individuals_to_keep)) {
        warning(paste("No individuals found with at least", cells_threshold, "cells for cluster", cluster_id))
        # Return empty list or empty data frames, consistent with the function's output format
        empty_df <- data.frame(matrix(ncol = 0, nrow = nrow(count_data_out)))
        rownames(empty_df) <- rownames(count_data_out)
        return(list(pseudobulk_tmm_matrix = empty_df, var_df = empty_df))
    }
    
    # Pre-allocate result matrices for efficiency
    n_genes <- nrow(count_data_out)
    n_individuals <- length(individuals_to_keep)
    gene_names <- rownames(count_data_out)
    
    pseudobulk_tmm_matrix <- matrix(NA, nrow = n_genes, ncol = n_individuals, 
                                   dimnames = list(gene_names, individuals_to_keep))
    var_df <- matrix(NA, nrow = n_genes, ncol = n_individuals, 
                     dimnames = list(gene_names, individuals_to_keep))

    # Loop through each individual, subset the sparse matrix, and calculate stats
    for (id in individuals_to_keep) {
        # Get column indices for cells belonging to the current individual
        cell_indices <- which(meta$orig.ident == id)
        
        # Subset the sparse matrix for the current individual. 
        # 'drop = FALSE' is crucial to maintain matrix structure
        individual_sparse_matrix <- count_data_out[, cell_indices, drop = FALSE]
        
        # Calculate mean per gene 
        pseudobulk_tmm_matrix[, id] <- Matrix::rowMeans(individual_sparse_matrix)
        
        # Calculate variance per gene 
        # For a single cell, variance is NA. rowVars handles this.
        var_df[, id] <- matrixStats::rowVars(as.matrix(individual_sparse_matrix))
    }
    
    # Final check for consistency 
    stopifnot(all(rownames(pseudobulk_tmm_matrix) == rownames(var_df)))
    stopifnot(all(colnames(pseudobulk_tmm_matrix) == colnames(var_df)))
    
    # Convert to data.frame 
    return(list(pseudobulk_tmm_matrix = as.data.frame(pseudobulk_tmm_matrix), 
                var_df = as.data.frame(var_df), counts = count_data_out))
}