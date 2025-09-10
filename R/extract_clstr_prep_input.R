#' Extract Cluster from Seurat Object and Return Pseudobulk and Variance Matrices
#'
#' This function processes a Seurat object to filter cells, select top genes, optionally normalize counts,
#' and compute pseudobulk data. It can generate pseudobulks per subject or create pseudoreplicates
#' from each subject's cells. Data is kept in a sparse matrix format for as long as possible to
#' improve performance and reduce memory usage.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param cluster_id The cluster ID or name to subset the Seurat object.
#' @param assay The assay to use (default is "RNA").
#' @param top_n_genes The number of top genes to retain based on mean expression (default is 10000).
#' @param cells_threshold The minimum number of cells required per individual (default is 10).
#' @param clip_sd The number of standard deviations to use for clipping cells based on total counts (default is 2).
#' @param normalize Logical indicating whether to normalize the data using edgeR's TMM method (default is TRUE).
#' @param num_replicates An integer specifying the number of pseudoreplicates to create per subject.
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#' @param random_seed An integer seed for reproducibility when creating pseudoreplicates.
#' @return A list containing pseudobulk normalized counts and variance matrices.
#'
#' @examples
#' # Example usage for pseudobulk per subject:
#' # pseudobulk_data <- extract_clstr_prep_input(seurat_obj, cluster_id = "2")
#'
#' # Example usage for 5 pseudoreplicates per subject:
#' # pseudorep_data <- extract_clstr_prep_input(seurat_obj, cluster_id = "2", num_replicates = 5)
#'
#' @importFrom Seurat GetAssayData
#' @importFrom SeuratObject Layers
#' @importFrom stats sd
#' @importFrom rlang is_empty
#' @importFrom Matrix rowMeans colSums
#' @importFrom matrixStats rowVars
#' @importFrom edgeR calcNormFactors DGEList
#' @export
extract_clstr_prep_input = function(seurat_obj, cluster_id, assay = "RNA", top_n_genes = 10000,
            cells_threshold = 10, clip_sd = NULL, normalize = TRUE, num_replicates = NULL,
            verbose = TRUE, random_seed = 42){
    # --- Input Validation for Pseudoreplicates ---
    if (!is.null(num_replicates) && num_replicates <= 0) {
        stop("Please specify 'num_replicates' > 0")
    }

    if(verbose) {
        message("Extracting cluster ", cluster_id, " from Seurat object...")
    }
    # Subset to the cluster of interest
    cluster <- subset(seurat_obj, subset = seurat_clusters == cluster_id)

    # Filter out cells with outlier total counts
    if(!is.null(clip_sd) && !is.na(clip_sd) && clip_sd > 0 && ncol(cluster) > 1) {
        log_counts <- log(cluster$nCount_RNA)
        max_cts_lim <- exp(mean(log_counts) + clip_sd * sd(log_counts))
        min_cts_lim <- exp(mean(log_counts) - clip_sd * sd(log_counts))
        cluster <- subset(cluster, subset = nCount_RNA > min_cts_lim & nCount_RNA < max_cts_lim)
        if(verbose) {
            message("Filtering cells with outlier total counts using ", clip_sd, " standard deviations...")
        }
    }
       
    counts <- cluster@assays[[assay]]$counts
    # Normalize data with edgeR's TMM method
    if(normalize){
        if(verbose) {
            message("Normalizing data using edgeR's TMM method...")
            message("Creating DGEList object and calculating normalization factors...")
        }
        y <- DGEList(counts = counts)
        y <- calcNormFactors(y, method = "TMM")
        edgeR_factors <- y$samples$norm.factors

        # calc normalized counts, keep the data on the raw count scale:

        #calculate the "biologically effective" library size for each cell 
        #by adjusting the raw library size (colSums(counts)) with the TMM composition factor. 
        effective_lib_size <- Matrix::colSums(counts) * edgeR_factors
        #calculate the average of these effective library sizes across all cells.
        #It represents the typical library size for a cell in this specific dataset.
        common_lib_size <- mean(effective_lib_size)
        #For each cell, calculate a final scaling factor. This factor represents how much
        # larger or smaller that cell's effective library is compared to the average.
        norm_factors_per_sample <- effective_lib_size / common_lib_size
        norm_factors_expanded <- rep.int(norm_factors_per_sample, diff(counts@p))
        count_data_out <- counts
        #divide every count value in a cell by that cell's specific norm_factors_per_sample.
        count_data_out@x <- counts@x / norm_factors_expanded

        #Above is in contrast to using edgeR's cpm(), which would return values on a much larger scale.
    } else {
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

    # Filter the final count matrix to match the individuals being kept
    # This ensures the returned counts object is consistent with the pseudobulk data.
    original_cell_count <- ncol(count_data_out)
    cells_to_keep_mask <- meta$orig.ident %in% individuals_to_keep
    
    meta <- meta[cells_to_keep_mask, , drop = FALSE]
    count_data_out <- count_data_out[, cells_to_keep_mask, drop = FALSE]

    if(verbose && (length(individuals_to_keep) < length(cell_counts))) {
        message("Removed ", original_cell_count - ncol(count_data_out), " cells belonging to individuals with fewer than ", cells_threshold, " cells.")
    }
    # --------------------------------------------------------------------------
    if (is_empty(individuals_to_keep)) {
        warning(paste("No individuals found with at least", cells_threshold, "cells for cluster", cluster_id))
        empty_df <- data.frame(matrix(ncol = 0, nrow = nrow(count_data_out)))
        rownames(empty_df) <- rownames(count_data_out)
        return(list(pseudobulk_tmm_matrix = empty_df, var_df = empty_df))
    }
    
    # Initialize lists to store results
    pseudobulk_list <- list()
    variance_list <- list()
    set.seed(random_seed)

    # Loop through each individual
    for (id in individuals_to_keep) {
        cell_indices <- which(meta$orig.ident == id)
        n_cells_subject <- length(cell_indices)

        # Handle pseudoreplicates or standard pseudobulk
        create_replicates <- !is.null(num_replicates)

        if (!create_replicates) {
            # Standard pseudobulk per subject
            individual_sparse_matrix <- count_data_out[, cell_indices, drop = FALSE]
            pseudobulk_list[[id]] <- Matrix::rowMeans(individual_sparse_matrix)
            variance_list[[id]] <- matrixStats::rowVars(as.matrix(individual_sparse_matrix))
        } else {
            # Pseudoreplicate generation
            n_reps <- num_replicates

            # Randomly assign cells to replicates
            shuffled_indices <- sample(cell_indices)
            rep_assignments <- cut(seq_along(shuffled_indices), breaks = n_reps, labels = FALSE)
            
            for (r in 1:n_reps) {
                replicate_cell_indices <- shuffled_indices[rep_assignments == r]
                replicate_name <- paste(id, r, sep = "_rep")
                
                replicate_sparse_matrix <- count_data_out[, replicate_cell_indices, drop = FALSE]
                
                pseudobulk_list[[replicate_name]] <- Matrix::rowMeans(replicate_sparse_matrix)
                variance_list[[replicate_name]] <- matrixStats::rowVars(as.matrix(replicate_sparse_matrix))
            }
        }
    }
    
    if (rlang::is_empty(pseudobulk_list)) {
        warning(paste("No pseudobulks or pseudoreplicates were generated for cluster", cluster_id))
        empty_df <- data.frame(matrix(ncol = 0, nrow = nrow(count_data_out)))
        rownames(empty_df) <- rownames(count_data_out)
        return(list(pseudobulk_tmm_matrix = empty_df, var_df = empty_df))
    }

    # Combine lists into final data frames
    pseudobulk_tmm_matrix <- as.data.frame(do.call(cbind, pseudobulk_list))
    var_df <- as.data.frame(do.call(cbind, variance_list))

    # Final check for consistency 
    stopifnot(all(rownames(pseudobulk_tmm_matrix) == rownames(var_df)))
    stopifnot(all(colnames(pseudobulk_tmm_matrix) == colnames(var_df)))
    
    return(list(pseudobulk_tmm_matrix = pseudobulk_tmm_matrix, 
                var_df = var_df, counts = count_data_out))
}