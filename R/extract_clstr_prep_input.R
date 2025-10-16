#' Prepare a Seurat Object for ORPHEUS Analysis
#'
#' This function performs global, one-time pre-processing on a Seurat object. It filters
#' genes based on expression across all cells and then normalizes the data using edgeR's
#' TMM method, also across all cells. This ensures all cells are on a common, comparable
#' scale before being split into groups for analysis. The normalized counts are stored in a
#' new assay.
#'
#' @param seurat_obj A Seurat object with raw counts in the specified assay.
#' @param assay The assay containing raw counts (default is "RNA").
#' @param new_assay_name The name for the new assay that will store the normalized counts (default is "ORPHEUS_norm").
#' @param min_prop The minimum proportion of cells expressing a gene to retain it, passed to edgeR::filterByExpr() (default is 0.01).
#' @param min_count The minimum count threshold for a gene to be considered expressed in a cell (default is 1).
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#' @return A Seurat object with a new assay containing globally filtered and normalized counts.
#'
#' @importFrom Seurat GetAssayData CreateAssayObject
#' @importFrom Matrix colSums
#' @importFrom edgeR DGEList filterByExpr calcNormFactors
#' @export
prepare_seurat_for_orpheus <- function(seurat_obj, assay = "RNA", new_assay_name = "ORPHEUS_norm",
                                     min_prop = 0.01, min_count = 1, verbose = TRUE) {

    if (verbose) message("--- Starting Global Pre-processing for ORPHEUS ---")

    # Global Gene Filtering
    if (verbose) message("Filtering genes based on expression across cells (min.prop = ", min_prop, ")...")
    
    counts <- GetAssayData(seurat_obj, assay = assay, slot = "counts")
    
    keep <- sparse_filter_by_expr(counts, min_prop = min_prop, min_count = min_count)
    
    # Subset the Seurat object to keep only these genes
    genes_to_keep <- rownames(counts)[keep]
    seurat_obj <- seurat_obj[genes_to_keep, ]
    
    if (verbose) message(length(genes_to_keep), " genes kept after global filtering.")

    # Global TMM Normalization
    if (verbose) message("Performing global TMM normalization across ALL cells...")
    
    # Get the filtered counts
    counts <- GetAssayData(seurat_obj, assay = assay, slot = "counts")
    
    y <- DGEList(counts = counts)
    y <- calcNormFactors(y, method = "TMM")
    
    effective_lib_size <- y$samples$lib.size * y$samples$norm.factors
    common_lib_size <- mean(effective_lib_size)
    norm_factors_per_sample <- effective_lib_size / common_lib_size
    
    # Use the sweep method to keep data on the raw count scale
    normalized_counts <- sweep(as.matrix(counts), 2, norm_factors_per_sample, FUN = "/")
    
    # Store Normalized Data in a New Assay
    if (verbose) message("Storing normalized data in new assay: '", new_assay_name, "'")
    
    # Convert back to sparse matrix for efficiency
    normalized_counts_sparse <- as(normalized_counts, "sparseMatrix")
    
    seurat_obj[[new_assay_name]] <- CreateAssayObject(counts = normalized_counts_sparse)
    
    if (verbose) message("--- Global Pre-processing Complete ---")
    
    return(seurat_obj)
}

#' Calculate Pseudobulk and Variance Matrices for ORPHEUS
#'
#' This function takes a pre-processed Seurat object (ideally from `prepare_seurat_for_orpheus`)
#' and calculates pseudobulk mean and variance matrices for specified groups (e.g., cell types).
#' It can also generate pseudoreplicates within each group.
#'
#' @param seurat_obj A pre-processed Seurat object with globally normalized counts.
#' @param groups_to_analyze A character vector of the group IDs to analyze (e.g., c("T cells", "B cells")).
#' @param grouping_var The column name in the Seurat metadata that defines the groups (e.g., "celltype").
#' @param assay_to_use The assay to use for calculations, which should contain the normalized counts
#'        (default is "ORPHEUS_norm").
#' @param cells_threshold The minimum number of cells required per individual within a group (default is 10).
#' @param num_replicates An integer specifying the number of pseudoreplicates to create per subject.
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#' @param random_seed An integer seed for reproducibility when creating pseudoreplicates.
#' @return A nested list where each element corresponds to a group and contains the
#'         pseudobulk matrix, variance matrix, and the filtered counts for that group.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom rlang is_empty
#' @importFrom Matrix rowMeans
#' @importFrom matrixStats rowVars
#' @export
calculate_orpheus_inputs <- function(seurat_obj, groups_to_analyze, grouping_var, 
                                     assay_to_use = "ORPHEUS_norm", cells_threshold = 10,
                                     num_replicates = NULL, verbose = TRUE, random_seed = 42) {
    
    if (verbose) message("--- Calculating ORPHEUS inputs for specified groups ---")
    
    # Initialize a list to hold all results
    results_list <- list()
    
    # Loop through each specified group (e.g., cell type)
    for (group_id in groups_to_analyze) {
        
        if (verbose) message("\nProcessing group: ", group_id)
        
        # logical vec of cells in the group of interest
        cells_to_keep_vector <- seurat_obj@meta.data[[grouping_var]] == group_id
        cluster_cells <- seurat_obj[, cells_to_keep_vector]
        
        # Get the already-normalized counts for this subset
        count_data_out <- GetAssayData(cluster_cells, assay = assay_to_use, slot = "counts")
        
        meta <- cluster_cells@meta.data
        cell_counts <- table(meta$orig.ident)
        individuals_to_keep <- names(cell_counts[cell_counts >= cells_threshold])
        
        if (is_empty(individuals_to_keep)) {
            warning(paste("For group '", group_id, "', no individuals found with at least ", cells_threshold, " cells. Skipping."))
            results_list[[group_id]] <- NULL # Explicitly mark as skipped
            next
        }
        
        # Filter metadata and counts to only the cells from individuals we are keeping
        cells_to_keep_mask <- meta$orig.ident %in% individuals_to_keep
        meta <- meta[cells_to_keep_mask, , drop = FALSE]
        count_data_out <- count_data_out[, cells_to_keep_mask, drop = FALSE]
        
        # Initialize lists for this group's results
        pseudobulk_list <- list()
        variance_list <- list()
        set.seed(random_seed)

        # Loop through each individual within the group
        for (id in individuals_to_keep) {
            cell_indices <- which(meta$orig.ident == id)
            
            create_replicates <- !is.null(num_replicates)
            
            if (!create_replicates) {
                individual_matrix <- as.matrix(count_data_out[, cell_indices, drop = FALSE])
                pseudobulk_list[[id]] <- rowMeans(individual_matrix)
                variance_list[[id]] <- matrixStats::rowVars(individual_matrix)
            } else {
                shuffled_indices <- sample(cell_indices)
                rep_assignments <- cut(seq_along(shuffled_indices), breaks = num_replicates, labels = FALSE)
                
                for (r in 1:num_replicates) {
                    replicate_cell_indices <- shuffled_indices[rep_assignments == r]
                    replicate_name <- paste(id, r, sep = "_rep")
                    
                    replicate_matrix <- as.matrix(count_data_out[, replicate_cell_indices, drop = FALSE])
                    
                    pseudobulk_list[[replicate_name]] <- rowMeans(replicate_matrix)
                    variance_list[[replicate_name]] <- matrixStats::rowVars(replicate_matrix)
                }
            }
        }
        
        if (is_empty(pseudobulk_list)) {
            warning(paste("No pseudobulks generated for group '", group_id, "'"))
            results_list[[group_id]] <- NULL
            next
        }
        
        # Combine lists into final data frames for this group
        pseudobulk_tmm_matrix <- as.data.frame(do.call(cbind, pseudobulk_list))
        var_df <- as.data.frame(do.call(cbind, variance_list))

        # Store results for this group
        results_list[[group_id]] <- list(
            pseudobulk_tmm_matrix = pseudobulk_tmm_matrix, 
            var_df = var_df, 
            counts = count_data_out
        )
    }
    
    if (verbose) message("\n--- All groups processed ---")
    return(results_list)
}

#' Lightweight Gene Filtering for Sparse Matrices
#'
#' This function replicates the core logic of edgeR::filterByExpr for ungrouped data,
#' but operates directly on sparse matrices to avoid a costly dense conversion.
#'
#' @param counts A sparse matrix (from the Matrix package) of counts.
#' @param min_prop The minimum proportion of cells that must express a gene.
#' @param min_count The minimum count required to be considered "expressed".
#' @return A logical vector indicating which genes (rows) to keep.
#'
#' @importFrom Matrix rowSums
#'
sparse_filter_by_expr <- function(counts, min_prop = 0.01, min_count = 10) {
    # Calculate the minimum number of cells that must express the gene
    min_n_cells <- ceiling(min_prop * ncol(counts))
    
    # The expression `counts > min_count` creates a sparse *logical* matrix,
    # and rowSums on that efficiently counts the TRUE values (i.e., cells > min_count).
    genes_to_keep <- Matrix::rowSums(counts > min_count) >= min_n_cells
    
    return(genes_to_keep)
}

utils::globalVariables(c(
  "as"
))