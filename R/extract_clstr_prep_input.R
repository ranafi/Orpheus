#' Extract Cluster from Seurat Object and Return Pseudobulk and Variance Matrices
#'
#' This function processes a Seurat object to filter cells based on thresholds. extract normalized counts, compute pseudobulk data, 
#' and compute variance matrices. It returns a list containing the pseudobulk normalized counts and variance matrices.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param cluster_id The cluster ID or name to subset the Seurat object.
#' @param assay The assay to use (default is "RNA").
#' @param top_n_genes The number of top genes to retain based on mean expression (default is 10000).
#' @param cells_threshold The minimum number of cells required per individual (default is 10).
#'
#' @return A list containing pseudobulk normalized counts and variance matrices.
#'
#' @examples
#' # Example usage:
#' # pseudobulk_data <- extract_clstr_prep_input(seurat_obj, cluster_id = "2")
#' @importFrom Seurat NormalizeData GetAssayData
#' @importFrom SeuratObject Layers
#' @importFrom dplyr %>% group_by summarise across everything
#' @importFrom Seurat NormalizeData GetAssayData
#' @importFrom stats sd var
#' @importFrom rlang is_empty
#' @importFrom Matrix rowMeans
#' @export
extract_clstr_prep_input = function(seurat_obj, cluster_id, assay = "RNA", top_n_genes = 10000, cells_threshold = 10, clip_sd = 2){
    cluster = subset(seurat_obj, subset = seurat_clusters == cluster_id)

    #filter out cells with log(cts) >10, and < 8
    max_cts_lim = exp(mean(log(cluster$nCount_RNA)) + clip_sd*sd(log(cluster$nCount_RNA)))
    min_cts_lim = exp(mean(log(cluster$nCount_RNA)) - clip_sd*sd(log(cluster$nCount_RNA)))
    # Subset the Seurat object
    cluster <- subset(cluster, subset = nCount_RNA > min_cts_lim & nCount_RNA < max_cts_lim)
    print(paste0("Removing cells with counts < ", min_cts_lim, " or > ", max_cts_lim))

    cluster = NormalizeData(cluster, assay = assay, normalization.method = "RC", scale.factor = mean(cluster$nCount_RNA))
    
    # Get all data layers
    layers <- Layers(cluster[[assay]], search = "data")
    # Extract normalized counts for each layer
    norm_counts_list <- lapply(layers, function(layer) {
    GetAssayData(cluster[[assay]], layer = layer)
    })

    # Named list with normalized counts per subject/layer
    names(norm_counts_list) <- layers

    # Here I make pseudobulk by summing NORMED cell counts:
    # Extract the RC-normalized values from the "data" slot
    rc_data <- do.call(cbind, norm_counts_list)

    #keep top genes by mean across all cells
    top_n_genes = min(top_n_genes, nrow(rc_data))
    gene_means = Matrix::rowMeans(rc_data)
    gene_mean_idx <- order(gene_means, decreasing = TRUE)[1:top_n_genes]  # Get idx of top n values
    rc_data = rc_data[gene_mean_idx,]


    # Extract metadata to get individual IDs
    meta <- cluster@meta.data

    # Convert to a dataframe (Cells x Genes)
    rc_data_df <- as.data.frame(t(as.matrix(rc_data)))
    rc_data_df$individualID <- meta$orig.ident  # Add individualID metadata


    # Compute the mean per gene for each individual (instead of sum)
    pseudobulk_rc <- rc_data_df %>%
    group_by(individualID) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  # Mean instead of sum

    # Convert back to genes x individuals format
    pseudobulk_rc_matrix <- as.data.frame(t(pseudobulk_rc[,-1]))  # Drop individualID column
    colnames(pseudobulk_rc_matrix) <- pseudobulk_rc$individualID  # Rename columns


    # Compute variance per gene for each individual
    var_df_tmp <- rc_data_df %>%
    group_by(individualID) %>%
    summarise(across(everything(), \(x) var(x, na.rm = TRUE)))  # Compute variance

    # Convert back to genes x individuals format
    var_df <- as.data.frame(t(var_df_tmp[,-1]))  # Drop individualID column
    colnames(var_df) <- var_df_tmp$individualID  # Rename columns


    stopifnot(all(rownames(pseudobulk_rc_matrix) == rownames(var_df)))
    stopifnot(all(rownames(pseudobulk_rc_matrix) == rownames(pseudobulk_rc_matrix)))

    #Exclude subjects with fewer than cells_threshold cells:
    rm_me = which(table(cluster$orig.ident) < cells_threshold)
    if(!is_empty(rm_me)){
        pseudobulk_rc_matrix = pseudobulk_rc_matrix[,-rm_me]
        var_df = var_df[,-rm_me]
    }
    return(list(pseudobulk_rc_matrix = pseudobulk_rc_matrix, var_df = var_df))
}