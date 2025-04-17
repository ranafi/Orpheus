library(Seurat)
simulate_scRNAseq_circadian <- function(n = 10000, m = 5000, p_circ = 0.1,
                                        phi_sd = 0.3, A_range = c(.1, 5), G_range = c(0, 5),
                                        nb_size = 1) {
  set.seed(123)
  
  # Time points
  time_points <- seq(0, 2*pi, length.out = 10)
  cell_time_indices <- sort(sample(1:10, n, replace = TRUE))
  t_vals <- time_points[cell_time_indices]
  
  # Assign circadian vs non-circadian genes
  n_circ <- round(p_circ * m)
  stopifnot(n_circ <= m)
  circadian_genes <- 1:n_circ
  
  # Initialize expression matrix
  expr_matrix <- matrix(0, nrow = m, ncol = n)
  rownames(expr_matrix) <- paste0("Gene", 1:m)
  colnames(expr_matrix) <- paste0("Cell", 1:n)
  
  # Simulate circadian genes
  for (gene in circadian_genes) {
    A <- runif(1, A_range[1], A_range[2])
    G <- runif(1, G_range[1], G_range[2])
    G = max(A, G)
    phases_w_noise = t_vals + rnorm(n, 0, phi_sd)
    if (gene %% 2){
      circ_counts <- sapply(phases_w_noise, function(t){rpois(1, lambda = (A * cos(t) + G ))})
    }else{
      circ_counts <- sapply(phases_w_noise, function(t){rpois(1, lambda = (A * cos(t - pi) + G ))})
      
    }
    circ_counts[is.na(circ_counts)] <- 0  # Avoid negative values
    expr_matrix[gene, ] <- circ_counts
  }
  
  # Simulate non-circadian genes
  non_circadian_genes <- setdiff(1:m, circadian_genes)
  for (gene in non_circadian_genes) {
    G <- runif(1, G_range[1], G_range[2])
    # cell_noise = rnorm(n, 0, epsilon_sd)
    mu <- rep(G, n) #+ cell_noise
    my_counts = sapply(1:n, function(t){rpois(1, lambda = (mu))})
    expr_matrix[gene, ] <- my_counts
  }
  
  # Metadata
  cell_metadata <- data.frame(
    cell_id = colnames(expr_matrix),
    time_point = t_vals
  )
  
  return(list(counts = expr_matrix, cell_metadata = cell_metadata, circadian_genes = circadian_genes))
}



sim <- simulate_scRNAseq_circadian(n = 1000, m = 5000, p_circ = 1, phi_sd = 0.9,
                                   A_range = c(0, 5), G_range = c(0, 5))
counts <- sim$counts
metadata <- sim$cell_metadata
circadian_genes <- sim$circadian_genes

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
seurat_obj$orig.ident = seurat_obj$time_point
seurat_obj$seurat_clusters = "0"
###############ORPHEUS#################


tmp = extract_clstr_prep_input(seurat_obj, cluster_id = "0")

pb_mat = tmp[[1]]
var_mat = tmp[[2]]

ZT = as.numeric(colnames(pb_mat) )

pb_cycling = is_cycling(ZT, pb_mat, time_units = "rad")
var_cyclers = var_cycling(ZT, var_mat, pb_mat, time_units = "rad")


# JTK_res = var_cycling_JTK(ZT*12/pi, var_mat, pb_mat, out_file_circ = "input_12hr.txt", out_file_mean = "input_24hr.txt", out_file_pb = "input_pb.txt",out_dir = "output/")
# 
# doub_cyclers_JTK = JTK_res[[1]]
# fund_cycling_JTK = JTK_res[[2]]
# pb_cyclers_JTK = JTK_res[[3]] %>% mutate(acrophase_bulk = acrophase_bulk*pi/12)

# var_cycling_JTK_res = merge(doub_cyclers_JTK, fund_cycling_JTK , by = "Gene_Symbols", all = TRUE) %>% mutate(double_freq_acrophase =double_freq_acrophase*pi/12 ) %>% mutate(fund_acro = fund_acro *pi/12)

phi = Orpheus::calc_phi(ZT, pb_cycling_results = pb_cycling, var_cycling_results = var_cyclers, var_p_thr = .3, expected_acrodiff_thr = .5)
hist(phi$phi)

input = var_cyclers %>% 
  left_join(pb_cycling, by = "Gene_Symbols") %>%
  mutate(across(-Gene_Symbols, as.numeric)) %>%
  filter(p_statistic < 0.05 & double_freq_p_statistic < .3) %>%
  mutate(acrodiff = (double_freq_acrophase - acrophase_bulk) %% pi) %>%
  mutate(keep_acrodiff = acrodiff < (pi / 2 + .5) & acrodiff > (pi / 2 - .5)) %>%
  filter(keep_acrodiff == TRUE) %>%
  dplyr::select(Gene_Symbols, fund_amp,var_mesor,double_freq_amplitude , everything()) %>%
  column_to_rownames("Gene_Symbols")%>%
  filter( amplitude_bulk > 0.5) %>% select(c(1:3)) %>%rename(Af = fund_amp, Mv = var_mesor, Av = double_freq_amplitude)

rm = which(is.na(phi$phi))

guess = estimate_single_phi(input[-rm,])
estimate_CI(input[-rm,], guess[1,3])



########################

