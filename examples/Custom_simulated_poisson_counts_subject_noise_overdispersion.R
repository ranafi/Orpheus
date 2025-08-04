library(Seurat)
library(dplyr)
# setwd("~/Box Sync/Henry_stuff/Orpheus/Orpheus/R")
devtools::check()
library(Orpheus)

simulate_scRNAseq_circadian <- function(n = 10000, m = 5000, p_circ = 0.1,n_time_points = 14, 
                                        phi_sd = 0.3, A_range = c(.1, 5), G_range = c(0.1, 5), epsilon_sd = 0.,
                                        theta_range = c(1,1)) { 
  set.seed(123)
  
  # Time points
  time_points <- seq(0, 2*pi, length.out = n_time_points) # Discrete time points for sampling
  cell_time_indices <- sort(sample(1:n_time_points, n, replace = TRUE))
  t_vals <- time_points[cell_time_indices] # Actual time for each cell
  
  # Assign circadian vs non-circadian genes
  n_circ <- round(p_circ * m)
  stopifnot(n_circ <= m)
  circadian_gene_indices <- 1:n_circ # Row indices for circadian genes
  
  # Initialize expression matrix
  expr_matrix <- matrix(0, nrow = m, ncol = n)
  gene_names <- paste0("Gene", 1:m)
  rownames(expr_matrix) <- gene_names
  colnames(expr_matrix) <- paste0("Cell", 1:n)
  
  # Initialize data frame to store gene parameters
  gene_parameters <- data.frame(
    gene_id = gene_names,
    theta = 0,
    A = 0, # Amplitude, will be 0 for non-circadian
    G = 0, # Baseline expression
    V_w_theoret = 0,
    V_2w_theoret = 0,
    V_DC_theoret = 0,
    S_amp_theoret = 0,
    is_circadian = FALSE,
    stringsAsFactors = FALSE
  )

  rownames(gene_parameters) <- gene_names
  
  
  # Simulate circadian genes
  if (n_circ > 0) {
    for (i in 1:n_circ) {
      gene_idx <- circadian_gene_indices[i] # Actual row index in expr_matrix
      
      A_val <- runif(1, A_range[1], A_range[2])
      G_val <- runif(1, G_range[1], G_range[2])
      G_val <- max(A_val, G_val) # Ensure baseline is at least amplitude
      theta <- runif(1, theta_range[1], theta_range[2])
      
      esp = exp(-phi_sd^2)
      es2p = exp(-2*phi_sd^2)
      esp_half = exp(-phi_sd^2/2)
      V_omega = A_val * esp_half * (1 + 2*theta*G_val)
      V_2omega = (A_val^2/2) * (esp - (1+theta)*es2p)
      V_DC = G_val + theta*G_val^2 + (1+theta)*epsilon_sd^2 + (A_val^2/2)*(1 + theta - esp)
      
      
      # Store parameters
      gene_parameters[gene_idx, "A"] <- A_val
      gene_parameters[gene_idx, "G"] <- G_val
      gene_parameters[gene_idx, "theta"] <- theta
      gene_parameters[gene_idx, "is_circadian"] <- TRUE
      gene_parameters[gene_idx, "V_w_theoret"] <- V_omega
      gene_parameters[gene_idx, "V_2w_theoret"] <- V_2omega
      gene_parameters[gene_idx, "V_DC_theoret"] <- V_DC
      gene_parameters[gene_idx, "S_amp_theoret"] <- A_val*exp(-phi_sd^2/2)
      
      
      phases_w_noise <- t_vals + rnorm(n, 0, phi_sd) # Cell-specific phase noise
      cell_specific_epsilon_noise <- rnorm(n, 0, epsilon_sd) # Cell-specific expression noise
      
      if (gene_idx %% 2) { # Odd genes peak at t=0 (cosine)
        lambdas <- A_val * cos(phases_w_noise) + G_val + cell_specific_epsilon_noise
      } else { # Even genes peak at t=pi (cosine shifted)
        lambdas <- A_val * cos(phases_w_noise - pi) + G_val + cell_specific_epsilon_noise
      }
      
      # Ensure lambdas are non-negative (rpois requires lambda >= 0)
      non_negative_lambdas <- pmax(0, lambdas) 
      
     
      # For NB, mu is the mean. Size controls dispersion.
      # var = mu + mu^2/size. If size=1, var = mu + mu^2. For large mu, var ~ mu^2.
      counts <- rnbinom(n, size = 1/theta, mu = non_negative_lambdas)

      
      expr_matrix[gene_idx, ] <- counts
    }
  }
  
  # Simulate non-circadian genes
  non_circadian_gene_indices <- setdiff(1:m, circadian_gene_indices)
  if (length(non_circadian_gene_indices) > 0) {
    for (gene_idx in non_circadian_gene_indices) {
      G_val <- runif(1, G_range[1], G_range[2])
      theta <- runif(1, theta_range[1], theta_range[2])
      
      # Store parameters (A is already 0, is_circadian is already FALSE)
      gene_parameters[gene_idx, "G"] <- G_val
      gene_parameters[gene_idx, "theta"] <- theta
      
      
      cell_specific_epsilon_noise <- rnorm(n, 0, epsilon_sd)
      mus <- G_val + cell_specific_epsilon_noise
      
      non_negative_mus <- pmax(0, mus)
      counts <- rnbinom(n, size = 1/theta, mu = non_negative_mus)

      
      expr_matrix[gene_idx, ] <- counts
    }
  }
  
  # Metadata for cells
  cell_metadata <- data.frame(
    cell_id = colnames(expr_matrix),
    time_point = t_vals # The "true" sampled time for each cell
  )
  
  # Identify which gene names are circadian (can be derived from gene_parameters too)
  circadian_gene_names <- gene_parameters$gene_id[gene_parameters$is_circadian]
  
  return(list(
    counts = expr_matrix, 
    cell_metadata = cell_metadata, 
    gene_parameters = gene_parameters,
    circadian_gene_ids = circadian_gene_names # More informative than just indices
  ))
}

phi_sd = 0.0
epsilon_sd = 0.4
# sim <- simulate_scRNAseq_circadian(n = 1000, m = 5000, p_circ = .1, phi_sd = phi_sd,
#                                    A_range = c(0.1, 2), G_range = c(0.1, 2))
sim <- simulate_scRNAseq_circadian(n = 1000, m = 5000, p_circ = .1, phi_sd = phi_sd, n_time_points = 6,
                                   A_range = c(.1, 2), G_range = c(.1, 2), epsilon_sd = epsilon_sd, theta_range = c(0.1, 1))
counts <- sim$counts
metadata <- sim$cell_metadata
gene_info <- sim$gene_parameters

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
seurat_obj$orig.ident = seurat_obj$time_point
seurat_obj$seurat_clusters = "0"

vars = sparseMatrixStats::rowVars(seurat_obj@assays$RNA$counts)
means = sparseMatrixStats::rowMeans2(seurat_obj@assays$RNA$counts)
plot(log(means),log(vars))
lines(c(0,100), c(0, 200), col = "red")
lines(c(0,100), c(0,100), col = "blue")
# FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
###############ORPHEUS#################


tmp = extract_clstr_prep_input(seurat_obj, cluster_id = "0", clip_sd = 10)

pb_mat = tmp[[1]]
var_mat = tmp[[2]]

ZT = as.numeric(colnames(pb_mat) )

pb_cycling = is_cycling(ZT, pb_mat, time_units = "rad")
circ_genes = pb_cycling %>% filter(as.numeric(p_statistic) < 0.05) %>% select(Gene_Symbols) %>% unname %>% unlist
circ_genes_phases = pb_cycling %>% filter(as.numeric(p_statistic) < 0.05) %>% select(acrophase_bulk)  %>% unname %>% unlist %>% as.numeric

var_cyclers = var_cycling(ZT,var_mat[circ_genes,], pb_mat[circ_genes,], phase_pb = circ_genes_phases, time_units = "rad")
var_cyclers2 = var_cycling_fix_acro(ZT, var_mat[circ_genes,], pb_mat[circ_genes,], phase_pb = circ_genes_phases, time_units = "rad")

test = merge(var_cyclers, var_cyclers2, by = "Gene_Symbols")
plot(as.numeric(test$double_freq_amplitude.x), as.numeric(test$double_freq_amplitude.y))
lines(c(-100, 100), c(-100, 100))

# JTK_res = var_cycling_JTK(ZT*12/pi, var_mat, pb_mat, out_file_circ = "input_12hr.txt", out_file_mean = "input_24hr.txt", out_file_pb = "input_pb.txt",out_dir = "output/")
# 
# doub_cyclers_JTK = JTK_res[[1]]
# fund_cycling_JTK = JTK_res[[2]]
# pb_cyclers_JTK = JTK_res[[3]] %>% mutate(acrophase_bulk = acrophase_bulk*pi/12)

# var_cycling_JTK_res = merge(doub_cyclers_JTK, fund_cycling_JTK , by = "Gene_Symbols", all = TRUE) %>% mutate(double_freq_acrophase =double_freq_acrophase*pi/12 ) %>% mutate(fund_acro = fund_acro *pi/12)

# phi = Orpheus::calc_phi(ZT, pb_cycling_results = pb_cycling, var_cycling_results = var_cyclers, var_p_thr = .3, expected_acrodiff_thr = .5)
# hist(phi$phi)

# expected_acrodiff_thr = 3*pi/12

input = var_cyclers %>% 
  left_join(pb_cycling, by = "Gene_Symbols", suffix= c("_var", "_pb")) %>%
  mutate(across(-Gene_Symbols, as.numeric)) %>%
  filter(p_statistic < 0.05 & double_freq_p_statistic < 1) %>%
  mutate(acrodiff = (double_freq_acrophase - acrophase_bulk) %% pi) %>%
  # mutate(keep_acrodiff = acrodiff < (pi / 2 + expected_acrodiff_thr) & acrodiff > (pi / 2 - expected_acrodiff_thr)) %>%
  # filter(keep_acrodiff == TRUE) %>%
  dplyr::select(Gene_Symbols, double_freq_amplitude,fund_amp, amplitude_bulk,mesor_bulk, everything()) %>%
  column_to_rownames("Gene_Symbols")%>%
  filter( amplitude_bulk >= 0.5) %>%
  # select(c(1:5)) %>%
  filter(partial_model_AIC_var > full_model_AIC_var)

find_analytic = function(V_2w, V_w, S_amp, S_DC){
  
  # Solve for G: trivial
  G = S_DC
  
  #solve for Theta
  theta = ((V_w / S_amp) - 1) / (2 * G)

  R = V_2w / S_amp^2
  log_argument = (1 - 2*R) / (1 + theta)
  
  phi_squared = -log(log_argument)
  phi = sqrt(phi_squared)
  
  # Solve for Ac 
  Ac = S_amp / exp(-phi_squared / 2) # More direct than sqrt(exp(-phi^2))
  
  return(c(G=G, Ac=Ac, theta=theta, phi=phi))
}

phi = input %>%
  apply(., 1, function(row) {
      find_analytic(row[1], row[2], row[3], row[4])
  }) %>%
  t %>%
  as.data.frame
colnames(phi) = c("G_pred", "Ac_pred", "theta_pred","phi_pred")

# phi2 = calc_phi( pb_cycling, var_cyclers )

plot_df = merge(input, phi, by = 0) %>%
  merge(., gene_info, by.x="Row.names", by.y = 0)%>%
  mutate(phase_error = (phi_pred -phi_sd ))

theta = 1
# Prepare data for the line
x_line_values <- seq(0, max(plot_df$amplitude_bulk), by = 0.1)
# y_line_values <- -x_line_values^2 * (exp(1)^9 - 1 - theta) / (2*exp(1)^9)
# y_lines_values2 <- -(theta/2)*x_line_values^2
# y_line_values = -.5*exp(-2^2)*x_line_values^2*(exp(2^2)-theta - 1)
y_line_lb <- -(theta/2)*x_line_values^2
y_line_ub <- .5*x_line_values^2
line_data <- data.frame(x = x_line_values, y1 = y_line_lb, y2 = y_line_ub)

# Create the ggplot
library(ggplot2)
ggplot() +
  geom_point(
    data = plot_df,
    aes(x = amplitude_bulk, y = double_freq_amplitude,
        shape = is_circadian, color = phase_error)
  ) +
  geom_line(
    data = line_data,
    aes(x = x, y = y1),
    color = "blue", 
    linewidth = 1   
  ) +
  geom_line(
    data = line_data,
    aes(x = x, y = y2),
    color = "blue", 
    linewidth = 1   
  ) +
  coord_cartesian() +
  labs(
    x = "Amplitude bulk",
    y = "Double Frequency Amplitude (input$double_freq_amplitude)"
  ) +
  theme_minimal()+
  scale_colour_gradient(low = "blue", high = "red", na.value = "black", limits = c(-.1, 1))

