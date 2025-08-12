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
    theta_true = 0,
    A_true = 0, # Amplitude, will be 0 for non-circadian
    G_true = 0, # Baseline expression
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
      gene_parameters[gene_idx, "A_true"] <- A_val
      gene_parameters[gene_idx, "G_true"] <- G_val
      gene_parameters[gene_idx, "theta_true"] <- theta
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
      gene_parameters[gene_idx, "G_true"] <- G_val
      gene_parameters[gene_idx, "theta_true"] <- theta
      
      
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

phi_sd = 0.3
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

phi2 = calc_phi( pb_cycling, var_cyclers )
sum(phi -  phi2, na.rm = T)

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

plot(plot_df$phase_error, plot_df$amplitude_bulk)

Orpheus::plot_rhythms(ZT, var_mat, pb_mat, plot_me = c("Gene163", "Gene220"), time_units = "rads")

############### Repeat above but across many choices of phi ####################
output = list()
for(i in seq(0, 1.3, by = 0.1)){
  print(paste0("phi = ", i))
  phi_sd = i
  epsilon_sd = 0.3
  # sim <- simulate_scRNAseq_circadian(n = 5000, m = 10000, p_circ = .2, phi_sd = phi_sd, n_time_points = 12,
                                     # A_range = c(.1, 2), G_range = c(.1, 2), epsilon_sd = epsilon_sd, theta_range = c(0.1, 1))
  sim <- simulate_scRNAseq_circadian(n = 5000, m = 10000, p_circ = .2, phi_sd = phi_sd, n_time_points = 12,
                                     A_range = c(.1, 5), G_range = c(.1, 5), epsilon_sd = epsilon_sd, theta_range = c(0.1, 1))
  counts <- sim$counts
  metadata <- sim$cell_metadata
  gene_info <- sim$gene_parameters
  
  seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
  seurat_obj$orig.ident = seurat_obj$time_point
  seurat_obj$seurat_clusters = "0"
  
  tmp = extract_clstr_prep_input(seurat_obj, cluster_id = "0", clip_sd = 100, normalize = T)
  
  pb_mat = tmp[[1]]
  var_mat = tmp[[2]]
  counts = tmp[[3]]
  
  #Variance of NB is var = mean + mean^2*theta
  #rearrange to solve for theta = (var-mean)/mean^2
  time_wise_theta_meas <- as.data.frame(mapply(function(v,m){(v-m)/m^2}, var_mat, pb_mat))
  rownames(time_wise_theta_meas) = rownames(pb_mat)
  colnames(time_wise_theta_meas) = paste0(colnames(time_wise_theta_meas), "_theta_meas")
  vars = sparseMatrixStats::rowVars(counts)
  means = sparseMatrixStats::rowMeans2(counts)
  agg_theta_meas = data.frame(coarse_theta = (vars-means)/means^2)
  rownames(agg_theta_meas) = rownames(pb_mat)
  
  
  ZT = as.numeric(colnames(pb_mat) )
  
  pb_cycling = is_cycling(ZT, pb_mat, time_units = "rad")
  circ_genes = pb_cycling %>% filter(as.numeric(p_statistic) < 0.05) %>% select(Gene_Symbols) %>% unname %>% unlist
  circ_genes_phases = pb_cycling %>% filter(as.numeric(p_statistic) < 0.05) %>% select(acrophase_bulk)  %>% unname %>% unlist %>% as.numeric
  
  var_cyclers = var_cycling(ZT, var_mat[circ_genes,], pb_mat[circ_genes,], phase_pb = circ_genes_phases, time_units = "rad")
  phi = calc_phi( pb_cycling, var_cyclers, bulk_amp_thr = 0 )
  input = var_cyclers %>% 
    left_join(pb_cycling, by = "Gene_Symbols", suffix= c("_var", "_pb")) %>%
    mutate(across(-Gene_Symbols, as.numeric)) %>%
    filter(p_statistic < 0.05 & double_freq_p_statistic < 1) %>%
    mutate(acrodiff = (double_freq_acrophase - acrophase_bulk) %% pi) %>%
    # mutate(keep_acrodiff = acrodiff < (pi / 2 + expected_acrodiff_thr) & acrodiff > (pi / 2 - expected_acrodiff_thr)) %>%
    # filter(keep_acrodiff == TRUE) %>%
    dplyr::select(Gene_Symbols, double_freq_amplitude,fund_amp, amplitude_bulk,mesor_bulk, everything()) %>%
    column_to_rownames("Gene_Symbols")%>%
    # filter( amplitude_bulk >= 0.5) %>%
    # select(c(1:5)) %>%
    filter(partial_model_AIC_var > full_model_AIC_var)
  
  plot_df = merge(input, phi, by = 0) %>%
    merge(., gene_info, by.x="Row.names", by.y = 0)%>%
    mutate(phase_error = (phi -phi_sd ))
  plot_df$phi_true = phi_sd
  out = merge(plot_df, time_wise_theta_meas, by.x = "Row.names", by.y = 0)
  out = merge(out, agg_theta_meas, by.x = "Row.names", by.y = 0)
  
  output <- append(output, list(out))
}

#convert list into dataframe
df = do.call(rbind, output)
#Keep rows of results where bulk amp > thresh & not NaN
df_filt = dplyr::filter(df, !is.na(phi) & amplitude_bulk >= 0.5)

#Plot: Phi sd pred acc vs true phi sd
df_filt$phi_true_factor = as.factor(df_filt$phi_true*12/pi)
phi_acc_plot = ggplot(df_filt, aes(x = phi_true*12/pi,y =  phi*12/pi))+
  geom_boxplot(aes(group = phi_true_factor))+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 0.3,\ntheta_range=c(0.1,1)")+
  theme_minimal()+
  xlab("True phase std. (hrs)")+
  ylab("Prediceted phase std. (hrs)")+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/Simulated_phi_acc.pdf", plot = phi_acc_plot, width = 5, height = 5, units = "in")  

#Plot below uses "df" (which contains genes with bulk_amp < 0.5) (and phi == NA,
#which will be not be plotted because phase_error will be NA):

# model <- lm(abs(phase_error) ~ S_amp_theoret, data = df)
# summary_model <- summary(model)
# slope <- coef(model)[["S_amp_theoret"]]
# p_value <- summary_model$coefficients["S_amp_theoret", "Pr(>|t|)"]
# r_squared <- summary_model$r.squared
# label_text <- sprintf("Slope = %.3f\np = %.3g\nR² = %.3f", slope, p_value, r_squared)
bulk_amp_vs_phase_err = ggplot(df, aes(S_amp_theoret, abs(phase_error*12/pi)))+
  geom_point(col = "grey32")+ 
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal()+
  # annotate("text", x = Inf, y = -Inf, label = label_text, 
  #          hjust = 1.1, vjust = -0.5, size = 4, color = "black") +
  # ylim(0, 2)+
  ylab("Absolute error of phase std. (hrs)")+
  xlab("Simulated bulk amplitude")+
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 0.3,\ntheta_range=c(0.1,1)")+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/bulk_amp_vs_phase_err.pdf", plot = bulk_amp_vs_phase_err, width = 5, height = 5, units = "in")  

bulk_amp_vs_phase_err2 = ggplot(df, aes(amplitude_bulk, abs(phase_error*12/pi)))+
  geom_point(col = "grey32")+ 
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal()+
  # annotate("text", x = Inf, y = -Inf, label = label_text, 
  #          hjust = 1.1, vjust = -0.5, size = 4, color = "black") +
  # ylim(0, 2)+
  ylab("Absolute error of phase std. (hrs)")+
  xlab("Simulated bulk amplitude")+
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 0.3,\ntheta_range=c(0.1,1)")+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/bulk_amp_measured_vs_phase_err.pdf", plot = bulk_amp_vs_phase_err2, width = 5, height = 5, units = "in")  


theta_pred_vs_true = ggplot(df_filt, aes(theta_true, theta))+
  geom_point(col = "grey32")+ 
  theme_minimal()+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  xlab("Simulated overdispersion param.")+
  ylab("Predicted overdispersion param.")+
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 0.3,\ntheta_range=c(0.1,1)")+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/theta_pred_vs_true.pdf", plot = theta_pred_vs_true, width = 5, height = 5, units = "in")  

bulk_amp_pred_vs_true = ggplot(df_filt, aes(S_amp_theoret, amplitude_bulk))+
  geom_point(col = "grey32")+ 
  theme_minimal()+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  xlab("Simulated bulk amplitude")+
  ylab("Predicted bulk amplitude")+
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 0.3,\ntheta_range=c(0.1,1)")+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/bulk_amp_pred_vs_true.pdf", plot = bulk_amp_pred_vs_true, width = 5, height = 5, units = "in")  

V2w_amp_pred_vs_true = ggplot(df_filt, aes(V_2w_theoret, double_freq_amplitude))+
  geom_point(col = "grey32")+ 
  theme_minimal()+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  xlab("Simulated 2nd harmonic of variance amp.")+
  ylab("Predicted 2nd harmonic of variance amp.")+
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 0.3,\ntheta_range=c(0.1,1)")+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/V2w_amp_pred_vs_true.pdf", plot = V2w_amp_pred_vs_true, width = 5, height = 5, units = "in")  

Vw_amp_pred_vs_true = ggplot(df_filt, aes(V_w_theoret, fund_amp))+
  geom_point(col = "grey32")+ 
  theme_minimal()+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  xlab("Simulated fundamental of variance amp.")+
  ylab("Predicted fundamental of variance amp.")+
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 0.3,\ntheta_range=c(0.1,1)")+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/Vw_amp_pred_vs_true.pdf", plot = Vw_amp_pred_vs_true, width = 5, height = 5, units = "in")  

custom_labels <- setNames(paste0("phi sd. = ", unique(round(df_filt$phi_true*12/pi, 2)) ), unique(df_filt$phi_true))
acrodiff_compare = ggplot(df_filt, aes(x = acrodiff*12/pi)) +
  geom_histogram(binwidth = 1, color = "black", fill = "steelblue") +
  facet_wrap(~ as.factor(phi_true), scales = "free_y", labeller = as_labeller(custom_labels)) +
  xlab("Fundamental vs 2nd harmonic of variance\n absolute acrophase diff. (hrs)")+
  ylab("Count")+
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 0.3,\ntheta_range=c(0.1,1)")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/acrodiff_compare.pdf", plot = acrodiff_compare, width = 7, height = 7, units = "in")  


# phi4 = filter(df, phi_true == 0.4)
# ggplot() +
#   geom_point(
#     data = phi4,
#     aes(x = S_amp_theoret, y = V_2w_theoret,
#         shape = is_circadian, color = phase_error)
#   ) +
#   labs(
#     x = "Amplitude bulk",
#     y = "Double Frequency Amplitude (input$double_freq_amplitude)"
#   ) +
#   theme_minimal()+
#   scale_colour_gradient(low = "blue", high = "red", na.value = "black", limits = c(-.1, 1))
