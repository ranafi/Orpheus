library(Seurat)
library(dplyr)
# setwd("~/Box Sync/Henry_stuff/Orpheus/Orpheus/R")
devtools::check()
library(Orpheus)

simulate_scRNAseq_circadian <- function(n = 10000, m = 5000, p_circ = 0.3,n_time_points = 6, 
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
      # repeat_nums <- table(t_vals)
      # cell_specific_epsilon_noise <- rep(rnorm(n_time_points, 0, epsilon_sd), times = repeat_nums)
      
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
      # repeat_nums <- table(t_vals)
      # cell_specific_epsilon_noise <- rep(rnorm(n_time_points, 0, epsilon_sd), times = repeat_nums)
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

output = list()
for(i in seq(0, 1.3, by = 0.1)){
  print(paste0("phi = ", i))
  phi_sd = i
  epsilon_sd = 0.3
  sim <- simulate_scRNAseq_circadian(n = 5000, m = 10000, p_circ = .2, phi_sd = phi_sd, n_time_points = 12,
      A_range = c(.1, 5), G_range = c(.1, 5), epsilon_sd = epsilon_sd, theta_range = c(0.1, 1))
  # sim <- simulate_scRNAseq_circadian(n = 50000, m = 2000, p_circ = .5, phi_sd = phi_sd, n_time_points = 24,
  #                                    A_range = c(.1, 5), G_range = c(.1, 5), epsilon_sd = epsilon_sd, theta_range = c(0.1, .2))
  counts <- sim$counts
  metadata <- sim$cell_metadata
  gene_info <- sim$gene_parameters
  
  seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
  seurat_obj$orig.ident = seurat_obj$time_point
  seurat_obj$seurat_clusters = "0"
  
  tmp = extract_clstr_prep_input(seurat_obj, cluster_id = "0", normalize = T)
  
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

save(df, file = "~/Desktop/df_50kcells10kgenes_12times_epsilon03.RData")

#Keep rows of results where bulk amp > thresh & not NaN
df_filt = dplyr::filter(df, !is.na(phi) & amplitude_bulk >= 0.5)

#Plot: Phi sd pred acc vs true phi sd
df_filt$phi_true_factor = as.factor(df_filt$phi_true*12/pi)
phi_acc_plot = ggplot(df_filt, aes(x = phi_true*12/pi,y =  phi*12/pi))+
  geom_boxplot(aes(group = phi_true_factor))+
  # geom_jitter(width = 0.1, alpha = 0.5) +
  # geom_violin(aes(group = phi_true_factor)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  ggtitle("n=5k,m=10k, p_circ=.2, phi_sd=phi_sd, n_time_points=12,
    A_range=c(.1,5),G_range=c(.1,5),epsilon_sd=0.3,theta_range=c(0.1,1)")+
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

helper = function(df){
  isnan = length(which(is.na(df$phi)))/ nrow(df)
  total = nrow(df)
  med_pred_phi = median(df$phi, na.rm = T)
  return(list(isnan = isnan, total = total, pred_phi = med_pred_phi))
}
frac_nan = list()
for(i in seq(0, 5, by = 0.1)){
  for(j in unique(df$phi_true)){
    res = df %>% filter(as.numeric(amplitude_bulk) >= i & phi_true == j) %>% helper()
    res = c(res, i, j)
    names(res) = c("isnan", "total","median_pred_phi", "amp_bulk", "phi_true")
    frac_nan = append(frac_nan, list(res))
  }
}
nan_plot <- do.call(rbind, lapply(seq_along(frac_nan), function(i) {
  current_df <- frac_nan[[i]]
  unlist(current_df)
})) %>% as.data.frame %>% drop_na()

nan_plot_long <- nan_plot %>%
  mutate(true_phi = (phi_true)) %>%
  pivot_longer(cols = c(total, isnan, median_pred_phi, true_phi),
               names_to = "metric", values_to = "value")%>%
  mutate(group = ifelse(metric %in% c("true_phi", "median_pred_phi"), "phi_plot", metric))

custom_labels <- setNames(paste0("phi sd. = ", unique(round(df$phi_true, 2)) ), unique(df$phi_true))
ggplot(nan_plot_long, aes(x = amp_bulk, y = value, color = metric)) +
  geom_line(size = 1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey")+
  xlim(c(0,2))+
  facet_grid(group ~ phi_true, labeller = as_labeller(custom_labels), scales = "free_y") +
  scale_color_manual(values = c("isnan" = "blue", "total" = "green", "median_pred_phi" = "red", "true_phi" = "black")) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 14, angle = 90),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )

ggsave("~/Desktop/percent_nans_vs_bulk_amp.pdf", device = "pdf", width = 15, height = 10)

# model <- lm(abs(phase_error) ~ S_amp_theoret, data = df)
# summary_model <- summary(model)
# slope <- coef(model)[["S_amp_theoret"]]
# p_value <- summary_model$coefficients["S_amp_theoret", "Pr(>|t|)"]
# r_squared <- summary_model$r.squared
# label_text <- sprintf("Slope = %.3f\np = %.3g\nR² = %.3f", slope, p_value, r_squared)
bulk_amp_vs_phase_err = ggplot(df, aes(amplitude_bulk, abs(phase_error*12/pi)))+
  geom_point(col = "grey32")+ 
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal()+
  # annotate("text", x = Inf, y = -Inf, label = label_text, 
  #          hjust = 1.1, vjust = -0.5, size = 4, color = "black") +
  # ylim(0, 2)+
  ylab("Absolute error of phase std. (hrs)")+
  xlab("bulk amplitude (measured)")+
  ggtitle("n=50000,m=2000, p_circ=.5, phi_sd=phi_sd, n_time_points=24,
    A_range=c(.1,5),G_range=c(.1,5),epsilon_sd=0.3,theta_range=c(0.1,.2)")+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/bulk_amp_vs_phase_err.pdf", plot = bulk_amp_vs_phase_err, width = 5, height = 5, units = "in")  

theta_pred_vs_true = ggplot(df_filt, aes(theta_true, theta))+
  geom_point(col = "grey32")+ 
  theme_minimal()+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  xlab("Simulated overdispersion param.")+
  ylab("Predicted overdispersion param.")+
  ggtitle("n=50000,m=2000, p_circ=.5, phi_sd=phi_sd, n_time_points=24,
    A_range=c(.1,5),G_range=c(.1,5),epsilon_sd=0.3,theta_range=c(0.1,.2)")+
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
  ggtitle("n=50000,m=2000, p_circ=.5, phi_sd=phi_sd, n_time_points=24,
    A_range=c(.1,5),G_range=c(.1,5),epsilon_sd=0.3,theta_range=c(0.1,.2)")+
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
  ggtitle("n=50000,m=2000, p_circ=.5, phi_sd=phi_sd, n_time_points=24,
    A_range=c(.1,5),G_range=c(.1,5),epsilon_sd=0.3,theta_range=c(0.1,.2)")+
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
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 2.0,\ntheta_range=c(0.1,1)")+
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
  ggtitle("n=5000,m=10000,p_circ=.2,n_time_points = 12,\nA_range=c(.1, 5),G_range = c(.1, 5),epsilon_sd = 2.0,\ntheta_range=c(0.1,1)")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12)    # Increase legend title size
  )
ggsave("~/Desktop/acrodiff_compare.pdf", plot = acrodiff_compare, width = 7, height = 7, units = "in")  

################ Testing how crude dispersion compares to my theta_estimate ###########
phi_sd = 0.7
epsilon_sd = 0.3
sim <- simulate_scRNAseq_circadian(n = 5000, m = 10000, p_circ = .2, phi_sd = phi_sd, n_time_points = 24,
                                   A_range = c(.1, 5), G_range = c(.1, 5), epsilon_sd = epsilon_sd, theta_range = c(0.1, 1))
counts <- sim$counts
metadata <- sim$cell_metadata
gene_info <- sim$gene_parameters

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
seurat_obj$orig.ident = seurat_obj$time_point
seurat_obj$seurat_clusters = "0"

tmp = extract_clstr_prep_input(seurat_obj, cluster_id = "0", normalize = T)

pb_mat = tmp[[1]]
var_mat = tmp[[2]]
counts = tmp[[3]]

#Variance of NB is var = mean + mean^2*theta
#rearrange to solve for theta = (var-mean)/mean^2
# time_wise_theta_meas <- as.data.frame(mapply(function(v,m){(v-m)/m^2}, var_mat, pb_mat))
# rownames(time_wise_theta_meas) = rownames(pb_mat)
# colnames(time_wise_theta_meas) = paste0(colnames(time_wise_theta_meas), "_theta_meas")
vars = sparseMatrixStats::rowVars(counts)
means = sparseMatrixStats::rowMeans2(counts)
agg_theta_meas = data.frame(coarse_theta = (vars-means)/means^2, var = vars, mean = means)
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
  filter( amplitude_bulk >= 0.5) %>%
  # select(c(1:5)) %>%
  filter(partial_model_AIC_var > full_model_AIC_var)

plot_df = merge(input, phi, by = 0) %>%
  merge(., gene_info, by.x="Row.names", by.y = 0)%>%
  mutate(phase_error = (phi -phi_sd ))
plot_df$phi_true = phi_sd
# out = merge(plot_df, time_wise_theta_meas, by.x = "Row.names", by.y = 0)
out = merge(plot_df, agg_theta_meas, by.x = "Row.names", by.y = 0)

plot(out$coarse_theta, out$theta_true)
lines(c(0,100), c(0,100))
plot(out$theta, out$theta_true)
lines(c(0,100), c(0,100))
out = out %>% mutate(theta_inflation_factor = (1+theta_true)*(A_true^2/2 + epsilon_sd^2)/(G_true)^2 )
plot(out$coarse_theta - out$theta_inflation_factor, out$theta_true)
lines(c(0,100), c(0,100))
############### Testing the how variance changes with num replicates ####################
sim <- simulate_scRNAseq_circadian(n = 5000, m = 10000, p_circ = .2, phi_sd = 0.7, n_time_points = 12,
                                   A_range = c(.1, 5), G_range = c(.1, 5), epsilon_sd = 0.5, theta_range = c(0.1, 1))
counts <- sim$counts
metadata <- sim$cell_metadata
gene_info <- sim$gene_parameters

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
seurat_obj$orig.ident = seurat_obj$time_point
seurat_obj$seurat_clusters = "0"

tmp = extract_clstr_prep_input2(seurat_obj, cluster_id = "0", clip_sd = 10, normalize = T)

pb_mat = tmp[[1]]
var_mat = tmp[[2]]

# Reshape ground truth data into a long format for easy comparison
truth_mean_long <- pb_mat %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "orig.ident", values_to = "truth_mean")

truth_var_long <- var_mat %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "orig.ident", values_to = "truth_var")

ground_truth_df <- full_join(truth_mean_long, truth_var_long, by = c("gene", "orig.ident"))


results_list <- list()
for(i in seq(2, 10, by = 2)){
  message(paste("Processing for min_cells_per_rep =", i))
  rep_tmp = extract_clstr_prep_input2(seurat_obj, cluster_id = "0", clip_sd = 10, normalize = T, num_replicates =i)
  
  pb_mat_rep = rep_tmp[[1]]
  var_mat_rep = rep_tmp[[2]]
  # Reshape replicate data into a long format
  rep_mean_long <- pb_mat_rep %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "replicate_id", values_to = "rep_mean") %>%
    separate(replicate_id, into = c("orig.ident", "rep_num"), sep = "_rep", remove = FALSE)
  
  rep_var_long <- var_mat_rep %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "replicate_id", values_to = "rep_var") %>%
    separate(replicate_id, into = c("orig.ident", "rep_num"), sep = "_rep", remove = FALSE)
  
  # Calculate average mean and variance across replicates for each subject
  avg_rep_stats <- full_join(rep_mean_long, rep_var_long, by = c("gene", "replicate_id", "orig.ident", "rep_num")) %>%
    group_by(gene, orig.ident) %>% #group together the pseudoreps
    summarise(
      avg_rep_mean = mean(rep_mean, na.rm = TRUE),
      avg_rep_var = mean(rep_var, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate Cells/Replicate and Quantify Departure 
  
  # Get the number of cells per replicate
  total_cells_per_subject <- table(seurat_obj$orig.ident)
  avg_cells_per_rep <- mean(total_cells_per_subject / i)
  
  # Join with ground truth to calculate metrics
  comparison_df <- left_join(avg_rep_stats, ground_truth_df, by = c("gene", "orig.ident")) %>%
    mutate(percent_err_var = 100*(avg_rep_var - truth_var)/truth_var) %>%
    mutate(percent_err_mean = 100*(avg_rep_mean - truth_mean)/truth_mean)
  
  # Bin genes into Low, Middle, and High expression groups based on the ground truth mean
  # We do this within each subject to account for subject-specific expression shifts
  binned_metrics <- comparison_df %>%
    filter(!is.na(truth_mean)) %>%
    group_by(orig.ident) %>%
    mutate(
      expression_bin = ntile(truth_mean, 3), # Create 3 quantile bins (tertiles)
      expression_bin = factor(expression_bin, levels = 1:3, labels = c("Low", "Middle", "High"))
    ) %>%
    ungroup() %>%
    # Now calculate metrics PER BIN for each subject
    group_by(orig.ident, expression_bin) %>%
    summarise(
      rmse_mean = sqrt(mean((avg_rep_mean - truth_mean)^2, na.rm = TRUE)),
      rmse_var = sqrt(mean((avg_rep_var - truth_var)^2, na.rm = TRUE)),
      cor_mean = cor(avg_rep_mean, truth_mean, use = "pairwise.complete.obs"),
      cor_var = cor(avg_rep_var, truth_var, use = "pairwise.complete.obs"),
      .groups = "drop"
    ) %>%
    # Finally, average the metrics across all subjects for each bin
    group_by(expression_bin) %>%
    summarise(
      mean_rmse_mean = mean(rmse_mean, na.rm = TRUE),
      mean_rmse_var = mean(rmse_var, na.rm = TRUE),
      mean_cor_mean = mean(cor_mean, na.rm = TRUE),
      mean_cor_var = mean(cor_var, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Store results
  results_list[[as.character(i)]] <- binned_metrics %>%
    mutate(
      num_replicates = i,
      avg_cells_per_replicate = avg_cells_per_rep)
  
}


# Combine all results into a single data frame
final_results_df <- bind_rows(results_list)



# Plot RMSE vs. Cells per Replicate, faceted by expression bin
plot_rmse_binned <- ggplot(final_results_df, aes(x = avg_cells_per_replicate)) +
  geom_line(aes(y = mean_rmse_mean, color = "Mean")) +
  geom_point(aes(y = mean_rmse_mean, color = "Mean")) +
  geom_line(aes(y = mean_rmse_var, color = "Variance")) +
  geom_point(aes(y = mean_rmse_var, color = "Variance")) +
  facet_wrap(~ expression_bin, scales = "free_y") + # Use free y-axis scales
  scale_x_reverse() +
  labs(
    title = "Departure of Pseudoreplicate Statistics by Expression Level (RMSE)",
    x = "min_cells_per_replicate",
    y = "Root Mean Squared Error (RMSE)",
    color = "Statistic"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
