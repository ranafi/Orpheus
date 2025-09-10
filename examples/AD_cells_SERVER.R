#Exc_all object too big for local computer, here's non-notebook version of file
#For server use.
## -----------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(Matrix)


file = "exc_all_ctl.h5Seurat"
  
seurat_obj_ctl = LoadH5Seurat(file,assays = c("RNA"), reductions = F, graphs = F, neighbors = F, images = F)
  
VlnPlot(seurat_obj_ctl, features = "nCount_RNA", group.by = "times")

seurat_obj_ctl$seurat_clusters = "0"
seurat_obj_ctl$orig.ident = seurat_obj_ctl$individualID

tmp = Orpheus::extract_clstr_prep_input(seurat_obj_ctl, cluster_id = "0", clip_sd = 15, top_n_genes = 10e3, cells_threshold = 10, normalize = T)


pb_mat_ctl = tmp[[1]]
var_mat_ctl = tmp[[2]]
counts_ctl = tmp[[3]]
#Variance of NB is var = mean + mean^2*theta
#rearrange to solve for theta = (var-mean)/mean^2
time_wise_theta_meas_ctl <- as.data.frame(mapply(function(v,m){(v-m)/m^2}, var_mat_ctl, pb_mat_ctl))
rownames(time_wise_theta_meas_ctl) = rownames(pb_mat_ctl)
colnames(time_wise_theta_meas_ctl) = paste0(colnames(time_wise_theta_meas_ctl), "_theta_meas")

vars_ctl = sparseMatrixStats::rowVars(counts_ctl)
means_ctl = sparseMatrixStats::rowMeans2(counts_ctl)
agg_theta_meas_ctl = data.frame(coarse_theta = (vars_ctl-means_ctl)/means_ctl^2)
rownames(agg_theta_meas_ctl) = rownames(pb_mat_ctl)


times_ctl = seurat_obj_ctl$times[match(colnames(pb_mat_ctl), seurat_obj_ctl$individualID)] %>% unname %>% unlist



## -----------------------------------------------------------------------------------------------------------------------------------

file = "exc_all_ad.h5Seurat"

seurat_obj_ad = LoadH5Seurat(file,assays = c("RNA"), reductions = F, graphs = F, neighbors = F, images = F)
  

seurat_obj_ad$seurat_clusters = "0"
seurat_obj_ad$orig.ident = seurat_obj_ad$individualID

tmp_ad = Orpheus::extract_clstr_prep_input(seurat_obj_ad, cluster_id = "0", clip_sd = 15, top_n_genes = 10e3, cells_threshold = 10, normalize = T)

pb_mat_ad = tmp_ad[[1]]
var_mat_ad = tmp_ad[[2]]
counts_ad = tmp_ad[[3]]
#Variance of NB is var = mean + mean^2*theta
#rearrange to solve for theta = (var-mean)/mean^2
time_wise_theta_meas_ad <- as.data.frame(mapply(function(v,m){(v-m)/m^2}, var_mat_ad, pb_mat_ad))
rownames(time_wise_theta_meas_ad) = rownames(pb_mat_ad)
colnames(time_wise_theta_meas_ad) = paste0(colnames(time_wise_theta_meas_ad), "_theta_meas")
vars_ad = sparseMatrixStats::rowVars(counts_ad)
means_ad = sparseMatrixStats::rowMeans2(counts_ad)
agg_theta_meas_ad = data.frame(coarse_theta = (vars_ad-means_ad)/means_ad^2)
rownames(agg_theta_meas_ad) = rownames(pb_mat_ad)


times_ad = seurat_obj_ad$times[match(colnames(pb_mat_ad), seurat_obj_ad$individualID)] %>% unname %>% unlist




## -----------------------------------------------------------------------------------------------------------------------------------

is_cycling_custom <- function(times, pb_mat,meta_data, time_units = "hrs") {
  `%do%` <- foreach::`%do%`
  cat(paste("\nRunning is_cycling()"))
  gene1 = t(pb_mat)  # get the transpose, subjects x genes and put in order of CYCLOPS order
   if(time_units == "hrs"){
   cat(paste("\nRunning regression, converting hours to radians"))
   times1 = times / 12 * pi
  } else{
    cat(paste("\nRunning regression assuming times are 0-2pi radians"))
    times1 = times
  }
  batches = meta_data$batch_bin[match(colnames(pb_mat), meta_data$individualID)] %>% unname %>% unlist %>% as.factor()
  pmi = meta_data$pmi[match(colnames(pb_mat), meta_data$individualID)] %>% unname %>% unlist %>% as.numeric
  sex = meta_data$sex[match(colnames(pb_mat), meta_data$individualID)] %>% unname %>% unlist %>% as.factor()

    #loop:
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    Gene_Symbols = colnames(gene1)[gene_i]
    # Check for all zeros or zero variance
    if (all(gexp1 == 0) || sd(gexp1) == 0) {
      # Return NA row for this gene
      return(c(Gene_Symbols, rep(NA, 7))) # 9 columns to match the number of outputs
    }
    
    partial_model = lm(gexp1 ~ batches+pmi+sex)
    partial_model_AIC = AIC(partial_model)
    partial_model_BIC = BIC(partial_model)
    
    full_model = lm(gexp1 ~ sin(times1) + cos(times1)+batches+pmi+sex )
    full_model_AIC = AIC(full_model)
    full_model_BIC = BIC(full_model)
    
    design_matrix <- model.matrix(gexp1 ~ sin(times1) + cos(times1)+batches+pmi+sex)
    
    anova_results = anova(partial_model, full_model)
    sin_coff = full_model[["coefficients"]][["sin(times1)"]]
    cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
    acrophase_bulk = atan2(sin_coff, cos_coeff) %% (2*pi)
    
    p_statistic = anova_results$`Pr(>F)`[2]
    amplitude_bulk = sqrt(sin_coff^2 + cos_coeff^2)
    
    mesor_bulk = mean(design_matrix[,-c(2, 3)] %*% full_model[["coefficients"]][-c(2, 3)]) 
    # mesor_bulk = full_model$coefficients[["(Intercept)"]]
    
    amp_ratio = amplitude_bulk / mesor_bulk
    
    
    
   gene_summary = cbind( Gene_Symbols, acrophase_bulk,amplitude_bulk,
              p_statistic, amp_ratio, sin_coff, cos_coeff, mesor_bulk,
              partial_model_AIC, full_model_AIC, partial_model_BIC, full_model_BIC)
    
    return(gene_summary)
    
  }
  all_genes = as_tibble(all_genes) %>% drop_na
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p_statistic), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p_statistic), "bonferroni")
  return(all_genes)
  
}

var_cycling_custom <- function(times, var_mat, mean_mat, phase_pb,meta_data, time_units = "hrs") {
  `%do%` <- foreach::`%do%`
  cat(paste("\nRunning var_cycling()"))
  gene1 = t(var_mat)  # get the transpose, subjects x genes and put in order of CYCLOPS order
  means1 = t(mean_mat)
  if(time_units == "hrs"){
    cat(paste("\nRunning regression, converting hours to radians"))
    times1 = times / 12 * pi
  } else{
    cat(paste("\nRunning regression assuming times are 0-2pi radians"))
    times1 = times
  }
  batches = meta_data$batch_bin[match(colnames(mean_mat), meta_data$individualID)] %>% unname %>% unlist %>% as.factor()
  pmi = meta_data$pmi[match(colnames(mean_mat), meta_data$individualID)] %>% unname %>% unlist %>% as.numeric
  sex = meta_data$sex[match(colnames(mean_mat), meta_data$individualID)] %>% unname %>% unlist %>% as.factor()
  
  #loop:
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    m1 = as.numeric(unlist(means1[, gene_i]))
    Gene_Symbols = colnames(gene1)[gene_i]
    phase_bulk = as.numeric(phase_pb[gene_i])


    # Check for all zeros or zero variance
    if (all(m1 == 0) || all(gexp1 == 0) || sd(m1) == 0 || sd(gexp1) == 0) {
      # Return NA row for this gene
      return(c(Gene_Symbols, rep(NA, 8))) # 9 columns to match the number of outputs
    }
    
    partial_model = lm(gexp1 ~ m1 + batches+sex+pmi)
    partial_model_AIC = AIC(partial_model)
    partial_model_BIC = BIC(partial_model)

    full_model = lm(gexp1 ~ sin(2*times1) + cos(2*times1) + m1 + batches+pmi+sex)
    full_model_AIC = AIC(full_model)
    full_model_BIC = BIC(full_model)

    design_matrix <- model.matrix(gexp1 ~ sin(2*times1) + cos(2*times1) + m1+ batches+pmi+sex)
    
    anova_results = anova(partial_model, full_model)
    sin_coff = full_model[["coefficients"]][["sin(2 * times1)"]]
    cos_coeff = full_model[["coefficients"]][["cos(2 * times1)"]]
    double_freq_acrophase = atan2(sin_coff, cos_coeff) %% (2*pi) / 2

    #What is the amplitude of the fundamental according to the double frequency model?
    find_Af = lm(full_model$coefficients[["m1"]]*m1 ~ cos(times1) + sin(times1) + batches+pmi+sex)
    fund_amp = sqrt(find_Af$coefficients[["sin(times1)"]]^2 + find_Af$coefficients[["cos(times1)"]]^2)
    fund_acro = atan2(find_Af$coefficients[["sin(times1)"]], find_Af$coefficients[["cos(times1)"]]) %% (2*pi)
    
    
    double_freq_p_statistic = anova_results$`Pr(>F)`[2]
    raw_2w_amp = sqrt(sin_coff^2 + cos_coeff^2)

    # #Calculate the 95% CI for double frequency phase
    # if(double_freq_amplitude > 1e-3){
    #   ## Extract coefficients and covariance ---------------
    #   b  <- coef(full_model)["cos(2 * times1)"]   # β̂
    #   g  <- coef(full_model)["sin(2 * times1)"]   # γ̂
    #   V  <- vcov(full_model)[c("cos(2 * times1)","sin(2 * times1)"),
    #                   c("cos(2 * times1)","sin(2 * times1)")]

    #   ## SE via delta method -------------------------------
    #   var.phi <- (g^2*V[1,1] + b^2*V[2,2] - 2*b*g*V[1,2]) / (b^2+g^2)^2
    #   se.phi  <- sqrt(var.phi)
    #   ## 95 % CI (in radians) ------------------------------
    #   phi.hat <- atan2(g, b)
    #   ci.phi  <- phi.hat + c(-1, 1)*qnorm(0.975)*se.phi
    # }
    
    #Amp of double frequency variance wave projected onto expected bulk wave.
    double_freq_amplitude = -cos_coeff * cos(2 * phase_bulk) - sin_coff * sin(2 * phase_bulk)

    #MESOR or mean of the fitted variance wave
    # var_mesor = mean( full_model$fitted.values)
    var_mesor = mean(design_matrix[,-c(2, 3)] %*% full_model[["coefficients"]][-c(2, 3)])

    gene_summary = cbind( Gene_Symbols, double_freq_acrophase, double_freq_amplitude,
          double_freq_p_statistic, fund_amp, fund_acro, var_mesor,
          partial_model_AIC, full_model_AIC, partial_model_BIC, full_model_BIC, raw_2w_amp)
    
    
    
    return(gene_summary)
    
  }
  all_genes = as_tibble(all_genes) %>% drop_na
  all_genes$BHQ = p.adjust(as.numeric(all_genes$double_freq_p_statistic), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$double_freq_p_statistic), "bonferroni")
  return(all_genes)
  
}


## -----------------------------------------------------------------------------------------------------------------------------------
meta_ad = seurat_obj_ad@meta.data
meta_ctl = seurat_obj_ctl@meta.data

pb_cycling_ad_custom = is_cycling_custom(times_ad, pb_mat_ad, meta_ad, time_units = "rads")
pb_cycling_ctl_custom = is_cycling_custom(times_ctl, pb_mat_ctl, meta_ctl, time_units = "rads")

all(pb_cycling_ad_custom$Gene_Symbols == rownames(var_mat_ad))
pb_phases_cust_ad = as.numeric(pb_cycling_ad_custom$acrophase_bulk)

all(pb_cycling_ctl_custom$Gene_Symbols == rownames(var_mat_ctl))
pb_phases_cust_ctl = as.numeric(pb_cycling_ctl_custom$acrophase_bulk)

#Test how different answers are
# test = merge(pb_cycling_ad_custom, pb_cycling_ad, by = "Gene_Symbols", suffixes = c("_cust", "_orig")) %>% filter(as.numeric(p_statistic_cust) < 0.05 & as.numeric(p_statistic_orig) < 0.05)
# plot(test$amp_ratio_cust, test$amp_ratio_orig)


var_cyclers_cust_ad = var_cycling_custom(times_ad, var_mat = var_mat_ad, mean_mat = pb_mat_ad, phase_pb = pb_phases_cust_ad, meta_data = meta_ad, time_units = "rads")

var_cyclers_cust_ctl = var_cycling_custom(times_ctl, var_mat = var_mat_ctl, mean_mat = pb_mat_ctl, phase_pb = pb_phases_cust_ctl,meta_data = meta_ctl,  time_units = "rads")

#Test how different answers are
# test = merge(var_cyclers_cust_ad, var_cyclers_ad, by = "Gene_Symbols", suffixes = c("_cust", "_orig")) %>% filter(double_freq_p_statistic_cust < 0.05 & double_freq_p_statistic_orig < 0.05)
# plot(as.numeric(test$double_freq_acrophase_cust), as.numeric(test$double_freq_acrophase_orig))

bulk_thr_cutoff = 0.5

input_cos_ad_cust = var_cyclers_cust_ad %>% 
    left_join(pb_cycling_ad_custom, by = "Gene_Symbols", suffix = c("_var", "_pb")) %>%
    mutate(across(-Gene_Symbols, as.numeric)) %>%
    filter(BHQ_pb < 0.1 & amp_ratio >= 0.2) %>%
    mutate(acrodiff = (double_freq_acrophase - acrophase_bulk) %% pi) %>%
    dplyr::select(Gene_Symbols, double_freq_amplitude,fund_amp,amplitude_bulk, mesor_bulk,  everything()) %>%
    filter(amplitude_bulk >= bulk_thr_cutoff) %>%
    column_to_rownames("Gene_Symbols")%>%
    filter(partial_model_AIC_var > full_model_AIC_var)

# Calculate phi parameters
phi_ad_cust = input_cos_ad_cust %>%
    apply(., 1, function(row) {
        suppressWarnings(
             Orpheus::find_analytic(row[1], row[2], row[3], row[4])
             )
    }) %>%
    t %>%
    as.data.frame
colnames(phi_ad_cust) = c("G", "Ac", "theta", "phi")


input_cos_ctl_cust = var_cyclers_cust_ctl %>% 
    left_join(pb_cycling_ctl_custom, by = "Gene_Symbols", suffix = c("_var", "_pb")) %>%
    mutate(across(-Gene_Symbols, as.numeric)) %>%
    filter(BHQ_pb < 0.1 & amp_ratio >= 0.2) %>%
    mutate(acrodiff = (double_freq_acrophase - acrophase_bulk) %% pi) %>%
    dplyr::select(Gene_Symbols, double_freq_amplitude,fund_amp,amplitude_bulk, mesor_bulk,  everything()) %>%
    filter(amplitude_bulk >= bulk_thr_cutoff) %>%
    column_to_rownames("Gene_Symbols")%>%
    filter(partial_model_AIC_var > full_model_AIC_var)

# Calculate phi parameters
phi_ctl_cust = input_cos_ctl_cust %>%
    apply(., 1, function(row) {
        suppressWarnings(
             Orpheus::find_analytic(row[1], row[2], row[3], row[4])
             )
    }) %>%
    t %>%
    as.data.frame
colnames(phi_ctl_cust) = c("G", "Ac", "theta", "phi")

png(filename="~/orpheus/Exc_all_ctl_cycling_amps.png")
pb_cycling_ctl_custom %>% filter(BHQ < 0.1 & amp_ratio >= 0.2) %>% select(amplitude_bulk) %>%
  unlist %>% unname %>% as.numeric %>%
  hist(, xlim = c(0,2), breaks = seq(0, 100, by = 0.2), main = "bulk amplitude Exc all ctl cells")
dev.off()

png(filename="~/orpheus/Exc_all_ad_cycling_amps.png")
pb_cycling_ad_custom %>% filter(BHQ < 0.1 & amp_ratio >= 0.2) %>%
  select(amplitude_bulk) %>% unlist %>% unname %>% as.numeric %>%
  hist(, xlim = c(0,2), breaks = seq(0, 100, by = 0.2), main = "bulk amplitude Exc all ad cells")
dev.off()

print("dim ctl cts:")
print(dim(counts_ctl))
print("dim ad cts:")
print(dim(counts_ad))

## -----------------------------------------------------------------------------------------------------------------------------------
wilcox.test(phi_ctl_cust$phi, phi_ad_cust$phi)
boxplot(phi_ctl_cust$phi, phi_ad_cust$phi)

merged_phi = merge(phi_ctl_cust, phi_ad_cust, by = 0, suffixes = c("_CTL", "_AD")) %>% left_join(select(pb_cycling_ad_custom, Gene_Symbols, amplitude_bulk), by= join_by("Row.names"=="Gene_Symbols")) %>% rename("amplitude_bulk_AD"="amplitude_bulk") %>% left_join(select(pb_cycling_ctl_custom, Gene_Symbols, amplitude_bulk), by= join_by("Row.names"=="Gene_Symbols")) %>% rename("amplitude_bulk_CTL"="amplitude_bulk")

plot_df_long = pivot_longer(merged_phi, cols = -Row.names, names_to = c(".value", "condition"), names_pattern = "(.*)_(AD|CTL)") %>% mutate(condition = as.factor(condition))

custom_colors <- c("AD" = "#cc3700", "CTL" = "#0094cc")

matched_phi_boxplot = ggplot(plot_df_long, aes(condition, phi*12/pi, fill = condition))+
 geom_boxplot(outliers = F)+
  theme_minimal()+
  xlab("Condition")+
  ylab("Prediceted phase std. (hrs)")+
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = custom_colors) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.6),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12),   # Increase legend title size
    legend.position = "none"
  )

ggsave("~/orpheus/matched_phi_boxplot.pdf",matched_phi_boxplot, device = "pdf", width = 5, height = 5 )

wilcox.test(merged_phi$phi_CTL, merged_phi$phi_AD, paired = T)
boxplot(merged_phi$phi_CTL, merged_phi$phi_AD)


matches_Ac_boxplot = ggplot(plot_df_long, aes(condition, log(Ac), fill = condition))+
 geom_boxplot(outlier.shape = NA)+
  theme_minimal()+
  xlab("Condition")+
  ylab("Log(cell rhythm amp.)")+
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = custom_colors) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.6),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12),   # Increase legend title size
    legend.position = "none"
  )

ggsave("~/orpheus/matched_Ac_boxplot.pdf",matches_Ac_boxplot, device = "pdf", width = 5, height = 5 )

matches_pbAmp_boxplot = ggplot(plot_df_long, aes(condition, log(as.numeric(amplitude_bulk)), fill = condition))+
 geom_boxplot(outliers  = F)+
  theme_minimal()+
  xlab("Condition")+
  ylab("Log(PB rhythm amp.)")+
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = custom_colors) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.6),    # Increase x-axis text size
    axis.text.y = element_text(size = 14),    # Increase y-axis text size
    axis.title.x = element_text(size = 16),   # Increase x-axis title size
    axis.title.y = element_text(size = 16),   # Increase y-axis title size
    legend.title = element_text(size = 12),   # Increase legend title size
    legend.position = "none"
  )

ggsave("~/orpheus/matched_pbAmp_boxplot.pdf",matches_pbAmp_boxplot, device = "pdf", width = 5, height = 5 )

save(merged_phi, file = "~/orpheus/Exc_all_merged_phi.RData")
save(phi_ctl_cust, file = "~/orpheus/Exc_all_phi_ctl_cust.RData")
save(phi_ad_cust, file = "~/orpheus/Exc_all_phi_ad_cust.RData")

## -----------------------------------------------------------------------------------------------------------------------------------
out_rnk = data.frame(Gene_Symbols = merged_phi$Row.names,
                     metric = merged_phi$phi_CTL - merged_phi$phi_AD) %>% drop_na() %>% arrange(metric)

write.table(out_rnk, file = "~/orpheus/rosmap_exc_all_ctlPhiMinusADPhi.rnk", 
                  row.names = FALSE,
                  col.names = FALSE,
                  sep = "\t",
                  quote = F)
      

