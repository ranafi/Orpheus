#' @importFrom foreach foreach `%do%`
#' @importFrom tibble as_tibble
#' @importFrom tidyr drop_na
#' @importFrom stats lm anova p.adjust
NULL

#' Detect Rhythmic Genes in Pseudobulk Data
#'
#' This function identifies cycling genes based on pseudobulk data using cosinor regression.
#'
#' @param times A numeric vector of time points (e.g., ZT times, or radians).
#' @param pb_mat A matrix of pseudobulk gene expression data (genes x samples).
#' @param time_units A string indicating the time units ("hrs" for hours or "radians" for radians). Default is "hrs".
#' @return A tibble summarizing rhythmicity metrics for each gene, including acrophase, amplitude, p-values, and adjusted p-values.

#' @export
is_cycling <- function(times, pb_mat, time_units = "hrs") {
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
  #loop:
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    
    partial_model = lm(gexp1 ~ 1)
    full_model = lm(gexp1 ~ sin(times1) + cos(times1) )
    design_matrix <- model.matrix(gexp1 ~ sin(times1) + cos(times1))
    
    anova_results = anova(partial_model, full_model)
    sin_coff = full_model[["coefficients"]][["sin(times1)"]]
    cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
    acrophase_bulk = atan2(sin_coff, cos_coeff) %% (2*pi)
    
    p_statistic = anova_results$`Pr(>F)`[2]
    Gene_Symbols = colnames(gene1)[gene_i]
    amplitude_bulk = sqrt(sin_coff^2 + cos_coeff^2)
    
    mesor_bulk = full_model$coefficients[["(Intercept)"]]
    
    amp_ratio = amplitude_bulk / mesor_bulk
    
    
    
    gene_summary = cbind( Gene_Symbols, acrophase_bulk,amplitude_bulk, p_statistic, amp_ratio, sin_coff, cos_coeff, mesor_bulk)
    
    return(gene_summary)
    
  }
  all_genes = as_tibble(all_genes) %>% drop_na
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p_statistic), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p_statistic), "bonferroni")
  return(all_genes)
  
}

#' Detect Double-Frequency Rhythms in Variance Data
#'
#' This function tests for double-frequency rhythms in gene expression variance data.
#'
#' @param times A numeric vector of time points (e.g., ZT times).
#' @param var_mat A matrix of variance data (genes x samples).
#' @param mean_mat A matrix of mean expression data (genes x samples).
#' @param time_units A string indicating the time units ("hrs" for hours or "radians" for radians). Default is "hrs".
#' @return A tibble summarizing double-frequency rhythmicity metrics for each gene, including acrophase, amplitude, p-values, and adjusted p-values.

#' @export
var_cycling <- function(times, var_mat, mean_mat, time_units = "hrs") {
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
  
  #loop:
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    m1 = as.numeric(unlist(means1[, gene_i]))
    
    partial_model = lm(gexp1 ~ m1)
    full_model = lm(gexp1 ~ sin(2*times1) + cos(2*times1) + m1)
    design_matrix <- model.matrix(gexp1 ~ sin(2*times1) + cos(2*times1) + m1)
    
    anova_results = anova(partial_model, full_model)
    sin_coff = full_model[["coefficients"]][["sin(2 * times1)"]]
    cos_coeff = full_model[["coefficients"]][["cos(2 * times1)"]]
    double_freq_acrophase = atan2(sin_coff, cos_coeff) %% (2*pi) / 2

    #What is the amplitude of the fundamental according to the double frequency model?
    find_Af = lm(full_model$coefficients[["m1"]]*m1 ~ cos(times1) + sin(times1))
    fund_amp = sqrt(find_Af$coefficients[["sin(times1)"]]^2 + find_Af$coefficients[["cos(times1)"]]^2)
    fund_acro = atan2(find_Af$coefficients[["sin(times1)"]], find_Af$coefficients[["cos(times1)"]]) %% (2*pi)
    
    
    double_freq_p_statistic = anova_results$`Pr(>F)`[2]
    Gene_Symbols = colnames(gene1)[gene_i]
    double_freq_amplitude = sqrt(sin_coff^2 + cos_coeff^2)
    
    #MESOR or mean of the fitted variance wave
    var_mesor = mean( full_model$fitted.values)

    gene_summary = cbind( Gene_Symbols, double_freq_acrophase, double_freq_amplitude, double_freq_p_statistic, fund_amp, fund_acro, var_mesor)
    
    
    
    return(gene_summary)
    
  }
  all_genes = as_tibble(all_genes) %>% drop_na
  all_genes$BHQ = p.adjust(as.numeric(all_genes$double_freq_p_statistic), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$double_freq_p_statistic), "bonferroni")
  return(all_genes)
  
}

#' Plot Rhythms for Selected Genes
#'
#' This function generates plots for pseudobulk expression and double-frequency variance for selected genes.
#'
#' @param times A numeric vector of time points (e.g., ZT times, or radians).
#' @param var_mat A matrix of variance data (genes x samples).
#' @param mean_mat A matrix of mean expression data (genes x samples).
#' @param plot_me A character vector of gene symbols to plot.
#' @param time_units A string indicating the time units ("hrs" for hours or "radians" for radians). Default is "hrs".
#' @return Generates side-by-side plots for pseudobulk expression and double-frequency variance for the selected genes.
#' @export
plot_rhythms <- function(times, var_mat, mean_mat, plot_me = c(), time_units = "hrs") {
  `%do%` <- foreach::`%do%`
  gene1 <- t(var_mat)  # get the transpose, subjects x genes and put in order of CYCLOPS order
  means1 <- t(mean_mat)
  
  if (time_units == "hrs") {
    cat(paste("\nRunning regression, converting hours to radians"))
    times1 <- times / 12 * pi
  } else {
    cat(paste("\nRunning regression assuming times are 0-2pi radians"))
    times1 <- times
  }
  
  # Subset to selected genes, ensuring matrices remain matrices
  gene1 <- gene1[, colnames(gene1) %in% plot_me, drop = FALSE]
  means1 <- means1[, colnames(means1) %in% plot_me, drop = FALSE]

  # Loop through genes
  foreach(gene_i = 1:ncol(gene1), .combine = rbind) %do% {
    gexp1 <- as.numeric(unlist(gene1[, gene_i]))
    m1 <- as.numeric(unlist(means1[, gene_i]))
    Gene_Symbols <- colnames(gene1)[gene_i]
    pb_model_partial = lm(m1~ 1)
    pb_model <- lm(m1 ~ sin(times1) + cos(times1))
    pb_anova_results = anova(pb_model_partial, pb_model)
    pb_acro = atan2(pb_model[["coefficients"]][["sin(times1)"]], pb_model[["coefficients"]][["cos(times1)"]]) %% (2*pi)

    var_partial <- lm(gexp1 ~ m1)
    var_model <- lm(gexp1 ~ sin(2 * times1) + cos(2 * times1) + m1)
    design_matrix <- model.matrix(gexp1 ~ sin(2 * times1) + cos(2 * times1) + m1)
    var_anova_results = anova(var_partial, var_model)
    sin_coff = var_model[["coefficients"]][["sin(2 * times1)"]]
    cos_coeff = var_model[["coefficients"]][["cos(2 * times1)"]]
    double_freq_acrophase = atan2(sin_coff, cos_coeff) %% (2*pi) / 2

    # 2 plots side by side
    par(mfrow = c(1, 2))  # 1 row, 2 columns

    # Plot pseudobulk expression
    plot(times1, m1, xlab = "Time (Rad)", ylab = "PB expression",
         main = paste("PB Expression for", Gene_Symbols),
         sub= paste("p-value:", round(pb_anova_results$`Pr(>F)`[2], 3), ", acrophase:", round(pb_acro, 3)),
         pch = 19, col = "blue")
    lines(sort(times1), pb_model$fitted.values[order(times1)], col = "red")

    # Plot double-frequency variance
    circ <- gexp1 - design_matrix[, -c(2, 3)] %*% var_model[["coefficients"]][-c(2, 3)]
    plot(times1, circ, xlab = "Time (Rad)", ylab = "Double freq component of var",
         main = paste("Double Frequency Variance for", Gene_Symbols),
         sub = paste("p-value:", round(var_anova_results$`Pr(>F)`[2], 3), ", acrophase:", round(double_freq_acrophase, 3)),
         pch = 19,
         col = "blue")
    circ_fit <- var_model$fitted.values - design_matrix[, -c(2, 3)] %*% var_model[["coefficients"]][-c(2, 3)]
    lines(sort(times1), circ_fit[order(times1)], col = "red")
  }
}