#' @importFrom foreach foreach `%do%`
#' @importFrom tibble as_tibble
#' @importFrom tidyr drop_na all_of
#' @importFrom stats lm anova p.adjust AIC qnorm vcov coef
#' @importFrom utils write.table
#' @importFrom MetaCycle meta2d
#' @importFrom tools file_path_sans_ext
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs
#' @importFrom patchwork wrap_plots
#' @importFrom stats BIC
#' @importFrom utils read.table
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
    Gene_Symbols = colnames(gene1)[gene_i]

      # Check for all zeros or zero variance
    if (all(gexp1 == 0) || sd(gexp1) == 0) {
      # Return NA row for this gene
      return(c(Gene_Symbols, rep(NA, 11))) #  columns to match the number of outputs
    }
    
    partial_model = lm(gexp1 ~ 1)
    partial_model_AIC = AIC(partial_model)
    partial_model_BIC = BIC(partial_model)

    full_model = lm(gexp1 ~ sin(times1) + cos(times1) )
    full_model_AIC = AIC(full_model)
    full_model_BIC = BIC(full_model)

    design_matrix <- model.matrix(gexp1 ~ sin(times1) + cos(times1))
    sv <- svd(design_matrix)$d
    tol <- 1e-14
    rankX <- sum(sv > tol)
    full <- ncol(design_matrix)

    if (rankX < full) {
      # design matrix is rank-deficient -> at least one coef would be NA
      return(c(Gene_Symbols, rep(NA, 11))) # columns to match the number of outputs
    }
    
    anova_results = anova(partial_model, full_model)
    sin_coeff = full_model[["coefficients"]][["sin(times1)"]]
    cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
    acrophase_bulk = atan2(sin_coeff, cos_coeff) %% (2*pi)
    
    p_statistic = anova_results$`Pr(>F)`[2]
    amplitude_bulk = sqrt(sin_coeff^2 + cos_coeff^2)
    
    mesor_bulk = full_model$coefficients[["(Intercept)"]]
    
    amp_ratio = amplitude_bulk / mesor_bulk
    
    
    
    gene_summary = cbind( Gene_Symbols, acrophase_bulk,amplitude_bulk,
              p_statistic, amp_ratio, sin_coeff, cos_coeff, mesor_bulk,
              partial_model_AIC, full_model_AIC, partial_model_BIC, full_model_BIC)
    
    return(gene_summary)
    
  }
  all_genes = as_tibble(all_genes)
  colnames(all_genes) <- c(
    "Gene_Symbols", "acrophase_bulk","amplitude_bulk",
    "p_statistic", "amp_ratio", "sin_coeff", "cos_coeff", "mesor_bulk",
    "partial_model_AIC", "full_model_AIC", "partial_model_BIC", "full_model_BIC")
  all_genes = drop_na(all_genes)
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
#' @param phase_pb A numeric vector of pseudobulk acrophases (in radians) for each gene, typically obtained from `is_cycling()`.
#' @param time_units A string indicating the time units ("hrs" for hours or "radians" for radians). Default is "hrs".
#' @return A tibble summarizing double-frequency rhythmicity metrics for each gene, including acrophase, amplitude, p-values, and adjusted p-values.

#' @export
var_cycling <- function(times, var_mat, mean_mat, phase_pb, time_units = "hrs") {
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
    Gene_Symbols = colnames(gene1)[gene_i]
    phase_bulk = as.numeric(phase_pb[gene_i])


    # Check for all zeros or zero variance
    if (all(m1 == 0) || all(gexp1 == 0) || sd(m1) == 0 || sd(gexp1) == 0) {
      # Return NA row for this gene
      return(c(Gene_Symbols, rep(NA, 11))) # columns to match the number of outputs
    }
    
    partial_model = lm(gexp1 ~ m1)
    partial_model_AIC = AIC(partial_model)
    partial_model_BIC = BIC(partial_model)

    full_model = lm(gexp1 ~ sin(2*times1) + cos(2*times1) + m1)
    full_model_AIC = AIC(full_model)
    full_model_BIC = BIC(full_model)

    design_matrix <- model.matrix(gexp1 ~ sin(2*times1) + cos(2*times1) + m1)
    sv <- svd(design_matrix)$d
    tol <- 1e-14
    rankX <- sum(sv > tol)
    full <- ncol(design_matrix)

    if (rankX < full) {
      # design matrix is rank-deficient -> at least one coef would be NA
      return(c(Gene_Symbols, rep(NA, 11))) # columns to match the number of outputs
    }
    
    anova_results = anova(partial_model, full_model)
    sin_coeff = full_model[["coefficients"]][["sin(2 * times1)"]]
    cos_coeff = full_model[["coefficients"]][["cos(2 * times1)"]]
    double_freq_acrophase = atan2(sin_coeff, cos_coeff) %% (2*pi) / 2

    #What is the amplitude of the fundamental according to the double frequency model?
    find_Af = lm(full_model$coefficients[["m1"]]*m1 ~ cos(times1) + sin(times1))
    fund_amp = sqrt(find_Af$coefficients[["sin(times1)"]]^2 + find_Af$coefficients[["cos(times1)"]]^2)
    fund_acro = atan2(find_Af$coefficients[["sin(times1)"]], find_Af$coefficients[["cos(times1)"]]) %% (2*pi)
    
    
    double_freq_p_statistic = anova_results$`Pr(>F)`[2]
    raw_2w_amp = sqrt(sin_coeff^2 + cos_coeff^2)

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
    double_freq_amplitude = -cos_coeff * cos(2 * phase_bulk) - sin_coeff * sin(2 * phase_bulk)

    #MESOR or mean of the fitted variance wave
    # var_mesor = mean( full_model$fitted.values)
    var_mesor = mean(design_matrix[,-c(2, 3)] %*% full_model[["coefficients"]][-c(2, 3)])

    gene_summary = cbind( Gene_Symbols, double_freq_acrophase, double_freq_amplitude,
          double_freq_p_statistic, fund_amp, fund_acro, var_mesor,
          partial_model_AIC, full_model_AIC, partial_model_BIC, full_model_BIC, raw_2w_amp)
    
    
    
    return(gene_summary)
    
  }
  all_genes = as_tibble(all_genes)
  colnames(all_genes) <- c(
    "Gene_Symbols", "double_freq_acrophase", "double_freq_amplitude",
    "double_freq_p_statistic", "fund_amp", "fund_acro", "var_mesor",
    "partial_model_AIC", "full_model_AIC", "partial_model_BIC", "full_model_BIC", "raw_2w_amp"
  )
  all_genes = drop_na(all_genes)
  all_genes$BHQ = p.adjust(as.numeric(all_genes$double_freq_p_statistic), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$double_freq_p_statistic), "bonferroni")
  return(all_genes)
  
}

#' Detect Double-Frequency Rhythms in Variance Data Given Bulk Acrophase from `is_cycling()`
#'
#' This function tests for double-frequency rhythms in gene expression variance data.
#'
#' @param times A numeric vector of time points (e.g., ZT times).
#' @param var_mat A matrix of variance data (genes x samples).
#' @param mean_mat A matrix of mean expression data (genes x samples).
#' @param phase_pb A numeric vector of pseudobulk acrophases (in radians) for each gene, typically obtained from `is_cycling()`.
#' @param time_units A string indicating the time units ("hrs" for hours or "radians" for radians). Default is "hrs".
#' @return A tibble summarizing double-frequency rhythmicity metrics for each gene, including acrophase, amplitude, p-values, and adjusted p-values.
#' @noRd
var_cycling_fix_acro <- function(times, var_mat, mean_mat, phase_pb, time_units = "hrs") {
  `%do%` <- foreach::`%do%`
  cat(paste("\nRunning var_cycling_fix_acro()"))
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
    phase_bulk = as.numeric(phase_pb[gene_i])
    Gene_Symbols = colnames(gene1)[gene_i]

    # Check for all zeros or zero variance
    if (all(m1 == 0) || all(gexp1 == 0) || sd(m1) == 0 || sd(gexp1) == 0) {
      # Return NA row for this gene
      return(c(Gene_Symbols, rep(NA, 6))) # 9 columns to match the number of outputs
    }


    partial_model = lm(gexp1 ~ m1)
    partial_model_AIC = AIC(partial_model)
    partial_model_BIC = BIC(partial_model)

    double_wave = cos(2 * (times1- pi/2 - phase_bulk))
    full_model = lm(gexp1 ~ double_wave + m1)
    full_model_AIC = AIC(full_model)
    full_model_BIC = BIC(full_model)

    design_matrix <- model.matrix(gexp1 ~ double_wave + m1)
    
    #What is the amplitude of the fundamental according to V2?
    find_Af = lm(full_model$coefficients[["m1"]]*m1 ~ cos(times1-phase_bulk))
    fund_amp = find_Af$coefficients[["cos(times1 - phase_bulk)"]]

    anova_results = anova(partial_model, full_model)
    double_freq_p_statistic = anova_results$`Pr(>F)`[2]
    

    #Amplitude of the double frequency component
    double_freq_amplitude = full_model$coefficients[["double_wave"]]
    
    #MESOR or mean of the fitted variance wave
    # var_mesor = mean( full_model$fitted.values)
    var_mesor = mean(design_matrix[,-2] %*% full_model[["coefficients"]][-2])

    gene_summary = cbind( Gene_Symbols, double_freq_amplitude, 
            double_freq_p_statistic, fund_amp, var_mesor,
            partial_model_AIC, full_model_AIC,
            partial_model_BIC, full_model_BIC)
    
    
    
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
#' @param line_thickness A numeric value specifying the thickness of the lines in the plots. Default is 1.
#' @param point_size A numeric value specifying the size of the points in the plots. Default is 3.
#' @param line_color A string specifying the color of the lines in the plots. Default is "red".
#' @param point_color A string specifying the color of the points in the plots. Default is "blue".
#' @return Generates side-by-side plots for pseudobulk expression and double-frequency variance for the selected genes.
#' @export
plot_rhythms <- function(times, var_mat, mean_mat, plot_me = c(), time_units = "hrs", line_thickness = 1, point_size = 3, line_color = "red", point_color = "blue") {
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
  times1_hrs = times1 * 12 / pi
  
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
    pb_acro_hrs = pb_acro * 12 / pi
    pb_amp = sqrt(pb_model[["coefficients"]][["sin(times1)"]]^2 + pb_model[["coefficients"]][["cos(times1)"]]^2)
    pb_intercept = pb_model[['coefficients']]["(Intercept)"]

    var_partial <- lm(gexp1 ~ m1)
    var_model <- lm(gexp1 ~ sin(2 * times1) + cos(2 * times1) + m1)
    design_matrix <- model.matrix(gexp1 ~ sin(2 * times1) + cos(2 * times1) + m1)
    var_anova_results = anova(var_partial, var_model)
    sin_coeff = var_model[["coefficients"]][["sin(2 * times1)"]]
    cos_coeff = var_model[["coefficients"]][["cos(2 * times1)"]]
    double_freq_acrophase = atan2(sin_coeff, cos_coeff) %% (2*pi) / 2
    double_freq_acrophase_hrs = double_freq_acrophase * 12 / pi
    var_double_freq_amp = sqrt(sin_coeff^2 + cos_coeff^2)

    my_theme <- theme_bw()+
      theme(
        axis.text.x = element_text(size = 14),    # Increase x-axis text size
        axis.text.y = element_text(size = 14),    # Increase y-axis text size
        axis.title.x = element_text(size = 16),   # Increase x-axis title size
        axis.title.y = element_text(size = 16),   # Increase y-axis title size
        legend.title = element_text(size = 12),    # Increase legend title size
        strip.text = element_text(size = 12)
      )

   # Data for plot 1 
    pb_data <- data.frame(times1_hrs = times1_hrs, m1 = m1) 
    times_ideal <- seq(times1_hrs[1], times1_hrs[length(times1_hrs)], by = 0.1) 
    pb_ideal_data <- data.frame(times_ideal = times_ideal, 
                            pb_ideal = pb_amp * cos((times_ideal - pb_acro_hrs)*(pi/12)) + pb_intercept) 

    # Create the ggplot object 
    p1 <- ggplot(pb_data, aes(x = times1_hrs, y = m1)) + 
    geom_point(color = point_color, size = point_size) + 
    geom_line(data = pb_ideal_data, aes(x = times_ideal, y = pb_ideal), color = line_color, linewidth = line_thickness) + 
    labs(title = paste("PB Expression for", Gene_Symbols), 
          subtitle = paste("p-value:", round(pb_anova_results$`Pr(>F)`[2], 3), 
                          ", acrophase:", round(pb_acro, 3)), 
          x = "Time (hrs)", 
          y = "Bulk expression")+
    my_theme

    # Data for plot 2 
    circ <- gexp1 - design_matrix[, -c(2, 3)] %*% var_model[["coefficients"]][-c(2, 3)] 
    var_data <- data.frame(times1_hrs = times1_hrs, circ = circ) 
    var_ideal_data <- data.frame(times_ideal = times_ideal, 
                              var_ideal = var_double_freq_amp * cos(2 *(pi/12)*(times_ideal - double_freq_acrophase_hrs))) 

    # Create the ggplot object 
    p2 <- ggplot(var_data, aes(x = times1_hrs, y = circ)) + 
    geom_point(color = point_color, size = point_size) + 
    geom_line(data = var_ideal_data, aes(x = times_ideal, y = var_ideal), color = line_color, linewidth = line_thickness) + 
    labs(title = paste("Double Frequency Variance for", Gene_Symbols), 
          subtitle = paste("p-value:", round(var_anova_results$`Pr(>F)`[2], 3), 
                          ", acrophase:", round(double_freq_acrophase, 3)), 
          x = "Time (hrs)", 
          y = "12-hr component of var.")+
    my_theme
    
    show(p1+p2)
      }
}


#' Run cycling Analysis with MetaCycle
#'
#' This function utilizes the `meta2d` function from the MetaCycle package to 
#' identify rhythms in the provided dataset using the JTK_CYCLE and Lomb-Scargle methods.
#'
#' @param df A data frame where rows represent genes and columns represent time points.
#'           The row names should correspond to gene identifiers.
#' @param timepoints A numeric vector specifying the time points of the measurements.
#' @param out_file A string specifying the file path to save the input data for MetaCycle.
#' @param out_dir A string specifying the directory path to save the output files generated by MetaCycle.
#' @param minper A numeric value specifying the minimum period for the analysis (default is 24).
#' @param maxper A numeric value specifying the maximum period for the analysis (default is 24).
#' @param ... Additional arguments to be passed to the `meta2d` function.
#'
#' @details The function first prepares the input data by adding a "Genes" column 
#'          containing the row names of the input data frame. It then writes the 
#'          prepared data to a text file specified by `out_file`. Finally, it calls 
#'          the `meta2d` function from the MetaCycle package to perform the JTK_CYCLE and Lomb-Scargle
#'          analysis, with a default period of 24 hrs.
#'
#' @return This function does not return a value. It writes the results of the 
#'         analysis to the directory specified by `out_dir`.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' df <- data.frame(
#'   time1 = c(1, 2, 3),
#'   time2 = c(4, 5, 6),
#'   time3 = c(7, 8, 9)
#' )
#' rownames(df) <- c("Gene1", "Gene2", "Gene3")
#' run_meta2d(df, timepoints = c(0, 4, 8), out_file = "input.txt", out_dir = "output/")
#' }
#' @noRd
run_meta2d <- function(df,
                      timepoints,
                      out_file,
                      out_dir,
                      minper = 24,
                      maxper = 24,
                      ...){
  cat("\nIdentifying rhythms using MetaCycle\n")

  # Write data files to text files for MetaCycle
  df = cbind("Genes"=rownames(df), df)
  utils::write.table(df, file=out_file,sep="\t", quote=FALSE, row.names=FALSE)

  MetaCycle::meta2d(infile=out_file,
                    filestyle="txt",
                    outdir=out_dir,
                    timepoints=timepoints,
                    minper = minper, maxper = maxper,
                    cycMethod = c("JTK", "LS"),
                    ...)
}

#' Analyze double-frequency and fundamental frequency rhythms using meta2D
#'
#' This function calculates the double-frequency component of variance, mean rhythm, and pseudobulk rhythm
#' for each gene, using `run_meta2d()` for rhythm analysis, and reads the MetaCycle results into dataframes.
#'
#' @param times A numeric vector of time points (e.g., ZT times in hours).
#' @param var_mat A matrix of variance data (genes x samples).
#' @param mean_mat A matrix of mean expression data (genes x samples).
#' @param out_file_circ A string specifying the file path to save the input data for MetaCycle for 12 hr rhythms.
#' @param out_file_mean A string specifying the file path to save the input data for MetaCycle for 24 hr rhythms.
#' @param out_file_pb A string specifying the file path to save the input data for MetaCycle for pseudobulk rhythms.
#' @param out_dir A string specifying the directory path to save the output files generated by MetaCycle.
#' @param ... Additional arguments to be passed to the `run_meta2d()` function.
#' @return A list containing three dataframes:
#' \itemize{
#'   \item `double_freq_results`: meta2d results for double-frequency rhythms, including `var_mesor`.
#'   \item `fund_results`: meta2d results for mean rhythms.
#'   \item `pb_results`: meta2d results for pseudobulk rhythms.
#' }
#' @noRd
cycling_meta2d <- function(times, var_mat, mean_mat, 
                            out_file_circ, out_file_mean, out_file_pb, out_dir,  ...) {
   cat(paste("\nRunning metacycle2d, input times must be regularly spaced, e.g. 0-24 and in hrs units\n"))
   if(length(unique(diff(times))) != 1){
     stop("Error: times must be regularly spaced in hrs units")
   }
  `%do%` <- foreach::`%do%`
  
  # Transpose matrices to have subjects x genes
  gene1 <- t(var_mat)
  means1 <- t(mean_mat)
  times1 <- times
  # Initialize matrices to store the double freq signals, mean rhythms, and pseudobulk rhythms
  double_freq_mat <- matrix(NA, nrow = ncol(gene1), ncol = nrow(gene1))
  mean_rhythm_matrix <- matrix(NA, nrow = ncol(gene1), ncol = nrow(gene1))
  pb_mat <- matrix(NA, nrow = ncol(gene1), ncol = nrow(gene1))
  var_mesor_values <- numeric(ncol(gene1))  # Store var_mesor for each gene
  rownames(double_freq_mat) <- colnames(gene1)  # Gene names
  colnames(double_freq_mat) <- rownames(gene1)  # Subject names
  rownames(mean_rhythm_matrix) <- colnames(gene1)  # Gene names
  colnames(mean_rhythm_matrix) <- rownames(gene1)  # Subject names
  rownames(pb_mat) <- colnames(gene1)  # Gene names
  colnames(pb_mat) <- rownames(gene1)  # Subject names
  
  # Loop through genes to calculate the double-frequency component, mean rhythm, and pseudobulk rhythm
  foreach(gene_i = 1:ncol(gene1), .combine = rbind) %do% {
    gexp1 <- as.numeric(unlist(gene1[, gene_i]))
    m1 <- as.numeric(unlist(means1[, gene_i]))

    # Fit the variance model
    var_model <- lm(gexp1 ~ sin(2 * (pi/12)*times1) + cos(2 * (pi/12)*times1) + m1)
    design_matrix <- model.matrix(gexp1 ~ sin(2 * (pi/12)*times1) + cos(2 * (pi/12)*times1) + m1)
    var_mesor <- mean(var_model$fitted.values)
    var_mesor_values[gene_i] <- var_mesor  # Store var_mesor for this gene

    # Calculate the double frequency signal (double-frequency component)
    circ <- gexp1 - design_matrix[, -c(2, 3)] %*% var_model[["coefficients"]][-c(2, 3)]
    double_freq_mat[gene_i, ] <- circ
    
    # Calculate the mean rhythm
    mean_rhythm <- var_model$coefficients[["m1"]] * m1
    mean_rhythm_matrix[gene_i, ] <- mean_rhythm

    # Calculate the pseudobulk rhythm (just m1)
    pb_mat[gene_i, ] <- m1
  }
  
  # Run meta2d analysis for double frequency rhythms
  run_meta2d(df = double_freq_mat, timepoints = times1, out_file = out_file_circ, out_dir = out_dir, minper = 12, maxper = 12,...)
  
  # Run meta2d analysis for mean rhythm
  run_meta2d(df = mean_rhythm_matrix, timepoints = times1, out_file = out_file_mean, out_dir = out_dir, minper = 24, maxper = 24, ...)
  
  # Run meta2d analysis for pseudobulk rhythm
  run_meta2d(df = pb_mat, timepoints = times1, out_file = out_file_pb, out_dir = out_dir, minper = 24, maxper = 24, ...)
  
  # Read the meta2d results for circ
  meta2d_results_circ_file <- file.path(out_dir, paste0("meta2d_", tools::file_path_sans_ext(out_file_circ), ".txt"))
  if (file.exists(meta2d_results_circ_file)) {
    circ_results_raw <- read.table(meta2d_results_circ_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    expected_cols <- c("CycID", "JTK_adjphase", "JTK_amplitude", "LS_adjphase", "meta2d_phase", "meta2d_AMP", "meta2d_period", "meta2d_pvalue", "meta2d_BH.Q")
    missing_cols <- setdiff(expected_cols, colnames(circ_results_raw))
    if (length(missing_cols) > 0) {
      stop("Check that data is suitable for JTK and LS, meta2d results file for circ is missing expected columns: ", paste(missing_cols, collapse = ", "))
    }
    circ_results <- circ_results_raw %>%
      select(all_of(expected_cols))
    colnames(circ_results) <- c("Gene_Symbols",
                                "JTK_double_freq_acrophase", "JTK_double_freq_amplitude",
                                "LS_double_freq_acrophase",
                                "double_freq_acrophase", "double_freq_amplitude",
                                "double_freq_period", "double_freq_p_statistic", "double_freq_BH.Q")
    circ_results$var_mesor <- var_mesor_values  # Add var_mesor column
  } else {
    stop("meta2d results file for circ not found: ", meta2d_results_circ_file)
  }
  
  # Read the meta2d results for mean rhythm
  meta2d_results_mean_file <- file.path(out_dir, paste0("meta2d_", tools::file_path_sans_ext(out_file_mean), ".txt"))
  if (file.exists(meta2d_results_mean_file)) {
    mean_results_raw <- read.table(meta2d_results_mean_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    expected_cols <- c("CycID", "JTK_adjphase", "JTK_amplitude", "LS_adjphase", "meta2d_phase", "meta2d_AMP", "meta2d_period", "meta2d_pvalue", "meta2d_BH.Q")
    missing_cols <- setdiff(expected_cols, colnames(mean_results_raw))
    if (length(missing_cols) > 0) {
      stop("Check that data is suitable for JTK and LS, meta2d results file for mean rhythm is missing expected columns: ", paste(missing_cols, collapse = ", "))
    }
    mean_results <- mean_results_raw %>%
      select(all_of(expected_cols))
    colnames(mean_results) <- c("Gene_Symbols",
                                "JTK_fund_acrophase", "JTK_fund_amplitude",
                                "LS_fund_acrophase",
                                "fund_acrophase", "fund_amplitude",
                                "fund_period", "fund_p_statistic", "fund_BH.Q")
  } else {
    stop("meta2d results file for mean rhythm not found: ", meta2d_results_mean_file)
  }

  # Read the meta2d results for pseudobulk rhythm
  meta2d_results_pb_file <- file.path(out_dir, paste0("meta2d_", tools::file_path_sans_ext(out_file_pb), ".txt"))
  if (file.exists(meta2d_results_pb_file)) {
    pb_results_raw <- read.table(meta2d_results_pb_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    expected_cols <- c("CycID", "JTK_adjphase", "JTK_amplitude", "LS_adjphase", "meta2d_phase", "meta2d_AMP", "meta2d_period", "meta2d_pvalue", "meta2d_BH.Q")
    missing_cols <- setdiff(expected_cols, colnames(pb_results_raw))
    if (length(missing_cols) > 0) {
      stop("Check that data is suitable for JTK and LS, meta2d results file for pseudobulk rhythm is missing expected columns: ", paste(missing_cols, collapse = ", "))
    }
    pb_results <- pb_results_raw %>%
      select(all_of(expected_cols))
    colnames(pb_results) <- c("Gene_Symbols",
                            "JTK_acrophase_bulk", "JTK_amplitude_bulk",
                            "LS_acrophase_bulk", 
                            "acrophase_bulk", "amplitude_bulk",
                            "period_bulk", "p_statistic", "BH.Q_bulk")
  } else {
    stop("meta2d results file for pseudobulk rhythm not found: ", meta2d_results_pb_file)
  }

  # Return all results as a list
  return(list(double_freq_results = circ_results, fund_results = mean_results, pb_results = pb_results))
}

utils::globalVariables(c(
  "gene_i", "all_of", "seurat_clusters", "nCount_RNA", "theme_bw", "theme", "element_text","pb_ideal", "var_ideal"
))
