#' Calculate Analytic Parameters
#'
#' This function calculates analytic parameters (G, Ac, and p) based on input amplitudes and mesor values.
#'
#' @param V_2w Double frequency variance amplitude.
#' @param V_w Fundamental frequency of variance amplitude.
#' @param S_amp pseudobulk cycling amplitude.
#' @param S_DC Mesor of the pseudobulk cycling.
#' @return A numeric vector containing G (Mesor), Ac (Cell rhythm amplitude), and p (phase shift_standard_dev).
#' @export

find_analytic = function(V_2w, V_w, S_amp, S_DC){
  
  # Solve for G: trivial
  G = S_DC
  
  #solve for Theta (NB dispersion)
  theta = ((V_w / S_amp) - 1) / (2 * G)

  R = V_2w / S_amp^2
  log_argument = (1 - 2*R) / (1 + theta)
  
  #solve for phi, sd of phase dispersion
  phi_squared = -log(log_argument)
  phi = sqrt(phi_squared)
  
  # Solve for Ac, amp of cell rhythm
  Ac = S_amp / exp(-phi_squared / 2) 
  
  return(c(G, Ac, theta, phi))
}

#' Calculate Phi for Rhythmic Genes
#'
#' This function calculates the phi parameter for rhythmic genes by combining results from pseudobulk cycling and double-frequency variance cycling.
#'
#' @param pb_cycling_results df of pseudobulk cycling results, as obtained from is_cycling().
#' @param var_cycling_results df of variance cycling results, as obtained from var_cycling().
#' @param cycling_p_thr A numeric threshold for cycling p-values to filter significant genes. Default is 0.05.
#' @param bulk_amp_thr A numeric threshold for the bulk amplitude. Default is 0.5.
#' @return A data frame summarizing the phi parameters (G, Ac, theta, and phi) for rhythmic genes, or an empty data frame if no genes meet the criteria.
#' @importFrom dplyr filter select mutate across left_join %>% everything
#' @importFrom tibble column_to_rownames
#' @importFrom stats lm anova p.adjust model.matrix sd var
#' @importFrom graphics points
#' @export
calc_phi = function(pb_cycling_results, var_cycling_results, cycling_p_thr = 0.05,  bulk_amp_thr = 0.5) {
    # Check if pb_cycling_results or var_cycling_results are empty
    if (nrow(pb_cycling_results) == 0) {
        stop("Error: pb_cycling_results is empty. No genes to process.")
    }
    if (nrow(var_cycling_results) == 0) {
        stop("Error: var_cycling_results is empty. No genes to process.")
    }

    # Get significant cycling genes
    circ_genes = filter(pb_cycling_results, as.numeric(p_statistic) < cycling_p_thr) %>% 
        select(Gene_Symbols) %>% 
        unname %>% 
        unlist
    print(paste("Number of cycling genes: ", length(circ_genes)))

    input = var_cycling_results %>% 
        left_join(pb_cycling_results, by = "Gene_Symbols") %>%
        mutate(across(-Gene_Symbols, as.numeric)) %>%
        filter(p_statistic < cycling_p_thr) %>%
        # mutate(acrodiff = (double_freq_acrophase - acrophase_bulk) %% pi) %>%
        # mutate(keep_acrodiff = acrodiff < (pi / 2 + expected_acrodiff_thr) & acrodiff > (pi / 2 - expected_acrodiff_thr)) %>%
        dplyr::select(Gene_Symbols, double_freq_amplitude,fund_amp,amplitude_bulk, mesor_bulk, partial_model_AIC, full_model_AIC,  everything()) %>%
        column_to_rownames("Gene_Symbols")%>%
        filter(partial_model_AIC > full_model_AIC)

    # hist(input$acrodiff, main = "Acrophase difference histogram", xlab = "Acrophase difference (radians)")
    # Check if input is empty after filtering
    if (nrow(input) == 0) {
        warning("Warning: No genes meet the filtering criteria. Returning an empty data frame.")
        return(data.frame(G = numeric(0), Ac = numeric(0), phi = numeric(0)))
    }

    # print(paste("Number of genes with acrophase difference within expected range:", nrow(filter(input, keep_acrodiff == TRUE))))

    # input = input %>%
    #     filter(keep_acrodiff == TRUE) %>%
    #     dplyr::select(Gene_Symbols, double_freq_amplitude, var_mesor, fund_amp, everything()) %>%
    #     column_to_rownames("Gene_Symbols")

    # Check if input is empty after filtering for PB amplitude
    if (nrow(filter(input, amplitude_bulk >= bulk_amp_thr)) == 0) {
        warning("Warning: No genes have PB amplitude greater than the threshold. Returning an empty data frame.")
        return(data.frame(G = numeric(0), Ac = numeric(0), phi = numeric(0)))
    }

    print(paste("Number of genes with PB amplitude >", bulk_amp_thr, ":", nrow(filter(input, amplitude_bulk >= bulk_amp_thr))))

    # Calculate phi parameters
    phi = input %>%
        filter(amplitude_bulk >= bulk_amp_thr) %>%
        apply(., 1, function(row) {
            Orpheus::find_analytic(row[1], row[2], row[3], row[4])
        }) %>%
        t %>%
        as.data.frame

    colnames(phi) = c("G", "Ac", "theta", "phi")
    return(phi)
}