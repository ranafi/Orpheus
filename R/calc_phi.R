#' Calculate Analytic Parameters
#'
#' This function calculates analytic parameters (G, Ac, and p) based on input amplitudes and mesor values.
#'
#' @param Av Double frequency variance amplitude.
#' @param Mv Fundamental frequency of variance mesor.
#' @param Af Fundamental frequency of variance amplitude.
#' @return A numeric vector containing G (Mesor), Ac (Cell rhythm amplitude), and p (phase shift_standard_dev).
#' @export
find_analytic = function(Av, Mv, Af){
  G = (-Af^2*Av+Af^2*Mv-2*Av*Mv)/(Af^2-2*Av)
  Ac = (Af^2/sqrt((Af^2-2*Av)))
  p = (sqrt(2)* sqrt( log( Af/sqrt((Af^2-2*Av)) ) ))
  return(c(G, Ac, p))
}

#' Calculate Phi for Rhythmic Genes
#'
#' This function calculates the phi parameter for rhythmic genes by combining results from pseudobulk cycling and double-frequency variance cycling.
#'
#' @param ZT A numeric vector of time points (e.g., ZT times).
#' @param pb_cycling_results df of pseudobulk cycling results, as obtained from is_cycling().
#' @param var_cycling_results df of variance cycling results, as obtained from var_cycling().
#' @param cycling_p_thr A numeric threshold for cycling p-values to filter significant genes. Default is 0.05.
#' @param expected_acrodiff_thr A numeric threshold for the expected acrophase difference. Default is 0.524.
#' @param bulk_amp_thr A numeric threshold for the bulk amplitude. Default is 0.5.
#' @return A data frame summarizing the phi parameters (G, Ac, and phi) for rhythmic genes.
#' @importFrom dplyr filter select mutate across left_join %>% everything
#' @importFrom tibble column_to_rownames
#' @importFrom stats lm anova p.adjust model.matrix sd var
#' @importFrom graphics points
#' @export
calc_phi = function(ZT, pb_cycling_results, var_cycling_results, cycling_p_thr = 0.05, expected_acrodiff_thr = 0.524, bulk_amp_thr = 0.5){
    #Get significant cycling genes
    circ_genes = filter(pb_cycling_results, ) %>% select(Gene_Symbols) %>% unname %>% unlist
    print(paste("Number of cycling genes: ", length(circ_genes)))

    #merge cycling and var_cycling results
    input = var_cycling_results %>% 
    left_join(pb_cycling_results, by = "Gene_Symbols") %>%
    filter(as.numeric(p_statistic) < cycling_p_thr) %>%
    mutate(across(-Gene_Symbols, as.numeric)) %>%
    mutate( acrodiff=(double_freq_acrophase - acrophase_bulk)%% pi) %>%
    mutate(keep_acrodiff = acrodiff < (pi/2+expected_acrodiff_thr) & acrodiff > (pi/2-expected_acrodiff_thr) )
    
    print(paste("Number of genes with acrophase difference within expected range:", nrow(filter(input, keep_acrodiff == T))))

    input = input %>%
    filter( keep_acrodiff == T) %>%
    dplyr::select(Gene_Symbols, double_freq_amplitude, var_mesor, fund_amp, everything()) %>% 
    column_to_rownames("Gene_Symbols") 

    print(paste("Number of genes with PB amplitude >", bulk_amp_thr, ":", nrow(filter(input, amplitude_bulk > bulk_amp_thr))))


    phi = input %>% 
    filter(amplitude_bulk > bulk_amp_thr) %>%
    apply(., 1, function(row) {
        Orpheus::find_analytic(row[1], row[2], row[3] )
    }) %>% 
    t %>%
    as.data.frame

    colnames(phi) = c("G", "Ac", "phi")
    return(phi)
}