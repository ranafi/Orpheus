#' Objective Function for Nonlinear Optimization
#'
#' This function defines the objective function to minimize during the optimization process.
#'
#' @param params A named vector of parameters to optimize, including `p`.
#' @param data A data frame containing the input data with columns `Af` (fundamental frequency amplitude) 
#' and `Av` (double frequency variance amplitude).
#' @return The sum of squared residuals between the predicted and observed values.
#' @keywords internal
objective <- function(params, data) {
  p <- params["p"]
  
  Af <- data$Af
  Av <- data$Av
  
  # Calculate the predicted p values for each observation
  pred_p <- sqrt(2) * sqrt(log(Af / sqrt(Af^2 - 2 * Av)))
  
  # We want all these to be close to the single p value
  residuals <- pred_p - p
  
  sum(residuals^2)  # Sum of squared residuals
}


#' Estimate Single Phi Parameters
#'
#' This function estimates the single phi parameter (`p`) for a given dataset using nonlinear least squares optimization.
#' It also calculates the `G` (mesor) and `Ac` (amplitude) values for each observation based on the optimized `p`.
#'
#' @param measurement_data A data frame containing the input data with columns `Af` (fundamental frequency amplitude), 
#' `Av` (double frequency variance amplitude), and `Mv` (variance mesor).
#' @param initial_guess A named vector providing the initial guess for the parameter `p`. Default is `c(p = 1.0)`.
#' @return A data frame with columns:
#' \itemize{
#'   \item `G`: The mesor for each observation.
#'   \item `Ac`: The cellular amplitude for each observation.
#'   \item `p`: The optimized single phi parameter.
#' }
#' @importFrom minpack.lm nls.lm
#' @noRd
estimate_single_phi <- function(measurement_data, initial_guess = c(p = 1.0)) {
  # Perform the optimization
  fit <- nls.lm(par = initial_guess, 
                fn = objective, 
                data = measurement_data)

  # Extract the optimal p value
  optimal_p <- coef(fit)["p"]

  # Calculate G and Ac for each observation
  results <- data.frame(
    G = (-measurement_data$Af^2 * measurement_data$Av + measurement_data$Af^2 * measurement_data$Mv - 
           2 * measurement_data$Av * measurement_data$Mv) / (measurement_data$Af^2 - 2 * measurement_data$Av),
    Ac = (measurement_data$Af^2 / sqrt((measurement_data$Af^2 - 2 * measurement_data$Av))),
    p = optimal_p
  )

  return(results)
}


#' Bootstrap Confidence Intervals for Phi Parameter
#'
#' This function estimates confidence intervals for the single phi parameter (`p`) using bootstrapping.
#'
#' @param measurement_data A data frame containing the input data with columns `Af` (fundamental frequency amplitude), 
#' `Av` (double frequency variance amplitude), and `Mv` (fundamental frequency mesor).
#' @param optimal_p phi, optimized by estimate_single_phi(). Suggested to use the output of estimate_single_phi()$p[1].
#' @param nPerms An integer specifying the number of bootstrap replicates. Default is 1000.

#' @return Prints the confidence intervals for the parameter `p` using different methods:
#' \itemize{
#'   \item Normal approximation CI
#'   \item Percentile CI
#'   \item Bias-corrected and accelerated (BCa) CI
#' }
#' @importFrom boot boot boot.ci
#' @noRd
estimate_CI <- function(measurement_data, optimal_p, nPerms = 1000) {

boot_fn <- function(data, indices) {
  d <- data[indices, ]
  fit <- nls.lm(
    par = list(p = optimal_p),  # Use original estimate as starting point
    fn = function(p, data) {
      pred_p <- sqrt(2) * sqrt(log(data$Af / sqrt(data$Af^2 - 2 * data$Av)))
      pred_p - p$p
    },
    data = d
  )
  coef(fit)
}
  # Run bootstrap with error handling
  boot_results <- boot(
    data = measurement_data,
    statistic = boot_fn,
    R = nPerms  # Number of bootstrap replicates
  )

  # Examine bootstrap distribution
  plot(boot_results, index = 1)  # Plot for parameter p

  # Calculate several types of CIs
  cat("Normal approximation CI:\n")
  print(boot.ci(boot_results, type = "norm", index = 1))

  cat("\nPercentile CI:\n")
  print(boot.ci(boot_results, type = "perc", index = 1))

  cat("\nBCa (bias-corrected) CI:\n")
  print(boot.ci(boot_results, type = "bca", index = 1))
}



