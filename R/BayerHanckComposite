#' Bayer–Hanck Composite Cointegration Test Function
#'
#' This function performs the Bayer–Hanck composite cointegration test by combining the p-values from:
#'   1. Engle–Granger test (via a unit root test on residuals using tseries::adf.test)
#'   2. Phillips–Ouliaris test (via aTSA::po.test)
#'   3. Optionally, the Johansen cointegration test (via urca::ca.jo)
#'
#' The p-values are combined using Fisher's method. A decision is made at the specified significance level.
#'
#' @param data A data frame containing the time series.
#' @param vars A character vector of column names. The first variable is taken as the dependent series.
#' @param lags Number of lags to use in the Johansen test if included (default is 2).
#' @param include_Johansen Logical indicating whether to include the Johansen test (default is FALSE).
#' @param alpha Significance level for the hypothesis test (default is 0.05).
#'
#' @return A list with the individual test p-values, Fisher test statistic, degrees of freedom,
#' the combined p-value, and the decision on the null hypothesis.
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   n <- 100
#'   df <- data.frame(
#'     y = cumsum(rnorm(n)),
#'     x1 = cumsum(rnorm(n)),
#'     x2 = cumsum(rnorm(n))
#'   )
#'   result <- BayerHanckComposite(data = df, vars = c("y", "x1", "x2"),
#'                                 lags = 2, include_Johansen = TRUE, alpha = 0.05)
#'   print(result)
#' }
BayerHanckComposite <- function(data, vars, lags = 2, include_Johansen = FALSE, alpha = 0.05) {
  
  # Check that all specified variables exist in the data frame
  if (!all(vars %in% names(data))) {
    stop("Not all specified variables are found in the data frame!")
  }
  
  # Subset the data to only the selected variables
  selected_data <- data[, vars]
  
  # Ensure that at least two variables are provided (dependent and independent variables)
  if (ncol(selected_data) < 2) {
    stop("Please select at least two variables (one dependent and at least one independent variable).")
  }
  
  # Assume the first variable is the dependent variable; the rest are independent variables.
  dep <- selected_data[, 1]
  indep <- selected_data[, -1, drop = FALSE]
  
  # -------------------------------
  # 1. Engle–Granger (EG) Test
  # -------------------------------
  # Run the cointegration regression; note that we use the dot notation for independent variables
  eg_model <- stats::lm(dep ~ ., data = as.data.frame(indep))
  eg_resids <- stats::resid(eg_model)
  
  # Use tseries::adf.test to test for a unit root in the residuals
  eg_test <- tseries::adf.test(eg_resids)
  p_EG <- eg_test$p.value
  cat("Engle–Granger test p-value:", p_EG, "\n")
  
  # -------------------------------
  # 2. Phillips–Ouliaris (PO) Test
  # -------------------------------
  # Use aTSA::po.test on the matrix of selected variables.
  po_test <- aTSA::po.test(as.matrix(selected_data))
  p_PO <- po_test$p.value
  cat("Phillips–Ouliaris test p-value:", p_PO, "\n")
  
  # Initialize the vector of p-values with the EG and PO test results.
  p_values <- c(Engle_Granger = p_EG, Phillips_Ouliaris = p_PO)
  
  # -------------------------------
  # 3. Optional: Johansen Cointegration Test
  # -------------------------------
  if (include_Johansen) {
    johansen_res <- urca::ca.jo(as.matrix(selected_data), type = "trace", ecdet = "const", K = lags)
    # Extract the trace statistic for r = 0 (i.e., no cointegration)
    stat <- johansen_res@teststat[1]
    # Use the 5% critical value from the test as a rough approximation for a p-value.
    crit <- johansen_res@cval[1, "5pct"]
    if (stat > crit) {
      p_J <- 0.05
    } else {
      p_J <- 1
    }
    cat("Johansen test (approximate) p-value:", p_J, "\n")
    p_values <- c(p_values, Johansen = p_J)
  }
  
  # -------------------------------
  # 4. Combine p-values using Fisher's Method
  # -------------------------------
  fisher_stat <- -2 * sum(log(p_values))
  df <- 2 * length(p_values)
  combined_p <- 1 - stats::pchisq(fisher_stat, df = df)
  
  cat("\nFisher test statistic:", fisher_stat, "\n")
  cat("Degrees of freedom:", df, "\n")
  cat("Combined p-value:", combined_p, "\n")
  
  # -------------------------------
  # 5. Decision on the Null Hypothesis
  # -------------------------------
  decision <- ifelse(combined_p < alpha, 
                     "Reject the null hypothesis of no cointegration.", 
                     "Fail to reject the null hypothesis of no cointegration.")
  cat("\nDecision (at alpha =", alpha, "):", decision, "\n")
  
  # Return the results as a list.
  return(list(
    individual_pvalues = p_values,
    fisher_statistic = fisher_stat,
    degrees_of_freedom = df,
    combined_pvalue = combined_p,
    decision = decision
  ))
}
