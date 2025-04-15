#' Run Unit Root Tests on Multiple Time Series
#'
#' This function runs several unit root tests (ADF, KPSS, Phillips-Perron, and NG-Perron)
#' on specified time series variables contained in a data frame. It extracts the test statistics,
#' p-values, and the critical values (1%, 5%, and 10%) into separate columns.
#' The results are printed to the console and exported to a Word document using \code{officer::read_docx},
#' \code{officer::body_add_par}, \code{flextable::flextable}, and \code{officer::print}.
#'
#' @param data A data frame containing at least one time series (of class \code{ts}).
#' @param varnames A character vector specifying the column names in \code{data} that are time series.
#'
#' @return Invisibly returns a combined data frame that contains the results for each variable.
#'
#' @examples
#' \dontrun{
#' # Simulate three monthly time series (10 years; 120 observations each)
#' set.seed(123)
#' ts1 <- ts(rnorm(120, mean = 0, sd = 1), start = c(2000, 1), frequency = 12)
#' set.seed(456)
#' ts2 <- ts(rnorm(120, mean = 5, sd = 2), start = c(2000, 1), frequency = 12)
#' set.seed(789)
#' ts3 <- ts(rnorm(120, mean = -3, sd = 1.5), start = c(2000, 1), frequency = 12)
#'
#' # Create a data frame with the simulated time series.
#' monthly <- data.frame(
#'   CPI = ts1,
#'   EXC_END = ts2,
#'   DEPRECIATION = ts3,
#'   stringsAsFactors = FALSE
#' )
#'
#' # Run the unit root tests on the selected variables.
#' my_results <- urtest(my_ts_df, c("CPI", "EXRATE", "REVENUE")) 
#'
#' # Print the combined results.
#' print(results_multi)
#' }
#'
#' @export
urtest <- function(data, varnames) {
  #--------------------------- 1.  Input checks ---------------------------#
  if (!is.character(varnames))
    stop("'varnames' must be a character vector of column names.")

  missing_vars <- setdiff(varnames, names(data))
  if (length(missing_vars))
    stop("These columns are not found in data: ",
         paste(missing_vars, collapse = ", "))

  # Required packages (only those truly needed now)
  required_pkgs <- c("tseries", "urca")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required but not installed.", pkg))
  }

  #--------------------------- 2.  Loop over vars -------------------------#
  results_list <- list()

  for (varname in varnames) {

    ts_data <- data[[varname]]

    if (!inherits(ts_data, "ts")) {
      warning(sprintf("'%s' is not a ts object. Skipping.", varname))
      next
    }

    ##--- 2.1  ADF --------------------------------------------------------##
    adf_ur      <- urca::ur.df(ts_data, type = "drift", selectlags = "AIC")
    adf_stat    <- adf_ur@teststat["tau2"]
    adf_cvals   <- adf_ur@cval["tau2", ]
    adf_pval    <- tseries::adf.test(ts_data)$p.value

    ##--- 2.2  KPSS -------------------------------------------------------##
    kpss_test   <- tseries::kpss.test(ts_data)
    kpss_stat   <- kpss_test$statistic
    kpss_pval   <- kpss_test$p.value
    kpss_cvals  <- c("1pct" = 0.739, "5pct" = 0.463, "10pct" = 0.347)

    ##--- 2.3  Phillips–Perron -------------------------------------------##
    pp_test     <- urca::ur.pp(ts_data, type = "Z-tau",
                               model = "constant", lags = "short")
    pp_stat     <- pp_test@teststat
    pp_cvals    <- pp_test@cval
    pp_pval     <- NA   # not returned by ur.pp

    ##--- 2.4  NG–Perron --------------------------------------------------##
    ng_test     <- urca::ur.ers(ts_data, model = "trend", lag.max = 4)
    ng_stat     <- ng_test@teststat[1]
    ng_cvals    <- ng_test@cval[1, ]
    ng_pval     <- NA   # not returned by ur.ers

    ##--- 2.5  Assemble ---------------------------------------------------##
    results_df <- data.frame(
      Variable   = varname,
      Test       = c("ADF", "KPSS", "Phillips–Perron", "NG–Perron"),
      Statistic  = c(adf_stat, kpss_stat, pp_stat, ng_stat),
      P_Value    = c(adf_pval, kpss_pval, pp_pval, ng_pval),
      Crit_1pct  = c(adf_cvals["1pct"], kpss_cvals["1pct"],
                     pp_cvals["1pct"],  ng_cvals["1pct"]),
      Crit_5pct  = c(adf_cvals["5pct"], kpss_cvals["5pct"],
                     pp_cvals["5pct"],  ng_cvals["5pct"]),
      Crit_10pct = c(adf_cvals["10pct"], kpss_cvals["10pct"],
                     pp_cvals["10pct"], ng_cvals["10pct"]),
      stringsAsFactors = FALSE
    )

    results_list[[varname]] <- results_df
  }

  #--------------------------- 3.  Return ---------------------------------#
  combined_results <- do.call(rbind, results_list)

  # Show in console for immediate feedback (optional)
  print(combined_results)

  invisible(combined_results)
}



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
