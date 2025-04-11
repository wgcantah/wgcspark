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

