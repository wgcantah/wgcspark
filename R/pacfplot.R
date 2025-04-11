#' Plot Partial Autocorrelation Function (PACF) with Stata Theme and Critical Bounds
#'
#' This function calculates and plots the partial autocorrelation function (PACF) for one or more numeric time
#' series variables contained in a data frame. The PACF values are displayed as a bar graph, and the function
#' overlays 5% critical bounds (calculated as ±1.96/√n). The Stata theme from ggthemes is applied.
#'
#' @param data A data frame containing the time series data.
#' @param vars A character vector of variable names to analyze.
#' @param lag.max Maximum lag at which to calculate the PACF. Default is 20.
#' @param title Title for the plot. Default is "Partial Autocorrelation Function".
#'
#' @return A ggplot object.
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   var1 = cumsum(stats::rnorm(100)),
#'   var2 = cumsum(stats::rnorm(100))
#' )
#' pacf_plot <- pacfplot(df, vars = c("var1", "var2"), lag.max = 20,
#'                       title = "PACF for Example Variables")
#' print(pacf_plot)
#'
#' @export
pacfplot <- function(data, vars, lag.max = 20, title = "Partial Autocorrelation Function") {
  # Check that each variable exists and is numeric.
  for (v in vars) {
    if (!v %in% names(data))
      stop(paste("Variable", v, "not found in the data frame."))
    if (!is.numeric(data[[v]]))
      stop(paste("Variable", v, "must be numeric."))
  }

  # Compute PACF for each variable and calculate the 5% critical bound.
  pacf_list <- lapply(vars, function(v) {
    x <- data[[v]]
    n <- length(x)
    crit <- 1.96 / sqrt(n)
    pacf_obj <- stats::pacf(x, plot = FALSE, lag.max = lag.max)
    # The pacf function does not provide explicit lag values; create them manually.
    lags <- seq_along(pacf_obj$acf) - 1
    pacf_values <- pacf_obj$acf
    df <- data.frame(variable = v, lag = lags, pacf = pacf_values)
    attr(df, "crit") <- crit
    return(df)
  })

  pacf_df <- do.call(rbind, pacf_list)

  # Create a data frame of critical bounds for each variable.
  crit_df <- do.call(rbind, lapply(pacf_list, function(df) {
    data.frame(variable = unique(df$variable), crit = attr(df, "crit"))
  }))

  # Build the bar graph.
  p <- ggplot2::ggplot(pacf_df, ggplot2::aes(x = lag, y = pacf)) +
    ggplot2::geom_col(fill = "navy", width = 0.1) +  # Bars for PACF values.
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
    ggplot2::geom_hline(data = crit_df, mapping = ggplot2::aes(yintercept = crit),
                        linetype = "dashed", color = "red", size = 1) +
    ggplot2::geom_hline(data = crit_df, mapping = ggplot2::aes(yintercept = -crit),
                        linetype = "dashed", color = "red", size = 1) +
    ggplot2::labs(title = title, x = "Lag", y = "PACF") +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggthemes::theme_stata(base_size = 12)

  return(p)
}
