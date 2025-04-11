#' Plot Autocorrelation Function (ACF) with Stata Theme and Critical Bounds
#'
#' This function calculates and plots the autocorrelation function (ACF) for one or more numeric time
#' series variables contained in a data frame. The ACF values are displayed as a bar graph, and the function
#' overlays 5% critical bounds (calculated as ±1.96/√n). The Stata theme from ggthemes is applied.
#'
#' @param data A data frame containing the time series data.
#' @param vars A character vector of variable names to analyze.
#' @param lag.max Maximum lag at which to calculate the ACF. Default is 20.
#' @param title Title for the plot. Default is "Autocorrelation Function".
#'
#' @return A ggplot object.
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   var1 = cumsum(rnorm(100)),
#'   var2 = cumsum(rnorm(100))
#' )
#' acf_plot <- acfplot(df, vars = c("var1", "var2"), lag.max = 20,
#'                            title = "ACF for Example Variables")
#' print(acf_plot)
#'
#' @import ggplot2 ggthemes
#' @export
acfplot <- function(data, vars, lag.max = 20, title = "Autocorrelation Function") {
  # Check if each variable exists and is numeric.
  for (v in vars) {
    if (!v %in% names(data))
      stop(paste("Variable", v, "not found in the data frame."))
    if (!is.numeric(data[[v]]))
      stop(paste("Variable", v, "must be numeric."))
  }

  # Compute ACF for each variable and calculate the 5% critical bound.
  acf_list <- lapply(vars, function(v) {
    x <- data[[v]]
    n <- length(x)
    crit <- 1.96 / sqrt(n)
    acf_obj <- stats::acf(x, plot = FALSE, lag.max = lag.max)
    lags <- as.vector(acf_obj$lag)
    acf_values <- as.vector(acf_obj$acf)
    df <- data.frame(variable = v, lag = lags, acf = acf_values)
    attr(df, "crit") <- crit
    return(df)
  })

  acf_df <- do.call(rbind, acf_list)

  # Create a data frame of critical bounds for each variable.
  crit_df <- do.call(rbind, lapply(acf_list, function(df) {
    data.frame(variable = unique(df$variable), crit = attr(df, "crit"))
  }))

  # Build the bar graph.
  p <- ggplot2::ggplot(acf_df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_col(fill = "navy", width = 0.1) +  # Bars for ACF values.
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
    ggplot2::geom_hline(data = crit_df, mapping = ggplot2::aes(yintercept = crit),
                        linetype = "dashed", color = "red", size = 1) +
    ggplot2::geom_hline(data = crit_df, mapping = ggplot2::aes(yintercept = -crit),
                        linetype = "dashed", color = "red", size = 1) +
    ggplot2::labs(title = title, x = "Lag", y = "ACF") +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggthemes::theme_stata(base_size = 12)

  return(p)
}
