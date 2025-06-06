
# wgcspark

**wgcspark** is an R package that provides a comprehensive toolkit for
time series analysis and visualization. The package includes functions
to convert data frames to time series objects, create combined and
individual time series plots with attractive themes, plot
autocorrelation (ACF) and partial autocorrelation (PACF) as bar graphs
(with critical bounds), and run a suite of unit root tests.

## Features

- **Data Conversion**  
  Convert numeric columns in a data frame to time series objects with
  corresponding Date values using:
  - `tsconvert()`
- **Time Series Plotting**  
  Visualize your data with different plotting styles:
  - `tsplot()` creates a single combined plot.
  - `tsplot_ind()` generates individual plots for each series.
  - `tspanel_plot()` arranges time series plots in a grid layout.
- **ACF and PACF Visualization**  
  Plot autocorrelation and partial autocorrelation functions as bar
  graphs (with ±5% critical bounds) using:
  - `acfplot()`
  - `pacfplot()`
- **Unit Root Testing**  
  Run multiple unit root tests (ADF, KPSS, Phillips–Perron, NG–Perron)
  on your series and output results to a formatted Word document:
  - `urtest()`

## Installation

You can install the development version of **wgcspark** from GitHub
using **devtools**:

### Converting Data to Time Series

``` r
# install.packages("devtools")  # if not already installed
devtools::install_github("wgcantah/wgcspark")

library(wgcspark)

# Create a sample data frame
df <- data.frame(
  sales = c(100, 150, 200, 250),
  profit = c(50, 60, 70, 80),
  category = c("A", "B", "A", "C")
)

# Convert numeric columns to time series with monthly frequency starting in January 2020
monthly_ts <- tsconvert(df, start = c(2020, 1), period = "monthly")
str(monthly_ts)
```

### Creating Combined and Individual Time Series Plots

``` r
# Combined plot
ts_combined_plot <- tsplot(
  data = monthly_ts,
  vars = c("sales", "profit"),
  freq = "monthly",
  start_date = "2020-01-01"
)
print(ts_combined_plot)

# Individual plots
ts_ind_plots <- tsplot_ind(
  data = monthly_ts,
  vars = c("sales", "profit"),
  freq = "monthly",
  start_date = "2020-01-01"
)
lapply(ts_ind_plots, print)

# Panel plot
panel_plot <- tspanel_plot(
  data = df,
  vars = c("sales", "profit"),
  freq = "monthly",
  start_date = "2020-01-01"
)
print(panel_plot)
```

### Plotting ACF and PACF

``` r
# Plot ACF for the sales and profit series
acf_plot <- acfplot(
  data = monthly_ts,
  vars = c("sales", "profit"),
  lag.max = 20,
  title = "ACF for Sales and Profit"
)
print(acf_plot)

# Plot PACF for the same series
pacf_plot <- pacfplot(
  data = monthly_ts,
  vars = c("sales", "profit"),
  lag.max = 20,
  title = "PACF for Sales and Profit"
)
print(pacf_plot)
```

### Running Unit Root Tests

``` r
# Assume monthly_ts has time series objects in the columns 'sales' and 'profit'
results <- urtest(
  data = monthly_ts,
  varnames = c("sales", "profit")
)
print(results)
```

### Documentation

After installing the package, you can access help pages for each
function using the standard R help system. For example:

``` r
?tsconvert
?tsplot
?acfplot
?urtest
```

### License

This package is licensed under the GPL-3 license.
