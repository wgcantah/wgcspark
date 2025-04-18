% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsplot.R
\name{tsplot}
\alias{tsplot}
\title{Create a Combined Time Series Plot}
\usage{
tsplot(data, vars, freq = "monthly", start_date = "2000-01-01")
}
\arguments{
\item{data}{A data.frame containing your time series data (without a date column).}

\item{vars}{A character vector of time series variable names to plot.}

\item{freq}{A character string specifying the frequency of the time series. It can be "monthly", "quarterly", or "annual". Default is "monthly".}

\item{start_date}{A string representing the first date (e.g. "2000-01-01"). Default is "2000-01-01".}
}
\value{
A ggplot object with all time series plotted on a single panel.
}
\description{
Generates a single-pane time series plot with multiple variables using a colorful theme.
}
\examples{
set.seed(123)
mydata <- data.frame(
  bioRxiv         = cumsum(rnorm(24, mean = 10, sd = 5)),
  arXiv_q_bio     = cumsum(rnorm(24, mean = 8,  sd = 3)),
  PeerJ_Preprints = cumsum(rnorm(24, mean = 6,  sd = 2))
)

ts_combined_plot <- tsplot(
  data       = mydata,
  vars       = c("bioRxiv", "arXiv_q_bio", "PeerJ_Preprints"),
  freq       = "monthly",
  start_date = "2014-01-01"
)

print(ts_combined_plot)

}
