% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/histogram_stats.R
\name{histogram_stats}
\alias{histogram_stats}
\alias{hist_mean}
\alias{hist_var}
\alias{hist_cov}
\alias{hist_cor}
\title{Statistics for Histogram Data}
\usage{
hist_mean(x, var_name, method = "BG", ...)

hist_var(x, var_name, method = "BG", ...)

hist_cov(x, var_name1, var_name2, method = "BG")

hist_cor(x, var_name1, var_name2, method = "BG")
}
\arguments{
\item{x}{histogram-valued data object.}

\item{var_name}{the variable name or the column location.}

\item{method}{methods to calculate statistics: mean and var: BG (default), L2W; cov and cor: BG (default), BD, B, L2W.}

\item{...}{additional parameters.}

\item{var_name1}{the variable name or the column location.}

\item{var_name2}{the variable name or the column location.}
}
\value{
A numeric value: the mean, variance, covariance, or correlation.
}
\description{
Functions to compute the mean, variance, covariance, and correlation of histogram-valued data.
}
\details{
...
}
\examples{
library(HistDAWass)
}
\seealso{
int_mean int_var int_cov int_cor
}
\author{
Po-Wei Chen, Han-Ming Wu
}
