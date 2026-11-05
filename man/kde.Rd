\name{fit_kde_ratio}
\alias{fit_kde_ratio}
\title{Fit a KDE density-ratio baseline}
\usage{
fit_kde_ratio(source_data, target_data, method = c("global", "local"),
  source_bandwidth = NULL, target_bandwidth = NULL, clip = TRUE)
}
\arguments{
\item{source_data}{A list of source-domain matrices.}
\item{target_data}{A numeric target matrix.}
\item{method}{Global pooled KDE or one denominator per source.}
\item{source_bandwidth}{Optional source bandwidths.}
\item{target_bandwidth}{Optional target bandwidth.}
\item{clip}{Whether to truncate negative predictions at zero.}
}
\value{A \code{dkre_fit} object or list of local fits.}
\description{Fits Gaussian product-kernel density-ratio baselines.}
