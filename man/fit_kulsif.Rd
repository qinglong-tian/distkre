\name{fit_kulsif}
\alias{fit_kulsif}
\alias{fit_local_kulsif}
\title{Fit kernel KuLSIF baselines}
\usage{
fit_kulsif(source_data, target_data, lambda, bandwidth = NULL,
  clip = FALSE, upper_bound = NULL)

fit_local_kulsif(source_data, target_data, lambda, bandwidth = NULL,
  clip = FALSE, upper_bound = NULL)
}
\arguments{
\item{source_data}{A source matrix, or a list of source-domain matrices
for \code{fit_local_kulsif}.}
\item{target_data}{A numeric target matrix.}
\item{lambda}{Positive regularization parameter or one per local fit.}
\item{bandwidth}{Optional Gaussian-kernel length scale.}
\item{clip}{Whether to truncate negative predictions at zero.}
\item{upper_bound}{Optional positive B for clipping to [0,B]; requires clip.}
}
\value{A \code{dkre_fit} object, or a list of such objects for local fits.}
\description{Implements the centralized kernel KuLSIF estimator.}
