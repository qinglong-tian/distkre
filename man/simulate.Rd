\name{simulate_gaussian_mixture}
\alias{simulate_gaussian_mixture}
\alias{true_density_ratio}
\title{Simulation utilities for Gaussian source mixtures}
\usage{
simulate_gaussian_mixture(n_source, source_means, source_covariances,
  target_mean, target_covariance, n_target, seed = NULL)

true_density_ratio(data, source_means, source_covariances,
  source_weights, target_mean, target_covariance)
}
\arguments{
\item{n_source}{Positive source-domain sample sizes.}
\item{source_means}{Source mean vectors.}
\item{source_covariances}{Source covariance matrices.}
\item{source_weights}{Mixture weights.}
\item{target_mean}{Target mean vector.}
\item{target_covariance}{Target covariance matrix.}
\item{n_target}{Target sample size.}
\item{seed}{Optional simulation seed.}
\item{data}{Evaluation observations.}
}
\value{A simulated data list, or a vector of true density ratios.}
\description{Generates samples and analytic ratios for Gaussian source
mixtures.}
