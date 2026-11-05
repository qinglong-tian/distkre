\name{rff_features}
\alias{rff_features}
\alias{evaluate_rff}
\title{Shared random Fourier features}
\usage{
rff_features(x, n_features = NULL, bandwidth = 1, seed = NULL,
  spec = NULL, frequencies = NULL, phases = NULL)

evaluate_rff(x, spec)
}
\arguments{
\item{x}{A numeric matrix or vector of observations.}
\item{n_features}{Number of random features.}
\item{bandwidth}{Positive Gaussian-kernel length scale.}
\item{seed}{Optional seed.}
\item{spec}{An existing feature specification.}
\item{frequencies}{Optional frequency matrix.}
\item{phases}{Optional phase vector.}
}
\value{A feature specification, or a feature matrix when \code{spec} is
supplied or \code{evaluate_rff} is called.}
\description{Creates and applies a reproducible shared feature map.}
