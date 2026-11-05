\name{dkre_rff}
\alias{dkre_rff}
\title{Distributed random-feature KuLSIF}
\usage{
dkre_rff(source_data, target_data, lambda, n_features,
  bandwidth = NULL, method = c("closed", "gd"), iterations = 1000L,
  tolerance = 0, seed = NULL, rff_spec = NULL, clip = FALSE,
  trace = FALSE, upper_bound = NULL, retain_data = FALSE)
}
\arguments{
\item{source_data}{A non-empty list of non-empty numeric source matrices.}
\item{target_data}{A non-empty numeric target matrix.}
\item{lambda}{Positive paper-scale regularization parameter.}
\item{n_features}{Positive number of random Fourier features.}
\item{bandwidth}{Gaussian length scale. Defaults to the supplied feature
specification's bandwidth, or the median target pairwise distance.}
\item{method}{\code{"closed"} for equation (13), or \code{"gd"} for Algorithm 1.}
\item{iterations}{Exact update count when tolerance is zero; otherwise a maximum.}
\item{tolerance}{Relative gradient tolerance. Zero disables early stopping.}
\item{seed}{Optional seed for shared random features.}
\item{rff_spec}{Optional specification from \code{rff_features}.}
\item{clip}{Whether to clip predictions below zero. Defaults to FALSE.}
\item{trace}{Whether to retain diagnostics evaluated before updates.}
\item{upper_bound}{Optional positive B for clipping to [0,B]; requires clip.}
\item{retain_data}{Whether to retain raw source and target observations.}
}
\value{A \code{dkre_fit} object with source and target predictions. Raw
observations are omitted by default. Predictions remain available at new points.}
\description{Local simulation wrapper for normalized random-feature KuLSIF.
The gradient route uses the same vector-only interface as \code{dkre_server_fit}.
The closed-form route intentionally aggregates second-moment matrices.}
\details{The step size is 1/(2+lambda). With zero tolerance, exactly T updates
and T client exchanges occur, even if the gradient is zero. Diagnostics describe
the last evaluated iterate, before the last update; its index is recorded in
\code{metadata$diagnostics_iteration}. No extra exchange is used for diagnostics.
With positive tolerance, the server checks the gradient at the current iterate
before updating. \code{metadata$exchanges} includes a final stopping check if any.
\code{metadata$communication_coordinates} counts broadcast and response vector
coordinates, excluding initialization, transport overhead and prediction.}
\seealso{\code{\link{dkre_server_fit}}, \code{\link{dkre_client}}}
