\name{dkre_server_fit}
\alias{dkre_server_fit}
\alias{dkre_client}
\alias{dkre_client_update}
\title{Client and server interfaces for Algorithm 1}
\usage{
dkre_client(source_data, rff_spec)

dkre_client_update(client, coefficients)

dkre_server_fit(target_data, source_sizes, lambda, rff_spec,
  client_update, iterations = 1000L, tolerance = 0,
  clip = FALSE, trace = FALSE, upper_bound = NULL)
}
\arguments{
\item{source_data}{Non-empty numeric source matrix held at one client.}
\item{rff_spec}{Identical Gaussian random-feature specification at all sites.}
\item{client}{Client-local state returned by \code{dkre_client}.}
\item{coefficients}{Current coefficient R-vector broadcast by the server.}
\item{target_data}{Non-empty numeric target matrix held at the server.}
\item{source_sizes}{Positive integer client counts, in callback response order.}
\item{lambda}{Positive paper-scale regularization parameter.}
\item{client_update}{Transport callback accepting the coefficient vector and
returning a list of finite numeric R-vectors A_c b, in source-count order.}
\item{iterations}{Exact update count, or maximum with positive tolerance.}
\item{tolerance}{Relative gradient tolerance; zero disables early stopping.}
\item{clip}{Whether to clip negative predictions to zero.}
\item{trace}{Whether to retain diagnostics evaluated before updates.}
\item{upper_bound}{Optional positive B; requires clip.}
}
\value{\code{dkre_client} returns a client-local state containing its sample
count and second moment. \code{dkre_client_update} returns only an R-vector.
\code{dkre_server_fit} returns a prediction-ready \code{dkre_fit} with no raw
observations or client matrices. Source fitted values are NULL; evaluate them
at each client using \code{predict}.}
\description{Implements Algorithm 1 using vector exchanges. Keep client states
at the source sites. The server receives public counts, the shared feature
specification, and vectors through a user-supplied transport callback.}
\details{Transport is supplied by the caller; the package does not launch remote
processes or provide network transport. The local wrapper \code{dkre_rff}
simulates the sites through this interface. See its documentation for stopping,
diagnostics and communication accounting. All sites must use the identical
feature specification and stable client ordering throughout the run.}
\examples{
sources <- list(matrix(c(-1, 0, 1), ncol = 1),
                matrix(c(0, 1, 2, 3), ncol = 1))
target <- matrix(c(0, 0.5, 1), ncol = 1)
spec <- rff_features(target, n_features = 8, bandwidth = 1, seed = 42)
# Local simulation: each state belongs to its respective client.
clients <- lapply(sources, dkre_client, rff_spec = spec)
exchange <- function(b) lapply(clients, dkre_client_update, coefficients = b)
fit <- dkre_server_fit(target, c(3, 4), lambda = 0.1, rff_spec = spec,
                       client_update = exchange, iterations = 20)
predict(fit, target)
}
