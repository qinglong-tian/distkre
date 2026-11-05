library(distkre)

set.seed(42)
sources <- list(matrix(rnorm(160), ncol = 2), matrix(rnorm(240, 0.4), ncol = 2))
target <- matrix(rnorm(160, 0.2), ncol = 2)

# Share this specification once. Each client constructs its own local state.
spec <- rff_features(target, n_features = 32, bandwidth = 1, seed = 42)
clients <- lapply(sources, dkre_client, rff_spec = spec)
sizes <- vapply(clients, function(x) x$n_source, integer(1))

# This in-process callback stands for transport to separate clients.
# A real transport must return vectors in the same order as sizes.
exchange <- function(coefficients) {
  lapply(clients, dkre_client_update, coefficients = coefficients)
}
fit <- dkre_server_fit(target, sizes, lambda = 0.1, rff_spec = spec,
                       client_update = exchange, iterations = 1000)
print(head(predict(fit, target)))
stopifnot(is.null(fit$source_data), is.null(fit$target_data))
