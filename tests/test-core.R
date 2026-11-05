library(distkre)

set.seed(10)
source <- list(matrix(rnorm(40), ncol = 1), matrix(rnorm(60), ncol = 1))
target <- matrix(rnorm(50), ncol = 1)
fit <- dkre_rff(source, target, lambda = 0.2, n_features = 20,
                bandwidth = 0.8, method = "closed", seed = 7, clip = FALSE)

spec <- fit$metadata$rff_spec
z <- lapply(source, evaluate_rff, spec = spec)
zt <- evaluate_rff(target, spec)
n <- c(40, 60)
weights <- n / sum(n)
a_bar <- Reduce("+", Map(function(zi, ni, wi) wi * crossprod(zi) / ni,
                         z, n, weights))
expected <- solve(a_bar + diag(0.2, 20), colMeans(zt))
stopifnot(isTRUE(all.equal(as.numeric(fit$metadata$coefficients),
                           as.numeric(expected), tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(predict(fit, target), as.numeric(zt %*% expected),
                           tolerance = 1e-10)))

set.seed(12)
source <- list(matrix(rnorm(80, -0.5), ncol = 1),
               matrix(rnorm(70, 0.5), ncol = 1))
target <- matrix(rnorm(100), ncol = 1)
closed <- dkre_rff(source, target, lambda = 0.3, n_features = 24,
                   bandwidth = 0.7, method = "closed", seed = 3, clip = FALSE)
gd <- dkre_rff(source, target, lambda = 0.3, n_features = 24,
               bandwidth = 0.7, method = "gd", iterations = 50000,
               tolerance = 1e-8, seed = 3, clip = FALSE)
stopifnot(isTRUE(all.equal(as.numeric(gd$metadata$coefficients),
                           as.numeric(closed$metadata$coefficients),
                           tolerance = 1e-5)))

sim <- simulate_gaussian_mixture(
  n_source = c(30, 40),
  source_means = rbind(c(-0.5), c(0.5)),
  source_covariances = matrix(1),
  target_mean = 0,
  target_covariance = matrix(1),
  n_target = 50,
  seed = 1
)
kde <- fit_kde_ratio(sim$source_data, sim$target_data, method = "global")
kulsif <- fit_kulsif(do.call(rbind, sim$source_data), sim$target_data,
                    lambda = 0.1, bandwidth = 1)
stopifnot(length(predict(kde, sim$target_data)) == 50)
stopifnot(length(predict(kulsif, sim$target_data)) == 50)

cat("Core tests passed.\n")

# Algorithm 1 must match independent finite-step calculations with unequal sizes.
set.seed(42)
xs <- list(matrix(rnorm(14), 7, 2), matrix(rnorm(26), 13, 2))
zt <- matrix(rnorm(18), 9, 2)
spec <- rff_features(zt, 11, bandwidth = 0.8, seed = 13)
phi <- function(x) sqrt(2 / 11) * cos(sweep(
  x %*% t(spec$frequencies), 2, spec$phases, "+"))
F <- do.call(rbind, lapply(xs, phi))
A <- crossprod(F) / nrow(F)
a <- colMeans(phi(zt))
lam <- 0.17
clients <- lapply(xs, dkre_client, rff_spec = spec)
calls <- 0L
exchange <- function(b) {
  calls <<- calls + 1L
  lapply(clients, dkre_client_update, coefficients = b)
}
for (T in c(1L, 2L, 7L, 50L)) {
  b <- numeric(11)
  for (t in seq_len(T)) b <- as.numeric(b - (A %*% b - a + lam * b) / (2 + lam))
  calls <- 0L
  server <- dkre_server_fit(zt, c(7, 13), lam, spec, exchange,
                            iterations = T, trace = TRUE)
  stopifnot(max(abs(server$metadata$coefficients - b)) < 1e-12,
            server$metadata$iterations == T, calls == T,
            server$metadata$exchanges == T,
            server$metadata$communication_coordinates == 2 * 2 * 11 * T,
            is.na(server$metadata$converged),
            is.null(server$source_data), is.null(server$target_data),
            is.null(server$source_fit),
            server$metadata$diagnostics_iteration == T - 1L,
            length(server$metadata$history$iteration) == T)
  wrapper <- dkre_rff(xs, zt, lam, 11, method = "gd", iterations = T,
                       rff_spec = spec)
  stopifnot(max(abs(predict(wrapper, zt) - predict(server, zt))) < 1e-12,
            is.null(wrapper$source_data), is.null(wrapper$target_data),
            length(wrapper$source_fit) == 2L,
            is.null(wrapper$metadata$history))
}
# Positive tolerance explicitly enables early stopping.
calls <- 0L
stopped <- dkre_server_fit(zt, c(7, 13), lam, spec, exchange,
                           iterations = 1000, tolerance = 1e-7)
stopifnot(isTRUE(stopped$metadata$converged),
          stopped$metadata$iterations < 1000,
          stopped$metadata$exchanges == stopped$metadata$iterations + 1L,
          stopped$metadata$diagnostics_iteration == stopped$metadata$iterations)
# Saved models remain useful without observations; opt-in retention is explicit.
retained <- dkre_rff(xs, zt, lam, 11, rff_spec = spec, retain_data = TRUE)
stopifnot(identical(retained$source_data, xs), identical(retained$target_data, zt))
# Default raw output and both clipping variants must agree across prediction paths.
spec1 <- rff_features(matrix(0, 1, 1), 1, bandwidth = 1,
                       frequencies = matrix(1, 1, 1), phases = 0)
x1 <- list(matrix(c(0, pi), ncol = 1))
z1 <- matrix(0, 2, 1)
raw <- dkre_rff(x1, z1, 0.1, 1, rff_spec = spec1)
lower <- dkre_rff(x1, z1, 0.1, 1, rff_spec = spec1, clip = TRUE)
bounded <- dkre_rff(x1, z1, 0.1, 1, rff_spec = spec1,
                     clip = TRUE, upper_bound = 0.2)
stopifnot(predict(raw, matrix(pi, 1, 1)) < 0,
          predict(lower, matrix(pi, 1, 1)) == 0,
          all(bounded$target_fit <= 0.2),
          identical(bounded$source_fit[[1]], predict(bounded, x1[[1]])),
          identical(bounded$target_fit, predict(bounded, z1)))
# Reject matrix payloads instead of silently accepting full client matrices.
bad <- try(dkre_server_fit(zt, c(7, 13), lam, spec,
          function(b) lapply(clients, function(client) client$second_moment)), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
# Equation (6), independently evaluated at new points.
x <- do.call(rbind, xs)
v <- zt[1:4, , drop = FALSE]
k <- function(x, y) {
  d <- outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)
  exp(-pmax(d, 0) / (2 * 0.8^2))
}
alpha <- solve(k(x, x) / nrow(x) + diag(lam, nrow(x)),
               -rowSums(k(x, zt)) / (nrow(x) * nrow(zt) * lam))
expected <- as.numeric(k(v, x) %*% alpha + rowSums(k(v, zt)) / (nrow(zt) * lam))
exact <- fit_kulsif(x, zt, lam, bandwidth = 0.8)
stopifnot(max(abs(predict(exact, v) - expected)) < 1e-12)
local <- fit_local_kulsif(xs, zt, lam, bandwidth = 0.8,
                          clip = TRUE, upper_bound = 0.2)
stopifnot(all(vapply(local, function(f) all(predict(f, zt) <= 0.2), logical(1))))
cat("Methodology and client/server checks passed.\n")

# Even an exactly stationary iterate runs T updates when tolerance is disabled.
zero_target <- matrix(c(0, pi), ncol = 1)
zero_clients <- list(dkre_client(zero_target, spec1))
zero_exchange <- function(b) lapply(zero_clients, dkre_client_update, coefficients = b)
stationary <- dkre_server_fit(zero_target, 2, 0.1, spec1, zero_exchange,
                              iterations = 9, trace = TRUE)
stopifnot(all(stationary$metadata$coefficients == 0),
          stationary$metadata$iterations == 9,
          stationary$metadata$exchanges == 9,
          all(stationary$metadata$history$gradient_norm == 0))
# The stored diagnostic belongs to the stated iterate, not the updated one.
b_prev <- numeric(11)
for (t in seq_len(6)) b_prev <- as.numeric(
  b_prev - (A %*% b_prev - a + lam * b_prev) / (2 + lam))
diag_fit <- dkre_server_fit(zt, c(7, 13), lam, spec, exchange, iterations = 7)
J <- as.numeric(0.5 * crossprod(b_prev, A %*% b_prev) -
                crossprod(b_prev, a) + lam / 2 * sum(b_prev^2))
stopifnot(abs(diag_fit$metadata$objective - J) < 1e-12,
          diag_fit$metadata$diagnostics_iteration == 6)
round_trip <- unserialize(serialize(diag_fit, NULL))
stopifnot(identical(predict(round_trip, zt), predict(diag_fit, zt)),
          is.null(round_trip$source_data), is.null(round_trip$target_data))
cat("Fixed-point, diagnostic and saved-model checks passed.\n")
