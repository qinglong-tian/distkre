#' Median pairwise distance
#'
#' Computes the median of all pairwise Euclidean distances, evaluated on the
#' upper triangle without constructing a dense distance matrix.
#'
#' @param x Numeric matrix or vector. Rows are observations.
#' @return A non-negative scalar.
#' @export
median_pairwise_distance <- function(x) {
  x <- .as_numeric_matrix(x, "x")
  n <- nrow(x)
  if (n < 2L) {
    stop("At least two observations are required.", call. = FALSE)
  }
  count <- n * (n - 1L) / 2L
  distances <- numeric(count)
  index <- 1L
  for (i in seq_len(n - 1L)) {
    js <- seq.int(i + 1L, n)
    delta <- x[js, , drop = FALSE] - matrix(x[i, ], nrow = length(js), ncol = ncol(x), byrow = TRUE)
    distances[index:(index + length(js) - 1L)] <- sqrt(rowSums(delta^2))
    index <- index + length(js)
  }
  stats::median(distances)
}

#' Fit centralized kernel KuLSIF
#'
#' Fits the kernel unconstrained least-squares importance estimator. The
#' returned function is the unclipped density-ratio estimate unless
#' `clip = TRUE`.
#'
#' @param source_data Numeric matrix of source observations.
#' @param target_data Numeric matrix of target observations.
#' @param lambda Positive regularization parameter.
#' @param bandwidth Gaussian-kernel length scale. If `NULL`, the median
#'   pairwise target distance is used.
#' @param clip Logical; whether to truncate negative fitted values at zero.
#' @param upper_bound Optional positive B for clipping to [0,B]; requires clip.
#' @return An object of class `dkre_fit`.
#' @export
fit_kulsif <- function(source_data,
                       target_data,
                       lambda,
                       bandwidth = NULL,
                       clip = FALSE,
                       upper_bound = NULL) {
  source_data <- .as_numeric_matrix(source_data, "source_data")
  target_data <- .validate_target(target_data, ncol(source_data))
  .check_positive_scalar(lambda, "lambda")
  .validate_clipping(clip, upper_bound)
  if (is.null(bandwidth)) {
    bandwidth <- median_pairwise_distance(target_data)
  }
  .check_positive_scalar(bandwidth, "bandwidth")
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    stop("A positive finite bandwidth is required.", call. = FALSE)
  }

  k_xx <- .gaussian_kernel(source_data, source_data, bandwidth)
  k_xz <- .gaussian_kernel(source_data, target_data, bandwidth)
  k_zx <- t(k_xz)
  k_zz <- .gaussian_kernel(target_data, target_data, bandwidth)
  n <- nrow(source_data)
  m <- nrow(target_data)

  # Empirical KuLSIF objective in the paper's equation (5), expressed in its
  # equivalent RKHS representation: (Kxx / n + lambda I) alpha = -Kxz 1 / (n m lambda).
  lhs <- k_xx / n + diag(lambda, nrow = n)
  rhs <- -k_xz %*% rep(1, m) / (n * m * lambda)
  alpha <- .solve_linear(lhs, rhs)

  source_fit <- .row_sums(k_xx %*% alpha) +
    .row_sums(k_xz) / (m * lambda)
  target_fit <- .row_sums(k_zx %*% alpha) +
    .row_sums(k_zz) / (m * lambda)

  object <- new_dkre_fit(
    method = "kulsif",
    source_fit = .clip_predictions(source_fit, clip, upper_bound),
    target_fit = .clip_predictions(target_fit, clip, upper_bound),
    source_data = source_data,
    target_data = target_data,
    source_sizes = n,
    source_weights = 1,
    lambda = lambda,
    bandwidth = bandwidth,
    clip = clip,
    upper_bound = upper_bound,
    metadata = list(alpha = alpha, n = n, m = m)
  )
  object
}

#' Fit local kernel KuLSIF models
#'
#' Fits one kernel KuLSIF density-ratio estimate for each source domain against
#' the same target sample.
#'
#' @param source_data List of source-domain matrices.
#' @param target_data Numeric matrix of target observations.
#' @param lambda Positive regularization parameter, scalar or one per source.
#' @param bandwidth Optional Gaussian-kernel length scale.
#' @param clip Logical; whether to truncate negative fitted values at zero.
#' @param upper_bound Optional positive B for clipping to [0,B]; requires clip.
#' @return A list of `dkre_fit` objects.
#' @export
fit_local_kulsif <- function(source_data,
                             target_data,
                             lambda,
                             bandwidth = NULL,
                             clip = FALSE,
                             upper_bound = NULL) {
  source_data <- .validate_sources(source_data)
  target_data <- .validate_target(target_data, ncol(source_data[[1]]))
  n_sources <- length(source_data)
  if (length(lambda) == 1L) {
    lambda <- rep(lambda, n_sources)
  }
  if (length(lambda) != n_sources) {
    stop("`lambda` must be a scalar or have one entry per source.", call. = FALSE)
  }
  lapply(seq_along(source_data), function(i) {
    fit_kulsif(
      source_data[[i]],
      target_data,
      lambda = lambda[[i]],
      bandwidth = bandwidth,
      clip = clip,
      upper_bound = upper_bound
    )
  })
}

.gaussian_kernel <- function(x, y, bandwidth) {
  x_sq <- rowSums(x^2)
  y_sq <- rowSums(y^2)
  squared <- outer(x_sq, y_sq, "+") - 2 * tcrossprod(x, y)
  squared[squared < 0] <- 0
  exp(-squared / (2 * bandwidth^2))
}

.solve_linear <- function(a, b) {
  chol_a <- tryCatch(chol(a), error = function(e) NULL)
  if (is.null(chol_a)) {
    stop("Linear system is not positive definite; increase `lambda`.", call. = FALSE)
  }
  backsolve(chol_a, forwardsolve(t(chol_a), b))
}
