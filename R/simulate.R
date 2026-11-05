#' Simulate a Gaussian source mixture and target sample
#'
#' @param n_source Integer vector of source-domain sample sizes.
#' @param source_means Numeric matrix or list of source mean vectors.
#' @param source_covariances Covariance matrix shared by all sources, or a list
#'   with one covariance per source.
#' @param target_mean Target mean vector.
#' @param target_covariance Target covariance matrix.
#' @param n_target Target sample size.
#' @param seed Optional seed.
#' @return A list with `source_data`, `target_data`, `source_means`,
#'   `source_covariances`, `target_mean`, `target_covariance`, and
#'   `source_weights`.
#' @export
simulate_gaussian_mixture <- function(n_source,
                                      source_means,
                                      source_covariances,
                                      target_mean,
                                      target_covariance,
                                      n_target,
                                      seed = NULL) {
  n_source <- as.integer(n_source)
  if (length(n_source) < 1L || anyNA(n_source) || any(n_source < 1L)) {
    stop("`n_source` must contain positive integers.", call. = FALSE)
  }
  if (is.matrix(source_means)) {
    source_means <- lapply(seq_len(nrow(source_means)), function(i) source_means[i, ])
  }
  if (!is.list(source_means) || length(source_means) != length(n_source)) {
    stop("`source_means` must match the number of sources.", call. = FALSE)
  }
  dimension <- length(source_means[[1]])
  source_means <- lapply(source_means, function(mu) {
    mu <- .as_numeric_vector(mu, "source_means element")
    if (length(mu) != dimension) stop("Source means have inconsistent dimensions.", call. = FALSE)
    mu
  })
  target_mean <- .as_numeric_vector(target_mean, "target_mean")
  if (length(target_mean) != dimension) stop("Target mean has the wrong dimension.", call. = FALSE)
  if (is.matrix(source_covariances)) {
    source_covariances <- rep(list(source_covariances), length(n_source))
  }
  if (!is.list(source_covariances) || length(source_covariances) != length(n_source)) {
    stop("`source_covariances` must match the number of sources.", call. = FALSE)
  }
  source_covariances <- lapply(source_covariances, .as_covariance)
  target_covariance <- .as_covariance(target_covariance)

  generate <- function() {
    source_data <- lapply(seq_along(n_source), function(i) {
      .rmvnorm(n_source[[i]], source_means[[i]], source_covariances[[i]])
    })
    target_data <- .rmvnorm(as.integer(n_target), target_mean, target_covariance)
    list(
      source_data = source_data,
      target_data = target_data,
      source_means = source_means,
      source_covariances = source_covariances,
      target_mean = target_mean,
      target_covariance = target_covariance,
      source_weights = n_source / sum(n_source)
    )
  }
  if (is.null(seed)) generate() else with_seed(seed, generate())
}

#' True Gaussian-mixture density ratio
#'
#' @param data Evaluation observations.
#' @param source_means List of source means.
#' @param source_covariances List of source covariance matrices.
#' @param source_weights Mixture weights.
#' @param target_mean Target mean.
#' @param target_covariance Target covariance.
#' @return A numeric vector of target-to-pooled-source density ratios.
#' @export
true_density_ratio <- function(data,
                               source_means,
                               source_covariances,
                               source_weights,
                               target_mean,
                               target_covariance) {
  data <- .as_numeric_matrix(data, "data")
  source_means <- lapply(source_means, .as_numeric_vector, name = "source mean")
  source_covariances <- lapply(source_covariances, .as_covariance)
  target_mean <- .as_numeric_vector(target_mean, "target_mean")
  target_covariance <- .as_covariance(target_covariance)
  source_weights <- .as_numeric_vector(source_weights, "source_weights")
  if (length(source_means) != length(source_weights) ||
      length(source_covariances) != length(source_weights)) {
    stop("Source means, covariances, and weights must have the same length.", call. = FALSE)
  }
  if (!isTRUE(all.equal(sum(source_weights), 1))) {
    stop("`source_weights` must sum to one.", call. = FALSE)
  }

  source_density <- vapply(seq_along(source_means), function(i) {
    .dmvnorm(data, source_means[[i]], source_covariances[[i]])
  }, numeric(nrow(data)))
  target_density <- .dmvnorm(data, target_mean, target_covariance)
  as.numeric(target_density / (source_density %*% source_weights))
}

.as_covariance <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) != ncol(x) || any(!is.finite(x))) {
    stop("A covariance must be a finite square numeric matrix.", call. = FALSE)
  }
  x
}

.rmvnorm <- function(n, mean, covariance) {
  dimension <- length(mean)
  if (dimension == 1L) {
    return(matrix(stats::rnorm(n, mean = mean, sd = sqrt(covariance[1, 1])), ncol = 1L))
  }
  factor <- tryCatch(chol(covariance), error = function(e) {
    stop("Covariance must be positive definite.", call. = FALSE)
  })
  z <- matrix(stats::rnorm(n * dimension), nrow = n, ncol = dimension)
  sweep(z %*% factor, 2L, mean, "+")
}

.dmvnorm <- function(x, mean, covariance) {
  dimension <- length(mean)
  if (dimension == 1L) {
    return(stats::dnorm(x[, 1], mean = mean, sd = sqrt(covariance[1, 1])))
  }
  inverse <- tryCatch(solve(covariance), error = function(e) {
    stop("Covariance must be invertible.", call. = FALSE)
  })
  determinant <- determinant(covariance, logarithm = TRUE)$modulus[[1]]
  centered <- sweep(x, 2L, mean, "-")
  exponent <- rowSums((centered %*% inverse) * centered)
  exp(-0.5 * (dimension * log(2 * pi) + determinant + exponent))
}
