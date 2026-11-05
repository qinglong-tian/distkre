#' Fit a KDE density-ratio baseline
#'
#' Uses Gaussian product kernels. For `method = "global"`, source kernels are
#' pooled with mixture weights `n_c / n`. For `method = "local"`, each source
#' domain receives its own denominator kernel.
#'
#' @param source_data List of source-domain matrices.
#' @param target_data Numeric target matrix.
#' @param method Either `"global"` or `"local"`.
#' @param source_bandwidth Optional source bandwidth, scalar or one per source.
#' @param target_bandwidth Optional target bandwidth.
#' @param clip Logical; whether to truncate negative fitted values at zero.
#' @return An object of class `dkre_fit`.
#' @export
fit_kde_ratio <- function(source_data,
                          target_data,
                          method = c("global", "local"),
                          source_bandwidth = NULL,
                          target_bandwidth = NULL,
                          clip = TRUE) {
  method <- match.arg(method)
  source_data <- .validate_sources(source_data)
  target_data <- .validate_target(target_data, ncol(source_data[[1]]))
  n_source <- .source_sizes(source_data)
  source_weights <- n_source / sum(n_source)

  if (is.null(target_bandwidth)) {
    target_bandwidth <- .kde_bandwidth(target_data)
  }
  .check_positive_scalar(target_bandwidth, "target_bandwidth")
  if (is.null(source_bandwidth)) {
    source_bandwidth <- if (method == "global") {
      .kde_bandwidth(do.call(rbind, source_data))
    } else {
      vapply(source_data, .kde_bandwidth, numeric(1))
    }
  }
  if (length(source_bandwidth) == 1L) {
    source_bandwidth <- rep(source_bandwidth, length(source_data))
  }
  if (length(source_bandwidth) != length(source_data) ||
      any(!is.finite(source_bandwidth)) || any(source_bandwidth <= 0)) {
    stop("`source_bandwidth` must be positive and match the number of sources.", call. = FALSE)
  }

  if (method == "local") {
    denominator <- lapply(seq_along(source_data), function(i) {
      .density_provider(source_data[[i]], source_bandwidth[[i]])
    })
    local_fits <- lapply(seq_along(source_data), function(i) {
      .kde_ratio_object(
        numerator = .density_provider(target_data, target_bandwidth),
        denominator = denominator[[i]],
        source_data = source_data,
        target_data = target_data,
        source_sizes = n_source,
        source_weights = source_weights,
        bandwidth = source_bandwidth[[i]],
        clip = clip
      )
    })
    return(local_fits)
  }

  numerator <- .density_provider(target_data, target_bandwidth)
  denominator <- function(points) {
    if (is.null(dim(points))) {
      points <- matrix(points, ncol = ncol(target_data))
    }
    values <- vapply(seq_along(source_data), function(i) {
      .density_provider(source_data[[i]], source_bandwidth[[i]])(points)
    }, numeric(nrow(points)))
    as.numeric(values %*% source_weights)
  }

  .kde_ratio_object(
    numerator = numerator,
    denominator = denominator,
    source_data = source_data,
    target_data = target_data,
    source_sizes = n_source,
    source_weights = source_weights,
    bandwidth = median(source_bandwidth),
    clip = clip
  )
}

.kde_ratio_object <- function(numerator,
                              denominator,
                              source_data,
                              target_data,
                              source_sizes,
                              source_weights,
                              bandwidth,
                              clip) {
  predict_fn <- function(points) {
    numerator(points) / denominator(points)
  }
  source_fit <- lapply(source_data, predict_fn)
  new_dkre_fit(
    method = "kde-ratio",
    source_fit = source_fit,
    target_fit = predict_fn(target_data),
    source_data = source_data,
    target_data = target_data,
    source_sizes = source_sizes,
    source_weights = source_weights,
    lambda = NA_real_,
    bandwidth = bandwidth,
    clip = clip,
    metadata = list(predict = predict_fn)
  )
}

.density_provider <- function(data, bandwidth) {
  data <- .as_numeric_matrix(data, "data")
  n <- nrow(data)
  dimension <- ncol(data)
  normalization <- (2 * pi * bandwidth^2)^(-dimension / 2)
  function(points) {
    points <- .as_numeric_matrix(points, "points")
    if (ncol(points) != dimension) {
      stop("Evaluation points have the wrong dimension.", call. = FALSE)
    }
    kernel_matrix <- .gaussian_kernel(points, data, bandwidth)
    normalization * .row_sums(kernel_matrix) / n
  }
}

.kde_bandwidth <- function(x) {
  x <- .as_numeric_matrix(x, "x")
  n <- nrow(x)
  dimension <- ncol(x)
  if (n < 2L) {
    stop("At least two observations are required for a KDE bandwidth.", call. = FALSE)
  }
  # A stable multivariate normal-reference rule. It reduces to the usual
  # n^(-1/(d+4)) scaling and avoids the inflated factors in the prototype.
  scale <- vapply(seq_len(dimension), function(j) stats::sd(x[, j]), numeric(1))
  scale[!is.finite(scale) | scale <= 0] <- 1
  mean(scale) * n^(-1 / (dimension + 4))
}
