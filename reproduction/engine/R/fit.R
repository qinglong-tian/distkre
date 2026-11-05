new_dkre_fit <- function(method,
                         source_fit,
                         target_fit,
                         source_data,
                         target_data,
                         source_sizes,
                         source_weights,
                         lambda,
                         bandwidth,
                         clip,
                         metadata = list(),
                         upper_bound = NULL,
                         dimension = ncol(target_data),
                         n_target = nrow(target_data)) {
  structure(
    list(
      method = method,
      source_fit = source_fit,
      target_fit = target_fit,
      source_data = source_data,
      target_data = target_data,
      source_sizes = source_sizes,
      source_weights = source_weights,
      lambda = lambda,
      bandwidth = bandwidth,
      clip = clip,
      upper_bound = upper_bound,
      dimension = dimension,
      n_target = n_target,
      metadata = metadata
    ),
    class = "dkre_fit"
  )
}

#' Predict density ratios from a fitted model
#'
#' @param object A fitted object returned by [dkre_rff()], [fit_kulsif()],
#'   [fit_local_kulsif()], or [fit_kde_ratio()].
#' @param newdata Numeric matrix or vector of evaluation points.
#' @param ... Unused.
#' @return A numeric vector of estimated density ratios.
#' @export
predict.dkre_fit <- function(object, newdata, ...) {
  newdata <- .as_numeric_matrix(newdata, "newdata")
  dimension <- if (is.null(object$dimension)) ncol(object$target_data) else object$dimension
  if (ncol(newdata) != dimension) {
    stop("`newdata` has the wrong dimension.", call. = FALSE)
  }

  if (object$method == "kulsif") {
    alpha <- object$metadata$alpha
    k_ns <- .gaussian_kernel(newdata, object$source_data, object$bandwidth)
    k_nt <- .gaussian_kernel(newdata, object$target_data, object$bandwidth)
    fitted <- .row_sums(k_ns %*% alpha) +
      .row_sums(k_nt) / (nrow(object$target_data) * object$lambda)
  } else if (startsWith(object$method, "distributed-rff-")) {
    features <- .rff_transform(newdata, object$metadata$rff_spec)
    fitted <- as.numeric(features %*% object$metadata$coefficients)
  } else if (object$method == "kde-ratio") {
    fitted <- object$metadata$predict(newdata)
  } else {
    stop("Unsupported fitted object.", call. = FALSE)
  }

  .clip_predictions(fitted, object$clip, object$upper_bound)
}

#' Print a density-ratio fit
#'
#' @param x A fitted object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dkre_fit <- function(x, ...) {
  cat(sprintf("<dkre_fit: %s>\n", x$method))
  cat(sprintf("  source domains: %d\n", length(x$source_sizes)))
  cat(sprintf("  source observations: %d\n", sum(x$source_sizes)))
  n_target <- if (is.null(x$n_target)) nrow(x$target_data) else x$n_target
  cat(sprintf("  target observations: %d\n", n_target))
  cat(sprintf("  lambda: %.6g\n", x$lambda))
  cat(sprintf("  bandwidth: %.6g\n", x$bandwidth))
  if (!is.null(x$metadata$iterations)) {
    cat(sprintf("  iterations: %d (%s)\n", x$metadata$iterations,
                if (isTRUE(x$metadata$converged)) "converged" else
                  if (is.na(x$metadata$converged)) "fixed iterations" else "iteration limit"))
  }
  invisible(x)
}

.validate_clipping <- function(clip, upper_bound) {
  if (!is.logical(clip) || length(clip) != 1L || is.na(clip)) {
    stop("`clip` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(upper_bound)) {
    .check_positive_scalar(upper_bound, "upper_bound")
    if (!clip) stop("`upper_bound` requires `clip = TRUE`.", call. = FALSE)
  }
}

.clip_predictions <- function(values, clip, upper_bound = NULL) {
  if (clip) {
    values <- pmax(values, 0)
    if (!is.null(upper_bound)) values <- pmin(values, upper_bound)
  }
  values
}
