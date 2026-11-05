#' Distributed random-feature KuLSIF
#'
#' Local simulation wrapper for equation (13) or Algorithm 1. For independent
#' sites, use [dkre_client()], [dkre_client_update()], and [dkre_server_fit()].
#' @param source_data List of numeric source matrices.
#' @param target_data Numeric target matrix.
#' @param lambda Positive paper-scale regularization parameter.
#' @param n_features Positive random-feature count.
#' @param bandwidth Gaussian length scale; defaults to median target distance.
#' @param method `"closed"` for equation (13), `"gd"` for Algorithm 1.
#' @param iterations Number of updates, or maximum updates with early stopping.
#' @param tolerance Relative gradient tolerance; zero disables early stopping.
#' @param seed Shared random-feature seed.
#' @param rff_spec Shared specification from [rff_features()].
#' @param clip Whether to clip predictions below zero.
#' @param trace Whether to retain diagnostics evaluated before each update.
#' @param upper_bound Optional positive B for clipping to [0,B]; requires clip.
#' @param retain_data Whether to retain raw observations in the returned fit.
#' @return A `dkre_fit` object.
#' @export
dkre_rff <- function(source_data, target_data, lambda, n_features,
                     bandwidth = NULL, method = c("closed", "gd"),
                     iterations = 1000L, tolerance = 0, seed = NULL,
                     rff_spec = NULL, clip = FALSE, trace = FALSE,
                     upper_bound = NULL, retain_data = FALSE) {
  method <- match.arg(method)
  source_data <- .validate_sources(source_data)
  target_data <- .validate_target(target_data, ncol(source_data[[1]]))
  .check_positive_scalar(lambda, "lambda")
  .check_positive_integer(n_features, "n_features")
  .check_positive_integer(iterations, "iterations")
  .check_nonnegative_scalar(tolerance, "tolerance")
  .validate_clipping(clip, upper_bound)
  if (is.null(bandwidth)) {
    bandwidth <- if (!is.null(rff_spec)) rff_spec$bandwidth else
      median_pairwise_distance(target_data)
  }
  .check_positive_scalar(bandwidth, "bandwidth")
  if (is.null(rff_spec)) {
    rff_spec <- rff_features(target_data, n_features, bandwidth, seed)
  } else {
    if (!inherits(rff_spec, "dkre_rff_spec")) {
      stop("`rff_spec` must be returned by `rff_features()`.", call. = FALSE)
    }
    if (rff_spec$n_features != n_features ||
        !isTRUE(all.equal(rff_spec$bandwidth, bandwidth))) {
      stop("Feature count or bandwidth does not match `rff_spec`.", call. = FALSE)
    }
  }
  clients <- lapply(source_data, dkre_client, rff_spec = rff_spec)
  n_source <- vapply(clients, function(client) client$n_source, integer(1))
  if (method == "gd") {
    # This callback represents transport. The server receives only R-vectors.
    exchange <- function(coefficients) {
      lapply(clients, dkre_client_update, coefficients = coefficients)
    }
    fit <- dkre_server_fit(target_data, n_source, lambda, rff_spec, exchange,
                           iterations, tolerance, clip, trace, upper_bound)
  } else {
    # Equation (13) intentionally aggregates matrices, unlike Algorithm 1.
    weights <- n_source / sum(n_source)
    moment <- Reduce("+", Map(function(client, weight) {
      weight * client$second_moment
    }, clients, weights))
    z_target <- evaluate_rff(target_data, rff_spec)
    solution <- .ridge_coefficients(moment, colMeans(z_target), lambda)
    fit <- new_dkre_fit(
      method = "distributed-rff-closed", source_fit = NULL,
      target_fit = .clip_predictions(as.numeric(z_target %*% solution$coefficients),
                                     clip, upper_bound),
      source_data = NULL, target_data = NULL, source_sizes = n_source,
      source_weights = weights, lambda = lambda, bandwidth = bandwidth,
      clip = clip, upper_bound = upper_bound,
      dimension = ncol(target_data), n_target = nrow(target_data),
      metadata = c(solution, list(rff_spec = rff_spec, iterations = 0L,
        converged = TRUE, stopping_reason = "closed_form", history = NULL))
    )
  }
  # Local evaluations are a convenience of this simulation wrapper.
  fit$source_fit <- lapply(source_data, function(x) predict(fit, x))
  if (isTRUE(retain_data)) {
    fit$source_data <- source_data
    fit$target_data <- target_data
  }
  fit
}

#' Initialize a source client
#' @param source_data A numeric source matrix held at this client.
#' @param rff_spec Shared random-feature specification.
#' @return A client-local state containing its sample count and second moment.
#' @export
dkre_client <- function(source_data, rff_spec) {
  source_data <- .as_numeric_matrix(source_data, "source_data")
  .check_positive_integer(nrow(source_data), "source sample size")
  features <- evaluate_rff(source_data, rff_spec)
  structure(list(n_source = nrow(source_data),
                 second_moment = crossprod(features) / nrow(source_data)),
            class = "dkre_client")
}

#' Compute the vector a source client sends to the server
#' @param client A client-local state from [dkre_client()].
#' @param coefficients The current R-vector broadcast by the server.
#' @return The R-vector A_c b. Keep the client state at the client.
#' @export
dkre_client_update <- function(client, coefficients) {
  if (!inherits(client, "dkre_client")) {
    stop("`client` must be returned by `dkre_client()`.", call. = FALSE)
  }
  coefficients <- .as_numeric_vector(coefficients, "coefficients")
  if (length(coefficients) != ncol(client$second_moment)) {
    stop("Coefficient dimension does not match the client.", call. = FALSE)
  }
  as.numeric(client$second_moment %*% coefficients)
}

#' Run Algorithm 1 at the target server
#' @param target_data Target observations held at the server.
#' @param source_sizes Public client sample counts, in callback response order.
#' @param lambda Positive paper-scale regularization parameter.
#' @param rff_spec Shared random-feature specification.
#' @param client_update Transport callback accepting the current coefficient
#'   vector and returning a list of client R-vectors A_c b, in count order.
#' @param iterations Number of updates, or maximum updates with early stopping.
#' @param tolerance Relative gradient tolerance; zero disables early stopping.
#' @param clip Whether to clip predictions below zero.
#' @param trace Retain diagnostics evaluated before updates.
#' @param upper_bound Optional positive B; requires clip.
#' @return A prediction-ready `dkre_fit` without raw observations or client
#'   matrices. Source fitted values are NULL; evaluate those at each client.
#' @export
dkre_server_fit <- function(target_data, source_sizes, lambda, rff_spec,
                            client_update, iterations = 1000L, tolerance = 0,
                            clip = FALSE, trace = FALSE, upper_bound = NULL) {
  target_data <- .validate_target(target_data)
  .check_positive_integer(nrow(target_data), "target sample size")
  .check_positive_scalar(lambda, "lambda")
  .check_positive_integer(iterations, "iterations")
  .check_nonnegative_scalar(tolerance, "tolerance")
  .validate_clipping(clip, upper_bound)
  source_sizes <- .as_numeric_vector(source_sizes, "source_sizes")
  if (!length(source_sizes) || any(source_sizes <= 0 | source_sizes != floor(source_sizes))) {
    stop("`source_sizes` must contain positive integer counts.", call. = FALSE)
  }
  if (!is.function(client_update)) {
    stop("`client_update` must be a transport callback.", call. = FALSE)
  }
  features <- evaluate_rff(target_data, rff_spec)
  target_mean <- colMeans(features)
  weights <- source_sizes / sum(source_sizes)
  dimension <- length(target_mean)
  coefficients <- numeric(dimension)
  step_size <- 1 / (2 + lambda)
  completed <- 0L
  exchanges <- 0L
  converged <- if (tolerance > 0) FALSE else NA
  history <- if (trace) list(iteration = integer(), objective = numeric(),
                             gradient_norm = numeric(), relative_gradient = numeric()) else NULL
  for (iteration in seq_len(iterations)) {
    products <- client_update(coefficients)
    exchanges <- exchanges + 1L
    if (!is.list(products) || length(products) != length(source_sizes) ||
        !all(vapply(products, function(x) {
          is.numeric(x) && is.null(dim(x)) && length(x) == dimension && all(is.finite(x))
        }, logical(1)))) {
      stop("Callback must return one finite numeric R-vector per client, in count order.",
           call. = FALSE)
    }
    product <- Reduce("+", Map(function(x, weight) weight * x, products, weights))
    gradient <- product - target_mean + lambda * coefficients
    gradient_norm <- sqrt(sum(gradient^2))
    relative_gradient <- gradient_norm /
      max(1, sqrt(sum(coefficients^2)) + sqrt(sum(target_mean^2)))
    objective <- as.numeric(0.5 * crossprod(coefficients, product) -
      crossprod(coefficients, target_mean) + 0.5 * lambda * sum(coefficients^2))
    diagnostics_iteration <- completed
    if (trace) {
      history$iteration <- c(history$iteration, completed)
      history$objective <- c(history$objective, objective)
      history$gradient_norm <- c(history$gradient_norm, gradient_norm)
      history$relative_gradient <- c(history$relative_gradient, relative_gradient)
    }
    # Check at the evaluated iterate, without an extra diagnostic exchange.
    if (tolerance > 0 && relative_gradient <= tolerance) {
      converged <- TRUE
      break
    }
    coefficients <- coefficients - step_size * gradient
    completed <- completed + 1L
  }
  new_dkre_fit(
    method = "distributed-rff-gd", source_fit = NULL,
    target_fit = .clip_predictions(as.numeric(features %*% coefficients), clip, upper_bound),
    source_data = NULL, target_data = NULL, source_sizes = source_sizes,
    source_weights = weights, lambda = lambda, bandwidth = rff_spec$bandwidth,
    clip = clip, upper_bound = upper_bound,
    dimension = ncol(target_data), n_target = nrow(target_data),
    metadata = list(coefficients = coefficients, rff_spec = rff_spec,
      iterations = completed, converged = converged,
      stopping_reason = if (isTRUE(converged)) "gradient_tolerance" else "iteration_limit",
      objective = objective, gradient_norm = gradient_norm,
      relative_gradient = relative_gradient, diagnostics_iteration = diagnostics_iteration,
      step_size = step_size, history = history, exchanges = exchanges,
      communication_coordinates = 2 * length(source_sizes) * dimension * exchanges)
  )
}

.ridge_coefficients <- function(second_moment, target_mean, lambda) {
  chol_a <- .solve_ridge(second_moment, lambda)
  coefficients <- as.numeric(backsolve(chol_a, forwardsolve(t(chol_a), target_mean)))
  gradient <- second_moment %*% coefficients - target_mean + lambda * coefficients
  list(coefficients = coefficients,
       objective = .quadratic_objective(second_moment, coefficients, target_mean, lambda),
       gradient_norm = sqrt(sum(gradient^2)))
}

.quadratic_objective <- function(second_moment, coefficients, target_mean, lambda) {
  as.numeric(0.5 * crossprod(coefficients, second_moment %*% coefficients) -
    crossprod(coefficients, target_mean) + 0.5 * lambda * sum(coefficients^2))
}
