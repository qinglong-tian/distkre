#' Random Fourier features for a Gaussian kernel
#'
#' Generates the shared random feature map used by every client and the server.
#' For a Gaussian kernel with length scale `bandwidth`, the spectral
#' distribution is normal with variance `bandwidth^-2` in every coordinate.
#'
#' @param x A numeric matrix or vector. Rows are observations.
#' @param n_features Positive integer number of random features.
#' @param bandwidth Positive Gaussian-kernel length scale.
#' @param seed Optional integer seed. For a distributed fit, generate the
#'   specification once and reuse it at every site.
#' @param spec Optional feature specification returned in the `spec` component
#'   of an earlier call. If supplied, `x` is mapped using the same random
#'   frequencies and phases.
#' @param frequencies Optional matrix of frequency vectors, with one row per
#'   feature.
#' @param phases Optional vector of phases, with one entry per feature.
#'
#' @return If `spec` is not supplied, a list with class `dkre_rff_spec`
#'   containing `frequencies`, `phases`, `n_features`, `bandwidth`, and `seed`.
#'   If `spec` is supplied, an `n x n_features` feature matrix.
#' @export
rff_features <- function(x,
                         n_features = NULL,
                         bandwidth = 1,
                         seed = NULL,
                         spec = NULL,
                         frequencies = NULL,
                         phases = NULL) {
  if (!is.null(spec)) {
    return(.rff_transform(x, spec))
  }

  x <- .as_numeric_matrix(x, "x")
  .check_positive_scalar(bandwidth, "bandwidth")
  .check_positive_integer(n_features, "n_features")
  if (!is.null(frequencies) != !is.null(phases)) {
    stop("Supply both `frequencies` and `phases`, or neither.", call. = FALSE)
  }

  if (!is.null(frequencies)) {
    frequencies <- .as_numeric_matrix(frequencies, "frequencies")
    if (ncol(frequencies) != ncol(x)) {
      stop("`frequencies` and `x` must have the same number of columns.", call. = FALSE)
    }
    phases <- .as_numeric_vector(phases, "phases")
    if (length(phases) != nrow(frequencies)) {
      stop("`phases` must have one entry per row of `frequencies`.", call. = FALSE)
    }
  } else {
    if (!is.null(seed)) {
      with_seed(seed, {
        frequencies <- matrix(
          rnorm(n_features * ncol(x), mean = 0, sd = 1 / bandwidth),
          nrow = n_features,
          ncol = ncol(x)
        )
        phases <- runif(n_features, min = 0, max = 2 * pi)
      })
    } else {
      frequencies <- matrix(
        rnorm(n_features * ncol(x), mean = 0, sd = 1 / bandwidth),
        nrow = n_features,
        ncol = ncol(x)
      )
      phases <- runif(n_features, min = 0, max = 2 * pi)
    }
  }

  structure(
    list(
      frequencies = frequencies,
      phases = phases,
      n_features = nrow(frequencies),
      bandwidth = bandwidth,
      seed = seed
    ),
    class = "dkre_rff_spec"
  )
}

#' Return random-feature values
#'
#' Low-level helper for mapping observations with an existing specification.
#' Most users should call [rff_features()] instead.
#'
#' @param x Numeric matrix or vector.
#' @param spec An object returned by [rff_features()].
#' @return A numeric feature matrix with one row per observation.
#' @export
evaluate_rff <- function(x, spec) {
  .rff_transform(x, spec)
}

.rff_transform <- function(x, spec) {
  if (!inherits(spec, "dkre_rff_spec")) {
    stop("`spec` must be an object returned by `rff_features()`.", call. = FALSE)
  }
  x <- .as_numeric_matrix(x, "x")
  if (ncol(x) != ncol(spec$frequencies)) {
    stop("`x` has a different dimension from the feature specification.", call. = FALSE)
  }
  projection <- x %*% t(spec$frequencies) + rep(spec$phases, each = nrow(x))
  sqrt(2 / spec$n_features) * cos(projection)
}

with_seed <- function(seed, expr) {
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be a single numeric value.", call. = FALSE)
    }
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  force(expr)
}
