.as_numeric_matrix <- function(x, name) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  } else if (length(dim(x)) != 2L) {
    stop(sprintf("`%s` must be a vector, matrix, or data frame.", name), call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` contains missing or non-finite values.", name), call. = FALSE)
  }
  matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
}

.as_numeric_vector <- function(x, name) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  }
  d <- dim(x)
  if (!is.null(d) && sum(d > 1L) != 1L) {
    stop(sprintf("`%s` must be a numeric vector.", name), call. = FALSE)
  }
  x <- as.numeric(x)
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` contains missing or non-finite values.", name), call. = FALSE)
  }
  x
}

.check_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a finite positive number.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.check_nonnegative_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop(sprintf("`%s` must be a finite non-negative number.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.check_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != floor(x)) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.validate_sources <- function(source_list) {
  if (!is.list(source_list) || length(source_list) == 0L) {
    stop("`source_data` must be a non-empty list of numeric matrices.", call. = FALSE)
  }
  source_list <- lapply(seq_along(source_list), function(i) {
    .as_numeric_matrix(source_list[[i]], sprintf("source_data[[%d]]", i))
  })
  if (any(vapply(source_list, nrow, integer(1)) == 0L)) {
    stop("Source samples must be non-empty.", call. = FALSE)
  }
  dimensions <- vapply(source_list, ncol, integer(1))
  if (length(unique(dimensions)) != 1L) {
    stop("All source data sets must have the same number of columns.", call. = FALSE)
  }
  source_list
}

.validate_target <- function(target_data, dimension = NULL) {
  target_data <- .as_numeric_matrix(target_data, "target_data")
  .check_positive_integer(nrow(target_data), "target sample size")
  if (!is.null(dimension) && ncol(target_data) != dimension) {
    stop("Source and target data must have the same number of columns.", call. = FALSE)
  }
  target_data
}

.source_sizes <- function(source_list) {
  vapply(source_list, nrow, integer(1))
}

.source_mixture_weights <- function(source_list) {
  source_list <- .validate_sources(source_list)
  .source_sizes(source_list) / sum(.source_sizes(source_list))
}

.row_sums <- function(x) {
  if (ncol(x) == 1L) {
    as.numeric(x)
  } else {
    rowSums(x)
  }
}

.column_means <- function(x) {
  if (nrow(x) == 1L) {
    as.numeric(x)
  } else {
    colMeans(x)
  }
}

.solve_ridge <- function(a, lambda) {
  d <- ncol(a)
  a <- a + diag(lambda, nrow = d)
  chol_a <- tryCatch(chol(a), error = function(e) NULL)
  if (is.null(chol_a)) {
    stop("The ridge system is not positive definite. Increase `lambda`.", call. = FALSE)
  }
  chol_a
}
