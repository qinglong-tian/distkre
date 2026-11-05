# Compare fresh computations with archived paper metrics, ignoring timings.
source('reproduction/common.R')
input <- absolute(arg('input', 'reproduction/output/full'))
scope <- arg('scope', 'full')
stopifnot(scope %in% c('full', 'smoke', 'applications'))
read_ref <- function(path) {
  con <- gzfile(path, 'rt'); on.exit(close(con)); read.csv(con)
}
compare <- function(actual, reference, keys, fields, complete) {
  key <- function(z) apply(z[, keys, drop = FALSE], 1, paste, collapse = '|')
  ka <- key(actual); kr <- key(reference)
  stopifnot(!anyDuplicated(ka), !anyDuplicated(kr))
  ii <- match(ka, kr)
  if (anyNA(ii) || (complete && !setequal(ka, kr))) stop('Job/metric collection differs from the paper.')
  for (field in fields) {
    a <- actual[[field]]; b <- reference[[field]][ii]
    stopifnot(identical(is.na(a), is.na(b)))
    ok <- !is.na(a)
    if (any(abs(a[ok] - b[ok]) > 1e-8 * pmax(1, abs(b[ok]))))
      stop('Numerical mismatch in ', field)
  }
  cat('Matched', nrow(actual), 'rows against archived metrics.\n')
}
if (scope == 'full') {
  actual <- read.csv(file.path(input, 'simulation/summary/production/replicate-metrics.csv'))
  reference <- read_ref('reproduction/reference/simulation/replicate-metrics.csv.gz')
  compare(actual, reference,
    c('case_id', 'replication', 'method', 'row_type', 'R', 'lambda', 'iterations', 'budget'),
    c('rmse', 'mse', 'validation', 'fit_coordinates', 'tuning_coordinates'), TRUE)
}
actual <- read.csv(file.path(input, 'applications/natural/metrics.csv'))
reference <- read_ref('reproduction/reference/natural/metrics.csv.gz')
compare(actual, reference, c('dataset', 'target', 'seed', 'method', 'task'), 'value', scope != 'smoke')
actual <- read.csv(file.path(input, 'applications/spam/metrics.csv'))
reference <- read_ref('reproduction/reference/spam/metrics.csv.gz')
compare(actual, reference, c('replication', 'allocation', 'method', 'metric'), 'value', scope != 'smoke')
cat('Verification passed; relative/absolute tolerance 1e-8. Timings are excluded.\n')
