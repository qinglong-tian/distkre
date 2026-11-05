# Integration checks against archived numerical results, not a second copy
# of the fitting implementation. Run from the repository root.
source('reproduction/common.R')
run_command(rscript, c('reproduction/run-simulations.R', '--mode=full',
  '--case=accuracy_N1000_M1000', '--replications=1',
  '--output=reproduction/output/replay/simulation'))
run_command(rscript, c('reproduction/run-applications.R', '--mode=smoke',
  '--output=reproduction/output/replay/applications'))
read_reference <- function(path) {
  con <- gzfile(path, 'rt'); on.exit(close(con)); read.csv(con)
}
expected <- read_reference('reproduction/reference/simulation/replicate-metrics.csv.gz')
expected <- subset(expected, case_id == 'accuracy_N1000_M1000' & replication == 1)
actual <- readRDS('reproduction/output/replay/simulation/production/accuracy_N1000_M1000/rep-0001.rds')$rows
ii <- match(actual$method, expected$method)
stopifnot(!anyNA(ii), nrow(actual) == nrow(expected))
for (field in c('rmse', 'mse', 'lambda', 'sigma', 'validation', 'iterations'))
  stopifnot(isTRUE(all.equal(actual[[field]], expected[[field]][ii], tolerance = 1e-8)))
expected <- read_reference('reproduction/reference/natural/metrics.csv.gz')
actual <- read.csv('reproduction/output/replay/applications/natural/metrics.csv')
key <- function(z) paste(z$dataset, z$target, z$seed, z$method, z$task)
ii <- match(key(actual), key(expected))
stopifnot(!anyNA(ii), nrow(actual) == 30L,
          isTRUE(all.equal(actual$value, expected$value[ii], tolerance = 1e-8)))
expected <- read_reference('reproduction/reference/spam/metrics.csv.gz')
actual <- read.csv('reproduction/output/replay/applications/spam/metrics.csv')
key <- function(z) paste(z$replication, z$allocation, z$method, z$metric)
ii <- match(key(actual), key(expected))
stopifnot(!anyNA(ii), nrow(actual) == 48L,
          isTRUE(all.equal(actual$value, expected$value[ii], tolerance = 1e-8)))
cat('Published simulation, hospital, station and spam fixtures replayed successfully.\n')
