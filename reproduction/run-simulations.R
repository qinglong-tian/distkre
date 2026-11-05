# Rscript reproduction/run-simulations.R --mode=full --workers=1
source('reproduction/common.R')
mode <- arg('mode', 'smoke')
stopifnot(mode %in% c('smoke', 'full', 'pilot'))
workers <- as.integer(arg('workers', '1'))
out <- absolute(arg('output', paste0('reproduction/output/', mode, '/simulation')))
plan <- dget('reproduction/config/simulation.R')
config <- plan$config
cases <- plan$cases
# Include the full-kernel extension in the original deterministic computation.
# No BLAS vendor is required; the original extension used Accelerate for speed.
config$exact_kernel_cap <- if (mode == 'pilot') 2000L else 8000L
phase <- if (mode == 'pilot') 'pilot' else 'production'
if (mode == 'smoke') {
  config$features <- 12L; config$mechanism_features <- c(8L, 12L)
  config$test_points <- 200L; config$lambdas <- .01
  config$bandwidth_multipliers <- 1; config$tuning_iterations <- 10L
  config$final_iterations <- 20L; config$checkpoints <- c(0L, 1L, 2L, 10L, 20L)
  config$mechanism_lambdas <- .01; config$communication_budgets <- c(128, 640)
  cases <- cases[c(1L, 8L, 9L, 10L, 12L)]
  cases <- lapply(cases, function(x) {x$N <- 100L; x$M <- 100L; x$replications <- 1L; x})
}
selected <- arg('case', 'all')
if (selected != 'all') cases <- Filter(function(x) x$id == selected, cases)
if (!length(cases)) stop('Unknown --case.')
reps <- arg('replications', 'all')
if (reps != 'all') {
  reps <- as.integer(strsplit(reps, ',')[[1]])
  stopifnot(all(is.finite(reps)), all(reps >= 1L), !anyDuplicated(reps))
}
sources <- c('reproduction/run-simulations.R', 'reproduction/common.R',
             list.files(engine, '[.]R$', full.names = TRUE, recursive = TRUE))
record <- record_run(out, list(mode = mode, config = config, cases = cases,
                              replications = reps), sources)
setwd(engine); load_engine()
manifest <- list(config = config, cases = cases, run_id = basename(out))
atomic_rds(manifest, file.path(out, 'manifest.rds'))
jobs <- list()
for (case in cases) {
  count <- if (mode == 'pilot') config$pilot_replications else case$replications
  rr <- if (identical(reps, 'all')) seq_len(count) else reps
  if (any(rr > count)) stop('Replication exceeds the frozen count for ', case$id)
  for (i in rr) {
    path <- file.path(out, phase, case$id, sprintf('rep-%04d.rds', i))
    if (!cached_job(path)) jobs[[length(jobs) + 1L]] <- list(case = case, rep = i, path = path)
  }
}
worker <- function(job) {
  cat(job$case$id, job$rep, '\n'); flush.console()
  value <- tryCatch(sim_run_replication(job$case, job$rep, phase, config,
                                       paste(record$sources, collapse = ':')),
                    error = function(e) list(status = 'failed', error = conditionMessage(e)))
  atomic_rds(value, job$path)
  if (!identical(value$status, 'ok')) stop('Numerical job failed: ', job$path)
  TRUE
}
result <- parallel_jobs(jobs, worker, workers)
if (any(!vapply(result, isTRUE, logical(1)))) stop('A worker failed. Inspect saved jobs.')
setwd(repo)
run_command(rscript, c('reproduction/engine/simulations/summarize.R',
                      paste0('--run-directory=', out), paste0('--phase=', phase),
                      paste0('--plots=', if (mode != 'smoke' && selected == 'all' &&
                                                identical(reps, 'all')) 'true' else 'false')))
cat('Simulation outputs:', out, '\n')
