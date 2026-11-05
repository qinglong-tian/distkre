# Frozen paper settings: 80 hospital/station fixtures + 20 spam fixtures.
source('reproduction/common.R')
mode <- arg('mode', 'smoke'); stopifnot(mode %in% c('smoke', 'full'))
workers <- as.integer(arg('workers', '1'))
output <- absolute(arg('output', paste0('reproduction/output/', mode, '/applications')))
datasets <- strsplit(arg('datasets', 'heart,beijing,spambase'), ',')[[1]]
stopifnot(length(datasets) > 0L, !anyDuplicated(datasets),
          all(datasets %in% c('heart', 'beijing', 'spambase')))
run_command(arg('python', 'python3'), c('reproduction/fetch-data.py', '--check', datasets))
sources <- c('reproduction/run-applications.R', 'reproduction/common.R',
             'reproduction/collect-applications.R', 'reproduction/fetch-data.py', 'reproduction/datasets.json',
             list.files(engine, '[.]R$', full.names = TRUE, recursive = TRUE))
config <- list(seeds = 1:5, features = 128L, iterations = 3000L,
               lambdas = c(.01, .1), bandwidths = c(.5, 1, 2), cap = 80L,
               predictor_ridge = .01)
spam_config <- list(phase = 'confirm', seeds = 201:220, features = 512L,
                    iterations = 5000L, lambdas = c(.001, .01, .1),
                    bandwidths = c(.5, 1, 2), practical_AUPRC_margin = .025,
                    version = 'pu-local-v1')
inputs <- list.files(file.path(engine, 'applications/data/raw'), full.names = TRUE, recursive = TRUE)
inputs <- inputs[!grepl('[.]zip$|download-provenance[.]json$', inputs)]
run_record <- record_run(output, list(mode = mode, datasets = datasets, config = config,
                       spam_config = spam_config, inputs = tools::md5sum(inputs)), sources)
setwd(engine)
for (f in list.files('R', '[.]R$', full.names = TRUE)) source(f)
for (f in c('applications/R/data.R', 'applications/R/protocol.R', 'applications/R/expansion.R',
            'applications/local-comparison/fixtures.R', 'applications/local-comparison/methods.R',
            'applications/local-comparison/pu.R', 'applications/natural-domains/data.R')) source(f)
natural <- intersect(datasets, c('heart', 'beijing'))
if (length(natural)) {
  data <- setNames(lapply(natural, app_dataset), natural)
  out <- file.path(output, 'natural'); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  source('natural-job.R')
  jobs <- list()
  for (name in natural) {
    targets <- sort(unique(data[[name]]$site))
    if (mode == 'smoke') targets <- head(targets, 1L)
    for (target in targets) for (seed in if (mode == 'smoke') 1L else 1:5) {
      job <- list(name = name, target = target, seed = seed)
      path <- file.path(out, paste0(paste(name, target, seed, sep = '-'), '.rds'))
      if (!cached_job(path)) jobs[[length(jobs) + 1L]] <- job
    }
  }
  result <- parallel_jobs(jobs, function(job) {
    run(job)
    path <- file.path(out, paste0(paste(job$name, job$target, job$seed, sep = '-'), '.rds'))
    cached_job(path)
  }, workers)
  if (any(!vapply(result, isTRUE, logical(1)))) stop('A natural-domain worker failed.')
}
if ('spambase' %in% datasets) {
  d <- app_dataset('spambase'); out <- file.path(output, 'spam')
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  seeds <- if (mode == 'smoke') 201L else 201:220
  pending <- Filter(function(seed) !cached_job(file.path(out, paste0('spam-', seed, '.rds'))), as.list(seeds))
  result <- parallel_jobs(pending, function(seed) {
    cat('spam', seed, '\n'); flush.console()
    value <- tryCatch(lc_compare_pu(lc_pu_fixture(d, seed), spam_config),
                      error = function(e) list(error = conditionMessage(e)))
    atomic_rds(value, file.path(out, paste0('spam-', seed, '.rds')))
    if (!is.null(value$error)) stop(value$error)
    TRUE
  }, workers)
  if (any(!vapply(result, isTRUE, logical(1)))) stop('A spam worker failed.')
}
setwd(repo)
run_command(rscript, c('reproduction/collect-applications.R', paste0('--input=', output)))
cat('Application outputs:', output, '\n')
