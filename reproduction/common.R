# Shared orchestration only; the numerical code is frozen under engine/.
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
args <- commandArgs(TRUE)
arg <- function(name, default) {
  x <- args[startsWith(args, paste0('--', name, '='))]
  if (length(x)) sub(paste0('^--', name, '='), '', tail(x, 1L)) else default
}
repo <- normalizePath('.', mustWork = TRUE)
stopifnot(file.exists('DESCRIPTION'), dir.exists('reproduction/engine'))
engine <- file.path(repo, 'reproduction', 'engine')
absolute <- function(path) {
  if (!grepl('^(/|[A-Za-z]:)', path)) path <- file.path(repo, path)
  normalizePath(path, mustWork = FALSE)
}
atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile('.saving-', tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop('Could not save ', path)
}
run_command <- function(command, arguments) {
  if (identical(command, file.path(R.home('bin'), 'Rscript')))
    arguments <- c('--vanilla', arguments)
  status <- system2(command, vapply(arguments, shQuote, character(1)))
  if (status != 0L) stop('Command failed: ', command, ' ', paste(arguments, collapse = ' '))
  invisible(status)
}
rscript <- file.path(R.home('bin'), 'Rscript')
parallel_jobs <- function(jobs, worker, workers) {
  stopifnot(length(workers) == 1L, is.finite(workers), workers >= 1L)
  if (.Platform$OS.type == 'windows' || workers == 1L) {
    if (workers > 1L) message('Using sequential execution on Windows.')
    lapply(jobs, worker)
  } else parallel::mclapply(jobs, worker, mc.cores = workers,
                           mc.set.seed = FALSE, mc.preschedule = FALSE)
}
record_run <- function(out, specification, sources) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  record <- list(specification = specification,
                 sources = tools::md5sum(sources))
  file <- file.path(out, 'configuration.rds')
  if (file.exists(file)) {
    if (!identical(readRDS(file), record))
      stop('Output contains a different configuration or source version. Use a new --output directory.')
  } else {
    atomic_rds(record, file)
    writeLines(c(paste('RNG:', paste(RNGkind(), collapse = ', ')), capture.output(sessionInfo())),
               file.path(out, 'session-info.txt'))
  }
  record
}
cached_job <- function(path) {
  if (!file.exists(path)) return(FALSE)
  value <- tryCatch(readRDS(path), error = function(e) stop('Incomplete job file: ', path))
  if (!is.null(value$error) || identical(value$status, 'failed') ||
      identical(value$status, 'method_failure'))
    stop('A failed job is saved at ', path, '. Inspect it before retrying in a new output directory.')
  TRUE
}
load_engine <- function() {
  for (f in list.files('R', '[.]R$', full.names = TRUE)) source(f, local = .GlobalEnv)
  source('simulations/configs/study.R', local = .GlobalEnv)
  source('simulations/R/core.R', local = .GlobalEnv)
  source('simulations/R/studies.R', local = .GlobalEnv)
}
