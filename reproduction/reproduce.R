# Single entry point. Run from the repository root.
source('reproduction/common.R')
mode <- arg('mode', 'reference'); stopifnot(mode %in% c('reference', 'smoke', 'full'))
workers <- arg('workers', '1')
out <- absolute(arg('output', paste0('reproduction/output/', mode)))
if (mode != 'reference') {
  run_command(arg('python', 'python3'), c('reproduction/fetch-data.py'))
  run_command(rscript, c('reproduction/run-simulations.R', paste0('--mode=', mode),
    paste0('--workers=', workers), paste0('--output=', file.path(out, 'simulation'))))
  run_command(rscript, c('reproduction/run-applications.R', paste0('--mode=', mode),
    paste0('--workers=', workers), paste0('--output=', file.path(out, 'applications')),
    paste0('--python=', arg('python', 'python3'))))
}
run_command(rscript, c('reproduction/render-results.R', paste0('--output=', out),
                       paste0('--input=', if (mode == 'reference') 'reference' else out)))
if (mode != 'reference') run_command(rscript, c('reproduction/verify-results.R',
  paste0('--input=', out), paste0('--scope=', mode)))
