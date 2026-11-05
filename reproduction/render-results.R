# Rebuild tables and figures from saved job metrics; no statistical fitting here.
source('reproduction/common.R')
input <- arg('input', 'reference')
out <- absolute(arg('output', 'reproduction/output/reference'))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_table <- function(path) {
  con <- if (grepl('[.]gz$', path)) gzfile(path, 'rt') else file(path, 'rt')
  on.exit(close(con)); read.csv(con, stringsAsFactors = FALSE)
}
plan <- dget('reproduction/config/simulation.R')
if (input == 'reference') {
  sim <- read_table('reproduction/reference/simulation/replicate-metrics.csv.gz')
  natural <- read_table('reproduction/reference/natural/metrics.csv.gz')
  spam <- read_table('reproduction/reference/spam/metrics.csv.gz')
  simdir <- file.path(out, 'simulation')
  config <- plan$config; config$exact_kernel_cap <- 8000L
  atomic_rds(list(config = config, cases = plan$cases, run_id = 'paper-reference'),
             file.path(simdir, 'manifest.rds'))
  # The historical analysis reads per-job files. Reconstruct only its metric
  # inputs, without publishing raw subjects, fitted models, or machine paths.
  groups <- split(seq_len(nrow(sim)), paste(sim$case_id, sim$replication, sep = '|'))
  for (ii in groups) {
    z <- sim[ii, , drop = FALSE]
    case <- Filter(function(c) c$id == z$case_id[1], plan$cases)[[1L]]
    atomic_rds(list(status = 'ok', rows = z, metadata = list(case = case, elapsed = 0)),
               file.path(simdir, 'production', case$id, sprintf('rep-%04d.rds', z$replication[1])))
  }
  run_command(rscript, c('reproduction/engine/simulations/summarize.R',
                         paste0('--run-directory=', simdir), '--phase=production'))
} else {
  input <- absolute(input); simdir <- file.path(input, 'simulation')
  natural_path <- file.path(input, 'applications/natural/metrics.csv')
  spam_path <- file.path(input, 'applications/spam/metrics.csv')
  natural <- if (file.exists(natural_path)) read_table(natural_path) else NULL
  spam <- if (file.exists(spam_path)) read_table(spam_path) else NULL
}
table_dir <- file.path(out, 'tables'); dir.create(table_dir, showWarnings = FALSE)
tables <- list()
if (file.exists(file.path(simdir, 'summary/production/method-summary.csv'))) {
  s <- read_table(file.path(simdir, 'summary/production/method-summary.csv'))
  cell <- function(z, method) {
    a <- z[z$method == method, , drop = FALSE]
    if (!nrow(a)) return('--')
    sprintf('%.4f (%.4f)', a$rmse_mean[1], a$rmse_mcse[1])
  }
  accuracy <- list()
  for (case in Filter(function(x) x$study == 'accuracy', plan$cases)) {
    z <- subset(s, case_id == case$id & row_type == 'estimate')
    if (!nrow(z)) next
    accuracy[[length(accuracy) + 1L]] <- data.frame(N = case$N, M = case$M,
      L = z$successful[z$method == 'distributed_gd'][1], GD = cell(z, 'distributed_gd'),
      RF = cell(z, 'rff_closed'), Kernel = cell(z, 'kernel_pooled'))
  }
  if (length(accuracy)) tables$simulation_accuracy <- do.call(rbind, accuracy)
  for (study in c('heterogeneity', 'imbalance')) {
    values <- list()
    for (id in vapply(Filter(function(c) c$study == study, plan$cases), `[[`, character(1), 'id')) {
      z <- subset(s, case_id == id & row_type == 'estimate')
      values[[length(values) + 1L]] <- data.frame(Setting = id,
        GD = cell(z, 'distributed_gd'), RF = cell(z, 'rff_closed'),
        Kernel = cell(z, 'kernel_pooled'), Local = cell(z, 'kernel_local_naive'),
        Uniform = cell(z, 'gd_uniform_weight_ablation'))
    }
    if (length(values)) {
      tab <- do.call(rbind, values)
      if (study == 'heterogeneity') {
        tab$Setting <- sub('heterogeneity_', '', tab$Setting)
        names(tab)[1] <- 'tau'
        tab$Uniform <- NULL
      }
      tables[[paste0('simulation_', study)]] <- tab
    }
  }
  values <- list()
  z <- subset(s, study == 'communication' & group_lambda == .01)
  for (features in sort(unique(z$R))) {
    ref <- subset(z, R == features & row_type == 'reference')
    path <- subset(z, R == features & row_type == 'checkpoint' & group_iterations > 0)
    path <- path[order(path$group_iterations), , drop = FALSE]
    good <- path[path$rmse_mean <= 1.02 * ref$rmse_mean[1], , drop = FALSE]
    if (!nrow(good)) next
    values[[length(values) + 1L]] <- data.frame(R = features,
      RF_RMSE = sprintf('%.4f', ref$rmse_mean[1]), T = good$group_iterations[1],
      Matrix = ref$fit_coordinates_mean[1], GD = good$fit_coordinates_mean[1],
      Ratio = sprintf('%.2f', good$fit_coordinates_mean[1] / ref$fit_coordinates_mean[1]))
  }
  if (length(values)) tables$simulation_communication <- do.call(rbind, values)
  if (input == 'reference') {
    expected <- read_table('reproduction/reference/simulation/method-summary.csv')
    keys <- c('case_id', 'method', 'row_type', 'R', 'group_lambda', 'group_iterations', 'group_budget')
    key <- function(a) apply(a[, keys], 1, paste, collapse = '|')
    stopifnot(setequal(key(s), key(expected)))
    expected <- expected[match(key(s), key(expected)), ]
    stopifnot(isTRUE(all.equal(s$rmse_mean, expected$rmse_mean, tolerance = 1e-10)),
              isTRUE(all.equal(s$rmse_mcse, expected$rmse_mcse, tolerance = 1e-10)),
              identical(s$successful, expected$successful))
  }
}
if (!is.null(natural)) {
  p <- subset(natural, task %in% c('target_Brier', 'target_MSE'))
  bytarget <- aggregate(value ~ dataset + target + method + task, p, mean)
  overall <- aggregate(value ~ dataset + method + task, bytarget, mean)
  overall$value[overall$task == 'target_MSE'] <- sqrt(overall$value[overall$task == 'target_MSE'])
  write.csv(overall, file.path(table_dir, 'natural-summary.csv'), row.names = FALSE)
  values <- lapply(intersect(c('heart', 'beijing'), unique(overall$dataset)), function(dataset) {
    z <- overall[overall$dataset == dataset, ]; digits <- if (dataset == 'heart') 4 else 2
    methods <- c('Unweighted', 'LocalKernel', 'Proposed', 'PooledRF', 'PooledKernel')
    v <- setNames(lapply(methods, function(m) formatC(z$value[z$method == m], format = 'f', digits = digits)),
                  c('Unweighted', 'Local', 'GD', 'RF', 'Kernel'))
    data.frame(Study = if (dataset == 'heart') 'Heart disease' else 'Beijing air quality',
               Metric = if (dataset == 'heart') 'Brier' else 'RMSE', v)
  })
  tables$applications_natural <- do.call(rbind, values)
  if (input == 'reference') {
    expected <- read_table('reproduction/reference/natural/overall-prediction.csv')
    ii <- match(paste(overall$dataset, overall$method), paste(expected$dataset, expected$method))
    stopifnot(!anyNA(ii), isTRUE(all.equal(overall$value, expected$value[ii], tolerance = 1e-10)))
  }
}
if (!is.null(spam)) {
  summary <- aggregate(value ~ allocation + method + metric, spam, mean)
  write.csv(summary, file.path(table_dir, 'spam-summary.csv'), row.names = FALSE)
  tables$applications_spam <- do.call(rbind, lapply(c('Random', 'Heterogeneous'), function(allocation) {
    z <- summary[summary$allocation == allocation & summary$metric == 'AUPRC', , drop = FALSE]
    v <- setNames(lapply(c('PositiveKDE', 'LocalKernel', 'Proposed', 'PooledRF', 'PooledKernel'),
      function(m) sprintf('%.4f', z$value[z$method == m])), c('KDE', 'Local', 'GD', 'RF', 'Kernel'))
    data.frame(Allocation = allocation, v)
  }))
  stopifnot(nrow(tables$applications_spam) == 2L)
  if (input == 'reference') {
    expected <- read_table('reproduction/reference/spam/summary.csv')
    ii <- match(paste(summary$allocation, summary$method, summary$metric),
                paste(expected$allocation, expected$method, expected$metric))
    stopifnot(!anyNA(ii), isTRUE(all.equal(summary$value, expected$statistic[ii], tolerance = 1e-10)))
  }
}
lines <- c('# Reconstructed paper tables', '',
           paste('Input:', input),
           'Simulation entries are mean RMSE (MCSE); application summaries use the manuscript averaging rules.',
           'Reference mode performs no fitting. Smoke outputs are not paper results.', '')
escape <- function(x) gsub('_', '\\_', x, fixed = TRUE)
for (name in names(tables)) {
  tab <- tables[[name]]
  write.csv(tab, file.path(table_dir, paste0(name, '.csv')), row.names = FALSE)
  latex <- c(paste0('\\begin{tabular}{', paste(rep('l', ncol(tab)), collapse = ''), '}'),
    '\\toprule', paste0(paste(escape(names(tab)), collapse = ' & '), ' \\\\'), '\\midrule',
    apply(tab, 1, function(x) paste0(paste(escape(x), collapse = ' & '), ' \\\\')),
    '\\bottomrule', '\\end{tabular}')
  writeLines(latex, file.path(table_dir, paste0(name, '.tex')))
  lines <- c(lines, paste0('## ', name), '', paste('|', paste(names(tab), collapse = ' | '), '|'),
    paste('|', paste(rep('---', ncol(tab)), collapse = ' | '), '|'),
    apply(tab, 1, function(x) paste('|', paste(x, collapse = ' | '), '|')), '')
}
writeLines(lines, file.path(out, 'PAPER_TABLES.md'))
writeLines(capture.output(sessionInfo()), file.path(out, 'analysis-session.txt'))
cat('Reconstructed', length(tables), 'tables at', out, '\n')
