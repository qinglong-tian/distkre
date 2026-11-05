# Reaggregate completed application jobs without refitting.
source('reproduction/common.R')
output <- absolute(arg('input', 'reproduction/output/full/applications'))
stopifnot(file.exists(file.path(output, 'configuration.rds')))
for (folder in intersect(c('natural', 'spam'), list.dirs(output, recursive = FALSE, full.names = FALSE))) {
  out <- file.path(output, folder)
  files <- list.files(out, '[.]rds$', full.names = TRUE)
  results <- lapply(files, readRDS)
  stopifnot(all(vapply(results, function(x) is.null(x$error), logical(1))))
  for (x in results) {
    stopifnot(!anyDuplicated(x$manifest$id), all(is.finite(x$metrics$value)))
    if (folder == 'natural') {
      used <- subset(x$manifest, role != 'unused')
      target <- used$role %in% c('q', 'qval', 'test')
      stopifnot(all(used$site[target] == x$natural_target),
                all(used$site[!target] %in% x$natural_sources))
      if (x$metrics$dataset[1] == 'beijing') {
        role <- used$role
        role[role == 'q'] <- 'fit'; role[role == 'qval'] <- 'val'; role[role == 'test'] <- 'evaluation'
        stopifnot(all(vapply(split(role, used$block), function(a) length(unique(a)) == 1L, logical(1))))
      }
    } else stopifnot(isTRUE(x$partition_invariance))
  }
  write.csv(do.call(rbind, lapply(results, `[[`, 'metrics')), file.path(out, 'metrics.csv'), row.names = FALSE)
}
cat('Application jobs audited and aggregated:', output, '\n')
