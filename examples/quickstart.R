library(distkre)

set.seed(42)
sources <- list(matrix(rnorm(200, -0.6), ncol = 2),
                matrix(rnorm(300, 0.6), ncol = 2))
target <- matrix(rnorm(240, 0.2), ncol = 2)

fit <- dkre_rff(sources, target, lambda = 0.1, n_features = 64,
                method = 'gd', iterations = 1000, seed = 42)
print(fit)
print(head(predict(fit, target)))

# Reuse the model after saving; raw observations are omitted by default.
path <- tempfile(fileext = '.rds')
saveRDS(fit, path)
restored <- readRDS(path)
stopifnot(identical(predict(fit, target), predict(restored, target)))
unlink(path)
