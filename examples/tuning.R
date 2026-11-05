library(distkre)
set.seed(51)
train <- list(matrix(rnorm(120, -0.4), ncol = 2), matrix(rnorm(160, 0.5), ncol = 2))
validation <- list(matrix(rnorm(60, -0.4), ncol = 2), matrix(rnorm(80, 0.5), ncol = 2))
target <- matrix(rnorm(120, 0.1), ncol = 2)
target_validation <- matrix(rnorm(60, 0.1), ncol = 2)
pi <- vapply(train, nrow, integer(1)); pi <- pi / sum(pi)
scale <- median_pairwise_distance(target)
grid <- expand.grid(lambda = c(.01, .1), bandwidth = scale * c(.5, 1))
models <- vector('list', nrow(grid)); losses <- numeric(nrow(grid))
for (i in seq_len(nrow(grid))) {
  fit <- dkre_rff(train, target, lambda = grid$lambda[i], n_features = 64,
                  bandwidth = grid$bandwidth[i], method = 'gd', iterations = 1000, seed = 42)
  losses[i] <- sum(pi * vapply(validation, function(x) mean(predict(fit, x)^2), numeric(1))) / 2 -
    mean(predict(fit, target_validation))
  models[[i]] <- fit
}
best <- which.min(losses)
print(cbind(grid, validation_loss = losses))
print(models[[best]])
