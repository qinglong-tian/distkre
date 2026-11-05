simulation_config <- function() list(
  version = 'three-studies-v1', workers = 10L, dimension = 5L,
  features = 200L, test_points = 20000L, batch_size = 500L,
  lambdas = c(1e-4, 1e-3, 1e-2, .1, 1), bandwidth_multipliers = c(.5,1,2),
  tuning_iterations = 5000L, final_iterations = 10000L,
  mechanism_features = c(50L,200L,400L), mechanism_lambdas = c(.001,.01,.1),
  checkpoints = c(0L,1L,2L,5L,10L,20L,50L,100L,200L,500L,1000L,2000L,5000L,10000L,20000L),
  communication_budgets = c(1e5,4e5,1.6e6,6.4e6),
  rho = c(.05,.15,.30,.50), pilot_replications = 20L,
  production_tuned = 200L, production_mechanism = 100L,
  seed = 20261102L, exact_kernel_cap = 2000L
)
simulation_cases <- function(config) {
  sizes <- unique(rbind(data.frame(N=c(1000,2000,4000,8000),M=1000),
                        data.frame(N=2000,M=c(250,500,1000,2000))))
  cases <- lapply(seq_len(nrow(sizes)), function(i) list(
    id=sprintf('accuracy_N%d_M%d',sizes$N[i],sizes$M[i]),study='accuracy',
    N=sizes$N[i],M=sizes$M[i],tau=1,pi=rep(.25,4),protocol='tuned',
    replications=config$production_tuned))
  cases[[length(cases)+1]] <- list(id='imbalance',study='imbalance',N=2000,M=1000,
    tau=1,pi=c(.7,.1,.1,.1),protocol='tuned',replications=config$production_tuned)
  for(tau in c(0,.5,1)) cases[[length(cases)+1]] <- list(
    id=paste0('heterogeneity_',tau),study='heterogeneity',N=2000,M=1000,
    tau=tau,pi=rep(.25,4),protocol='fixed',replications=config$production_mechanism)
  cases[[length(cases)+1]] <- list(id='communication',study='communication',
    N=2000,M=1000,tau=1,pi=rep(.25,4),protocol='path',replications=config$production_mechanism)
  for(i in seq_along(cases)) cases[[i]]$index <- i
  cases
}
