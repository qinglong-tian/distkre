for(f in list.files('R',full.names=TRUE,pattern='[.]R$')) source(f)
source('applications/R/data.R'); source('applications/R/protocol.R'); source('applications/R/expansion.R')
source('applications/local-comparison/fixtures.R'); source('applications/local-comparison/methods.R')
source('applications/local-comparison/pu.R')
d<-app_dataset('spambase'); p<-exp_partition(d,31)
changed<-d; changed$y[p$evaluation]<-1-changed$y[p$evaluation]
config<-list(features=32L,iterations=1000L,lambdas=.01,bandwidths=1)
a<-lc_compare_pu(lc_pu_fixture(d,31,180L),config)
b<-lc_compare_pu(lc_pu_fixture(changed,31,180L),config)
stopifnot(identical(a$selected_parameters,b$selected_parameters),
          identical(a$diagnostics,b$diagnostics),!anyDuplicated(a$manifest$id),
          isTRUE(a$partition_invariance))
cat('PU target/evaluation-label isolation and client-partition checks passed.\n')
