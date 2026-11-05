# Positive-unlabeled spam ranking: verified positive reference at the server,
# unlabeled source mail at clients. Only the final evaluator sees test labels.
lc_pu_fixture <- function(d,seed,n=600L) {
  set.seed(seed); p<-exp_partition(d,seed)
  source<-app_take(p$ratio,n); evaluation<-app_take(p$evaluation,2*n)
  target<-app_take(p$target_fit[d$y[p$target_fit]==1],n)
  stopifnot(length(target)>=30)
  sp<-lc_partition(length(source)); sq<-lc_partition(length(target))
  # Constructed provider partitions use capitalization style, never outcomes.
  owner<-function(ii) vapply(log1p(d$x[ii,55]),function(v) {
    t<-tanh(v-log(3)); weight<-exp(c(-3*t,0,3*t))
    # Full overlap: every record has at least 3% probability at every client.
    sample.int(3,1,prob=.03+.91*weight/sum(weight))
  },integer(1))
  alloc<-list(Random=list(train=lc_random_owner(length(sp$train)),
    validation=lc_random_owner(length(sp$validation)),evaluation=lc_random_owner(length(evaluation))),
    Heterogeneous=list(train=owner(source[sp$train]),validation=owner(source[sp$validation]),evaluation=owner(evaluation)))
  list(case='positive_unlabeled',dataset='spambase',scenario='case_only_spam',seed=seed,
    train=d$x[source[sp$train],,drop=FALSE],validation=d$x[source[sp$validation],,drop=FALSE],
    target=d$x[target[sq$train],,drop=FALSE],target_validation=d$x[target[sq$validation],,drop=FALSE],
    evaluation=d$x[evaluation,,drop=FALSE],values=d$y[evaluation],allocations=alloc,
    manifest=lc_manifest(d$id[c(source,target,evaluation)],
      c(rep('source',length(source)),rep('positive_reference',length(target)),rep('evaluation',length(evaluation)))),
    qualification='Real email features with constructed, probabilistic provider allocation and random positive-only labeling. Every provider has full population support. Representative verified positives are assumed; source/test class labels are unavailable to fitting.')
}
lc_compare_pu <- function(fixture,config) {
  prep<-app_preprocess(list(fixture$target)); tx<-function(x) app_transform(x,prep)
  x<-tx(fixture$train); v<-tx(fixture$validation); q<-tx(fixture$target)
  qv<-tx(fixture$target_validation); xe<-tx(fixture$evaluation)
  metrics<-diagnostics<-selected<-list(); predictions<-list()
  for(allocation in names(fixture$allocations)) {
    a<-fixture$allocations[[allocation]]
    xs<-lapply(1:3,function(c) x[a$train==c,,drop=FALSE])
    xv<-lapply(1:3,function(c) v[a$validation==c,,drop=FALSE])
    stopifnot(min(vapply(xs,nrow,integer(1)))>=5,min(vapply(xv,nrow,integer(1)))>=2)
    global<-lc_fit_global(xs,xv,q,qv,config,fixture$seed+7000)
    local<-lc_fit_locals(xs,xv,q,qv,config,fixture$seed+7000,'kernel')
    localrf<-lc_fit_locals(xs,xv,q,qv,config,fixture$seed+7000,'rf')
    best<- -Inf; med<-median_pairwise_distance(q)
    for(mult in c(.5,1,2)) {
      density<-.density_provider(q,med*mult); val<-mean(log(pmax(1e-300,density(qv))))
      if(val>best) {best<-val; density_selected<-density}
    }
    selected[[allocation]]<-list(global=lapply(global,function(f) f[c('lambda','bandwidth','validation')]),
      LocalKernel=lapply(local,function(f) f[c('lambda','bandwidth','validation')]),
      LocalRF=lapply(localrf,function(f) f[c('lambda','bandwidth','validation')]))
    for(method in c(names(global),'LocalKernel','LocalRF','PositiveKDE')) {
      raw<-switch(method,LocalKernel=lc_site_prediction(local,xe,a$evaluation),
        LocalRF=lc_site_prediction(localrf,xe,a$evaluation),PositiveKDE=density_selected(xe),
        global[[method]]$predict(xe))
      score<-app_rank_metrics(fixture$values,raw)
      for(metric in names(score)) metrics[[length(metrics)+1L]]<-data.frame(
        case=fixture$case,dataset=fixture$dataset,scenario=fixture$scenario,replication=fixture$seed,
        allocation=allocation,method=method,metric=metric,value=unname(score[metric]),estimate=NA_real_,truth=NA_real_)
      pi<-vapply(xs,nrow,integer(1)); pi<-pi/sum(pi)
      for(metric in names(score)) {
        local_score<-vapply(1:3,function(c) {
          ii<-which(a$evaluation==c)
          if(length(unique(fixture$values[ii]))<2) return(NA_real_)
          unname(app_rank_metrics(fixture$values[ii],raw[ii])[metric])
        },numeric(1))
        ok<-is.finite(local_score)
        metrics[[length(metrics)+1L]]<-data.frame(case=fixture$case,dataset=fixture$dataset,
          scenario=fixture$scenario,replication=fixture$seed,allocation=allocation,method=method,
          metric=paste0('within_client_',metric),value=sum(pi[ok]*local_score[ok])/sum(pi[ok]),
          estimate=NA_real_,truth=NA_real_)
      }
      diagnostics[[length(diagnostics)+1L]]<-data.frame(allocation=allocation,method=method,
        negative_fraction=mean(raw<0),max_prediction=max(raw))
      if(method=='Proposed') predictions[[allocation]]<-raw
    }
  }
  stopifnot(max(abs(predictions[[1]]-predictions[[2]]))<1e-7)
  list(metrics=do.call(rbind,metrics),diagnostics=do.call(rbind,diagnostics),
    manifest=fixture$manifest,allocations=fixture$allocations,selected_parameters=selected,
    partition_invariance=TRUE,qualification=fixture$qualification)
}
