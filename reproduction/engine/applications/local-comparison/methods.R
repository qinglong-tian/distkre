# All methods receive the same representation and target-only standardization.
# Local methods receive one client's source data and the common target sample.
# Evaluation outcomes are never accepted by the fitting interfaces.
lc_function <- function(fit) local({ f<-fit; function(x) as.numeric(predict(f,x)) })
lc_fit_global <- function(xs,xv,q,qv,config,seed) {
  bandwidth<-median_pairwise_distance(q)
  stopifnot(is.finite(bandwidth),bandwidth>0)
  nv<-vapply(xv,nrow,integer(1)); wv<-nv/sum(nv)
  validation<-function(fun) sum(wv*vapply(xv,function(x) mean(fun(x)^2)/2,numeric(1)))-mean(fun(qv))
  grid<-expand.grid(lambda=config$lambdas,mult=config$bandwidths)
  output<-list()
  for(method in c('Proposed','PooledRF','PooledKernel')) {
    best<-Inf; costs<-0
    for(i in seq_len(nrow(grid))) {
      bw<-bandwidth*grid$mult[i]; lambda<-grid$lambda[i]
      f<-if(method=='PooledKernel') fit_kulsif(do.call(rbind,xs),q,lambda,bw) else
        dkre_rff(xs,q,lambda,config$features,bw,method=if(method=='Proposed') 'gd' else 'closed',
                 iterations=config$iterations,seed=seed)
      fun<-lc_function(f); loss<-validation(fun)
      if(method=='Proposed') costs<-costs+f$metadata$communication_coordinates+
        length(xs)*(config$features+1)
      if(loss<best) {best<-loss; chosen<-list(predict=fun,fit=f,
        validation=loss,lambda=lambda,bandwidth=bw)}
    }
    chosen$ratio_tuning_coordinates<-costs
    output[[method]]<-chosen
  }
  output
}

lc_fit_one_client <- function(x,xv,q,qv,config,seed,method=c('kernel','rf')) {
  method<-match.arg(method); med<-median_pairwise_distance(q)
  grid<-expand.grid(lambda=config$lambdas,mult=config$bandwidths)
  best<-Inf
  for(i in seq_len(nrow(grid))) {
    bw<-med*grid$mult[i]; lambda<-grid$lambda[i]
    f<-if(method=='kernel') fit_kulsif(x,q,lambda,bw) else
      dkre_rff(list(x),q,lambda,config$features,bw,method='closed',seed=seed)
    fun<-lc_function(f); loss<-mean(fun(xv)^2)/2-mean(fun(qv))
    if(loss<best) {best<-loss; selected<-list(predict=fun,validation=loss,
      lambda=lambda,bandwidth=bw,source_n=nrow(x),ratio_source_coordinates=0L)}
  }
  selected
}

lc_fit_locals <- function(xs,xv,q,qv,config,seed,method) {
  # Separate calls: no other client's rows, moments, coefficients or validation
  # losses are inputs. Shared map generation information is public, not learned.
  Map(function(x,v) lc_fit_one_client(x,v,q,qv,config,seed,method),xs,xv)
}
lc_site_prediction <- function(models,x,owner) {
  out<-rep(NA_real_,nrow(x))
  for(c in 1:3) { ii<-which(owner==c); if(length(ii)) out[ii]<-models[[c]]$predict(x[ii,,drop=FALSE]) }
  stopifnot(all(is.finite(out))); out
}
lc_local_estimates <- function(models,x,values,owner,pi) {
  site<-vapply(1:3,function(c) {
    ii<-which(owner==c); if(!length(ii)) return(NA_real_)
    w<-pmax(0,models[[c]]$predict(x[ii,,drop=FALSE])); mass<-sum(w)
    if(mass<=0) return(NA_real_)
    sum(w*values[ii])/mass
  },numeric(1))
  list(value=if(all(is.finite(site))) sum(pi*site) else NA_real_,site=site)
}

lc_compare <- function(fixture,config) {
  # Public target moments provide an identical coordinate system without any
  # source-to-source preprocessing exchange in the local-only baseline.
  prep<-app_preprocess(list(fixture$target))
  tx<-function(x) app_transform(x,prep)
  x<-tx(fixture$train); v<-tx(fixture$validation)
  q<-tx(fixture$target); qv<-tx(fixture$target_validation); xe<-tx(fixture$evaluation)
  metric_rows<-list(); diagnostic_rows<-list(); selections<-list(); invariance<-list()
  for(allocation in names(fixture$allocations)) {
    own<-fixture$allocations[[allocation]]
    stopifnot(all(own$train %in% 1:3),all(own$validation %in% 1:3),all(own$evaluation %in% 1:3))
    xs<-lapply(1:3,function(c) x[own$train==c,,drop=FALSE])
    xv<-lapply(1:3,function(c) v[own$validation==c,,drop=FALSE])
    if(min(vapply(xs,nrow,integer(1)))<5 || min(vapply(xv,nrow,integer(1)))<2)
      stop('Client lacks fitting or validation observations')
    pi<-vapply(xs,nrow,integer(1)); pi<-pi/sum(pi)
    global<-lc_fit_global(xs,xv,q,qv,config,fixture$seed+7000L)
    localK<-lc_fit_locals(xs,xv,q,qv,config,fixture$seed+7000L,'kernel')
    localRF<-lc_fit_locals(xs,xv,q,qv,config,fixture$seed+7000L,'rf')
    selections[[allocation]]<-list(global=lapply(global,function(f) f[c('lambda','bandwidth','validation')]),
      LocalKernel=lapply(localK,function(f) f[c('lambda','bandwidth','validation','source_n')]),
      LocalRF=lapply(localRF,function(f) f[c('lambda','bandwidth','validation','source_n')]))
    invariance[[allocation]]<-list(proposed=global$Proposed$predict(xe),
      parameters=c(global$Proposed$lambda,global$Proposed$bandwidth))
    methods<-c(names(global),'LocalKernel','LocalRF',if(fixture$case=='ecology') 'PositiveKDE' else 'Unweighted')
    if(fixture$case=='ecology') {
      best<- -Inf; med<-median_pairwise_distance(q)
      for(mult in c(.5,1,2)) {
        density<-.density_provider(q,med*mult); val<-mean(log(pmax(1e-300,density(qv))))
        if(val>best) {best<-val; positive_density<-density}
      }
    }
    add<-function(method,metric,value,estimate=NA_real_,truth=NA_real_) {
      metric_rows[[length(metric_rows)+1L]]<<-data.frame(case=fixture$case,dataset=fixture$dataset,
        scenario=fixture$scenario,replication=fixture$seed,allocation=allocation,method=method,
        metric=metric,value=value,estimate=estimate,truth=truth)
    }
    for(method in methods) {
      local<-method %in% c('LocalKernel','LocalRF')
      if(local) {
        models<-if(method=='LocalKernel') localK else localRF
        raw<-lc_site_prediction(models,xe,own$evaluation)
      } else if(method=='Unweighted') raw<-rep(1,nrow(xe)) else if(method=='PositiveKDE')
        raw<-positive_density(xe) else raw<-global[[method]]$predict(xe)
      if(fixture$case=='ecology') {
        pooledAUC<-app_rank_metrics(fixture$values,raw)['AUROC']
        add(method,'global_AUC',unname(pooledAUC))
        within<-vapply(1:3,function(c) {
          ii<-which(own$evaluation==c)
          if(length(unique(fixture$values[ii]))<2) return(NA_real_)
          unname(app_rank_metrics(fixture$values[ii],raw[ii])['AUROC'])
        },numeric(1))
        # Some regions contain no held-out presence. Save the evaluable count;
        # do not invent an AUC for a one-class region.
        ok<-is.finite(within)
        add(method,'within_client_AUC',if(any(ok)) sum(pi[ok]*within[ok])/sum(pi[ok]) else NA_real_)
        add(method,'evaluable_clients',sum(ok))
      } else {
        w<-pmax(0,raw)
        if(local) {
          a<-lc_local_estimates(models,xe,fixture$values,own$evaluation,pi)
          estimate<-a$value
          add(method,'individual_client_MSE',sum(pi*(a$site-fixture$truth)^2))
        } else estimate<-if(sum(w)>0) sum(w*fixture$values)/sum(w) else NA_real_
        if(!is.finite(estimate)) stop('Undefined downstream estimate: ',method)
        add(method,'squared_error',(estimate-fixture$truth)^2,estimate,fixture$truth)
        add(method,'absolute_error',abs(estimate-fixture$truth),estimate,fixture$truth)
      }
      diag<-data.frame(case=fixture$case,dataset=fixture$dataset,scenario=fixture$scenario,
        replication=fixture$seed,allocation=allocation,method=method,
        negative_fraction=mean(raw<0),max_prediction=max(raw),
        source_fit_n=sum(vapply(xs,nrow,integer(1))),target_fit_n=nrow(q),
        fitting_coordinates=if(method=='Proposed') global[[method]]$fit$metadata$communication_coordinates else
          if(local) 0 else NA_real_,
        tuning_coordinates=if(method=='Proposed') global[[method]]$ratio_tuning_coordinates else if(local) 0 else NA_real_)
      diagnostic_rows[[length(diagnostic_rows)+1L]]<-diag
    }
  }
  # Reallocation must not change the pooled task, validation rule, or fit.
  first<-invariance[[1]]
  for(a in invariance[-1]) stopifnot(identical(first$parameters,a$parameters),
                                     max(abs(first$proposed-a$proposed))<1e-7)
  list(metrics=do.call(rbind,metric_rows),diagnostics=do.call(rbind,diagnostic_rows),
       selected_parameters=selections,manifest=fixture$manifest,
       allocations=fixture$allocations,preprocessing=prep,
       qualification=fixture$qualification,partition_invariance=TRUE)
}
