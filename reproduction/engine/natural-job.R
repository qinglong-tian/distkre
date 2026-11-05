# Extracted unchanged job computation from applications/run-natural-domains.R.
run<-function(job) {
  key<-paste(job$name,job$target,job$seed,sep='-'); path<-file.path(out,paste0(key,'.rds'))
  if(file.exists(path)) return(NULL)
  cat(key,'\n'); flush.console()
  result<-tryCatch({
    d<-data[[job$name]]; sp<-natural_split(d,job$target,job$seed,config$cap)
    prep<-app_preprocess(list(d$x[sp$q,,drop=FALSE]))
    tx<-function(ii) app_transform(d$x[ii,,drop=FALSE],prep)
    xs<-lapply(sp$fit,tx); xv<-lapply(sp$val,tx); q<-tx(sp$q); qv<-tx(sp$qval)
    fitted<-lc_fit_global(xs,xv,q,qv,config,job$seed+7000L)
    local<-lc_fit_locals(xs,xv,q,qv,config,job$seed+7000L,'kernel')
    fitted$Unweighted<-list(predict=function(x) rep(1,nrow(x)))
    xm<-lapply(sp$model,tx); ym<-lapply(sp$model,function(ii) d$y[ii])
    common<-app_predictor(xm,ym,lapply(ym,function(y) rep(1,length(y))),d$type)
    xe<-lapply(sp$evaluation,tx); ye<-lapply(sp$evaluation,function(ii) d$y[ii])
    test<-tx(sp$test); ytest<-d$y[sp$test]
    loss<-function(y,p) if(d$type=='binary') (y-p)^2 else (y-p)^2
    truth<-mean(loss(ytest,common$predict(test)))
    evloss<-Map(function(x,y) loss(y,common$predict(x)),xe,ye)
    pi<-lengths(sp$fit)/sum(lengths(sp$fit)); rows<-list()
    add<-function(method,task,value,reference=NA_real_) {
      rows[[length(rows)+1L]]<<-data.frame(dataset=job$name,target=job$target,seed=job$seed,
        method=method,task=task,value=value,reference=reference,
        error=if(is.finite(reference)) value-reference else NA_real_)
    }
    for(method in c(names(fitted),'LocalKernel')) {
      funs<-if(method=='LocalKernel') lapply(local,`[[`,'predict') else
        rep(list(fitted[[method]]$predict),length(xs))
      # Primary audit: same fixed model for every ratio comparison.
      if(method=='LocalKernel') {
        estimates<-vapply(seq_along(xe),function(c) {
          w<-pmax(0,funs[[c]](xe[[c]])); if(sum(w)<=0) stop('Zero local weight mass')
          sum(w*evloss[[c]])/sum(w)
        },numeric(1))
        estimate<-sum(pi*estimates)
      } else {
        weights<-Map(function(f,x) pmax(0,f(x)),funs,xe)
        estimate<-sum(unlist(weights)*unlist(evloss))/sum(unlist(weights))
      }
      # Risk ratio is dimensionless, so regression units do not dominate summaries.
      add(method,'relative_risk_estimate',estimate/truth,1)
      add(method,'target_risk_estimate',estimate,truth)
      # Secondary: genuine target prediction using observed source outcomes.
      if(method=='LocalKernel') {
        each<-lapply(seq_along(xm),function(c) app_predictor(list(xm[[c]]),list(ym[[c]]),
          list(pmax(0,funs[[c]](xm[[c]]))),d$type))
        targetloss<-sum(pi*vapply(each,function(f) mean(loss(ytest,f$predict(test))),numeric(1)))
      } else {
        w<-Map(function(f,x) pmax(0,f(x)),funs,xm)
        f<-app_predictor(xm,ym,w,d$type)
        targetloss<-mean(loss(ytest,f$predict(test)))
      }
      add(method,if(d$type=='binary') 'target_Brier' else 'target_MSE',targetloss)
    }
    list(metrics=do.call(rbind,rows),manifest=sp$manifest,
      parameters=list(global=lapply(fitted,function(f) f[c('lambda','bandwidth')]),
                      local=lapply(local,function(f) f[c('lambda','bandwidth')])),
      counts=c(source_fit=sum(lengths(sp$fit)),target_fit=length(sp$q),target_test=length(sp$test)),
      natural_sources=sp$sources,natural_target=job$target,
      caveat='Observed domains/outcomes only. No constructed policies, label masking, contamination or covariate reweighting. Conditional invariance is not established. Secondary local prediction uses separate local models, unlike the common-predictor audit.')
  },error=function(e) list(error=conditionMessage(e)))
  atomic_rds(result,path); if(!is.null(result$error)) cat('FAILED',result$error,'\n'); NULL
}
