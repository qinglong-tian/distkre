# Fixed observations and role splits are constructed before client allocation.
# No target evaluation outcome enters a ratio fit or parameter selection.
lc_partition <- function(n) {
  ii<-sample.int(n); k<-floor(.8*n)
  list(train=ii[seq_len(k)],validation=ii[-seq_len(k)])
}
lc_random_owner <- function(n) sample(rep(1:3,length.out=n))
lc_manifest <- function(ids,roles) {
  stopifnot(!anyDuplicated(ids)); data.frame(id=ids,role=roles)
}

lc_audit_fixture <- function(d,seed,scenario,n=600L) {
  p<-exp_partition(d,seed); model<-exp_policy_model(d,p$policy)
  propensity<-switch(scenario,
    selective_risk=function(z) .05+.9*plogis(-2*z),
    selective_confidence=function(z) .05+.9*exp(-z^2/.6),
    no_shift=function(z) rep(.5,length(z)))
  choose<-function(ii) app_take(ii[runif(length(ii))<propensity(model$score(ii))],n)
  src<-choose(p$ratio); ev<-choose(p$evaluation)
  target<-app_take(p$target_fit,n); truthids<-app_take(p$target_truth,2*n)
  sp<-lc_partition(length(src)); sq<-lc_partition(length(target))
  z<-matrix(model$score(src),ncol=1); q<-matrix(model$score(target),ncol=1)
  xe<-matrix(model$score(ev),ncol=1)
  loss<-function(ii) (model$probability(ii)-d$y[ii])^2
  strata<-function(z) as.integer(cut(z,c(-Inf,-.5,.5,Inf)))
  alloc<-list(Random=list(train=lc_random_owner(length(sp$train)),
    validation=lc_random_owner(length(sp$validation)),evaluation=lc_random_owner(length(ev))),
    Heterogeneous=list(train=strata(z[sp$train,1]),validation=strata(z[sp$validation,1]),
                       evaluation=strata(xe[,1])))
  list(case='risk_audit',dataset='adult',seed=seed,scenario=scenario,
    train=z[sp$train,,drop=FALSE],validation=z[sp$validation,,drop=FALSE],
    target=q[sq$train,,drop=FALSE],target_validation=q[sq$validation,,drop=FALSE],
    evaluation=xe,values=loss(ev),truth=mean(loss(truthids)),allocations=alloc,
    manifest=lc_manifest(c(model$training_ids,d$id[c(src,ev,target,truthids)]),
      c(rep('policy',length(model$training_ids)),rep('source',length(src)),
        rep('source_evaluation',length(ev)),rep('target',length(target)),rep('target_truth',length(truthids)))),
    qualification='Real census-derived records with constructed X-only selection. Hard confidence strata may not individually cover the target; pooled support remains the intended requirement.')
}

lc_gas_fixture <- function(d,seed,n=600L) {
  set.seed(seed); p<-exp_partition(d,seed)
  development<-app_take(p$policy,3000L)
  prep<-app_preprocess(list(d$x[development,,drop=FALSE]))
  xx<-app_transform(d$x[development,,drop=FALSE],prep)
  models<-lapply(1:6,function(k) app_predictor(list(xx),
    list(as.numeric(d$y[development]==k)),list(rep(1,length(development))),'binary',rho=.01))
  policy<-function(ii) {
    x<-app_transform(d$x[ii,,drop=FALSE],prep)
    a<-sapply(models,function(f) f$predict(x)); a<-a/rowSums(a)
    .05/6+.95*a
  }
  base<-matrix(.2/6,3,6); base[cbind(1:3,c(1,3,6))]<-.8+.2/6
  adaptive<-c(.1,.3,.6)
  draw<-function(ii,target=FALSE) {
    ii<-app_take(ii,n); pp<-policy(ii); client<-lc_random_owner(length(ii))
    a<-vapply(seq_along(ii),function(i) sample.int(6,1,prob=if(target) pp[i,] else
      adaptive[client[i]]*pp[i,]+(1-adaptive[client[i]])*base[client[i],]),integer(1))
    list(ids=ii,owner=client,x=cbind(pp[cbind(seq_along(ii),a)],
      sapply(1:6,function(k) as.numeric(a==k))),reward=as.numeric(a==d$y[ii]))
  }
  src<-draw(p$ratio); ev<-draw(p$evaluation); q<-draw(p$target_fit,TRUE)
  sp<-lc_partition(length(src$ids)); sq<-lc_partition(length(q$ids))
  truthids<-app_take(p$target_truth,2*n); pp<-policy(truthids)
  alloc<-list(Random=list(train=lc_random_owner(length(sp$train)),
    validation=lc_random_owner(length(sp$validation)),evaluation=lc_random_owner(length(ev$ids))),
    Heterogeneous=list(train=src$owner[sp$train],validation=src$owner[sp$validation],evaluation=ev$owner))
  list(case='policy_evaluation',dataset='gas',seed=seed,scenario='six_actions',
    train=src$x[sp$train,,drop=FALSE],validation=src$x[sp$validation,,drop=FALSE],
    target=q$x[sq$train,,drop=FALSE],target_validation=q$x[sq$validation,,drop=FALSE],
    evaluation=ev$x,values=ev$reward,truth=mean(pp[cbind(seq_along(truthids),d$y[truthids])]),
    allocations=alloc,
    manifest=lc_manifest(d$id[c(development,src$ids,ev$ids,q$ids,truthids)],
      c(rep('policy',length(development)),rep('source',length(src$ids)),
        rep('source_evaluation',length(ev$ids)),rep('target',length(q$ids)),rep('target_truth',length(truthids)))),
    qualification='Actual six-gas measurements with simulated context-dependent logging policies. Every logger assigns positive probability to every action; local transport is identifiable in the population.')
}

lc_ecology_fixture <- function(d,species,fold) {
  set.seed(2000+fold)
  columns<-setdiff(names(d),c('species','fold','role','cell','lon','lat'))
  x<-as.matrix(d[,columns]); background<-which(d$role=='background')
  pos<-which(d$species==species & d$fold==fold & d$role=='train')
  test<-which(d$species==species & d$fold==fold & d$role=='test')
  pos<-pos[!d$cell[pos] %in% d$cell[test]]
  pos<-sample(pos[!duplicated(d$cell[pos])]); test<-test[!duplicated(d$cell[test])]
  bg<-sample(background); src<-bg[1:600]; testbg<-bg[601:2600]
  sp<-lc_partition(length(src)); sq<-lc_partition(length(pos))
  # Public regional boundaries are determined from the fixed source geography;
  # they define ownership, rather than fitted preprocessing or tuned parameters.
  breaks<-quantile(d$lon[src],c(1/3,2/3))
  owner<-function(ii) as.integer(cut(d$lon[ii],c(-Inf,breaks,Inf)))
  eval<-c(test,testbg)
  alloc<-list(Random=list(train=lc_random_owner(length(sp$train)),
    validation=lc_random_owner(length(sp$validation)),evaluation=lc_random_owner(length(eval))),
    Heterogeneous=list(train=owner(src[sp$train]),validation=owner(src[sp$validation]),evaluation=owner(eval)))
  list(case='ecology',dataset=species,seed=fold,scenario='published_fold',
    train=x[src[sp$train],,drop=FALSE],validation=x[src[sp$validation],,drop=FALSE],
    target=x[pos[sq$train],,drop=FALSE],target_validation=x[pos[sq$validation],,drop=FALSE],
    evaluation=x[eval,,drop=FALSE],values=c(rep(1,length(test)),rep(0,length(testbg))),
    truth=NA_real_,allocations=alloc,
    manifest=lc_manifest(as.character(d$cell[c(src,pos,eval)]),
      c(rep('source',length(src)),rep('positive_train',length(pos)),
        rep('positive_test',length(test)),rep('background_test',length(testbg)))),
    qualification='Actual occurrences and environmental grids; geographic clients are constructed longitude bands. Background cells are not verified absences. Published folds reuse occurrences across folds.')
}
