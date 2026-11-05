# Additional application families; development and confirmation use distinct seeds.
exp_dataset <- function(name) {
  if(!name %in% c('adult','magic')) return(app_dataset(name))
  if(name=='adult') {
    path <- 'applications/data/raw/adult'
    a <- read.csv(file.path(path,'adult.data'),header=FALSE,strip.white=TRUE)
    b <- read.csv(file.path(path,'adult.test'),header=FALSE,strip.white=TRUE,skip=1)
    d <- rbind(a,b)
    x <- cbind(as.matrix(d[,c(1,5,11,12,13)]),male=as.numeric(d[,10]=='Male'),
               married=as.numeric(d[,6]=='Married-civ-spouse'))
    x[,c(3,4)] <- log1p(x[,c(3,4)])
    y <- as.numeric(gsub('[.]$','',d[,15])=='>50K')
  } else {
    d <- read.csv('applications/data/raw/magic/magic04.data',header=FALSE)
    x <- as.matrix(d[,1:10]); y <- as.numeric(d[,11]=='g')
  }
  list(x=x,y=y,name=name,type='binary',id=paste0(name,':',seq_along(y)),
       block=apply(cbind(x,y),1,paste,collapse='|'),site=rep('constructed',length(y)))
}

exp_partition <- function(d,seed) {
  set.seed(seed)
  blocks <- sample(unique(d$block))
  r <- cut(seq_along(blocks)/length(blocks),c(0,.2,.4,.65,.82,1),
           labels=c('policy','ratio','evaluation','target_fit','target_truth'))
  roles <- as.character(r[match(d$block,blocks)])
  split(seq_len(nrow(d$x)),roles)
}

exp_policy_model <- function(d,ii) {
  ii <- app_take(ii,3000L)
  prep <- app_preprocess(list(d$x[ii,,drop=FALSE]))
  fit <- app_predictor(list(app_transform(d$x[ii,,drop=FALSE],prep)),
                       list(d$y[ii]),list(rep(1,length(ii))),'binary')
  raw_score <- function(jj) qlogis(pmin(.999,pmax(.001,
                          fit$predict(app_transform(d$x[jj,,drop=FALSE],prep)))))
  mean <- mean(raw_score(ii)); sd <- sd(raw_score(ii)); if(sd<1e-8) sd <- 1
  list(score=function(jj) (raw_score(jj)-mean)/sd,
       probability=function(jj) plogis(raw_score(jj)),training_ids=d$id[ii],
       coefficients_preprocessing=prep,convergence=fit$convergence)
}

exp_clients <- function(z) {
  # Complementary confidence strata: local supports need not cover all of Q.
  split(seq_along(z),cut(z,c(-Inf,-.5,.5,Inf),labels=FALSE))
}
exp_validate <- function(clients) {
  fit <- lapply(clients,function(ii) ii[seq_len(floor(.8*length(ii)))])
  val <- Map(setdiff,clients,fit)
  if(length(fit)!=3 || min(lengths(fit))<8 || min(lengths(val))<2)
    stop('Insufficient client support/sample sizes')
  list(fit=fit,val=val)
}
exp_weighted_mean <- function(values,weights,clients,pi) {
  a <- sum(pi*vapply(clients,function(ii) mean(weights[ii]*values[ii]),numeric(1)))
  b <- sum(pi*vapply(clients,function(ii) mean(weights[ii]),numeric(1)))
  if(b<=0 || !is.finite(b)) stop('Invalid weight normalization')
  a/b
}
exp_diagnostics <- function(fits,x) {
  lapply(fits,function(f) {
    raw <- f$predict(x); w <- pmax(0,raw)
    list(validation=f$validation,negative_fraction=mean(raw<0),
         lambda=if(is.null(f$lambda)) NA_real_ else f$lambda,
         bandwidth=if(is.null(f$bandwidth)) NA_real_ else f$bandwidth,
         max_weight=max(w),ESS=if(sum(w^2)>0) sum(w)^2/sum(w^2) else 0,
         gradient=if(!is.null(f$fit)) f$fit$metadata$gradient_norm else NA_real_)
  })
}

exp_spline_ratio <- function(xs,q,xv,qv) {
  # Stronger centralized reference: action interactions and nonlinear score
  # effects prevent comparing the proposed method only with a misspecified
  # main-effects linear domain classifier.
  train <- do.call(rbind,xs)
  basis <- splines::ns(train[,1],df=4)
  feature <- function(x) {
    b <- predict(basis,x[,1])
    if(ncol(x)==1) b else cbind(b,x[,2],b*x[,2])
  }
  pi <- vapply(xs,nrow,integer(1)); pi <- pi/sum(pi)
  target_features <- feature(q); best <- Inf
  for(rho in c(.001,.01,.1)) {
    fit <- app_predictor(list(rbind(feature(train),target_features)),
      list(c(rep(0,nrow(train)),rep(1,nrow(q)))),
      list(rep(1,nrow(train)+nrow(q))),'binary',rho=rho)
    fun <- local({ f<-fit; prior<-nrow(train)/nrow(q)
      function(x) { p<-pmin(1-1e-8,pmax(1e-8,f$predict(feature(x)))); prior*p/(1-p) } })
    loss <- sum(pi*vapply(xv,function(x) mean(fun(x)^2)/2,numeric(1)))-mean(fun(qv))
    if(loss<best) {best<-loss; chosen<-list(predict=fun,validation=loss)}
  }
  chosen
}

exp_tuned_classifier <- function(xs,q,xv,qv) {
  train <- do.call(rbind,xs)
  pi <- vapply(xs,nrow,integer(1)); pi <- pi/sum(pi); best <- Inf
  for(rho in c(.001,.01,.1,1)) {
    f <- app_predictor(list(rbind(train,q)),
      list(c(rep(0,nrow(train)),rep(1,nrow(q)))),
      list(rep(1,nrow(train)+nrow(q))),'binary',rho)
    fun <- local({ fit<-f; prior<-nrow(train)/nrow(q)
      function(x) { p<-pmin(1-1e-8,pmax(1e-8,fit$predict(x))); prior*p/(1-p) } })
    loss <- sum(pi*vapply(xv,function(x) mean(fun(x)^2)/2,numeric(1)))-mean(fun(qv))
    if(loss<best) {best<-loss; chosen<-list(predict=fun,validation=loss)}
  }
  chosen
}

exp_reward_model <- function(source,sp) {
  # Direct reward-model reference, tuned with independent logged validation
  # rewards. No target outcomes or known logging propensities are used.
  train <- unlist(sp$fit); validation <- unlist(sp$val)
  basis <- splines::ns(source$x[train,1],df=4)
  feature <- function(x) { b<-predict(basis,x[,1]); cbind(b,x[,2],b*x[,2]) }
  best <- Inf
  for(rho in c(.001,.01,.1)) {
    fit <- app_predictor(lapply(sp$fit,function(i) feature(source$x[i,,drop=FALSE])),
      lapply(sp$fit,function(i) source$reward[i]),
      lapply(sp$fit,function(i) rep(1,length(i))),'regression',rho)
    fun <- local({f<-fit; function(x) pmin(1,pmax(0,f$predict(feature(x))))})
    loss <- mean((source$reward[validation]-fun(source$x[validation,,drop=FALSE]))^2)
    if(loss<best) {best<-loss; selected<-list(predict=fun,validation=loss,convergence=fit$convergence)}
  }
  selected
}

exp_audit <- function(d,seed,config,scenario='selective_risk') {
  p <- exp_partition(d,seed); model <- exp_policy_model(d,p$policy)
  # Source record retention is X-only and bounded. It represents selective
  # labeling/ascertainment; target samples represent the unfiltered population.
  propensity <- switch(scenario,
    selective_risk=function(z) .05+.9*plogis(-2*z),
    selective_confidence=function(z) .05+.9*exp(-z^2/.6),
    no_shift=function(z) rep(.5,length(z)))
  choose <- function(ii) {
    ii <- ii[runif(length(ii))<propensity(model$score(ii))]
    app_take(ii,config$sample_size)
  }
  source <- choose(p$ratio); evaluation <- choose(p$evaluation)
  target <- app_take(p$target_fit,config$sample_size)
  truth_ids <- app_take(p$target_truth,2*config$sample_size)
  zz <- model$score(source); clients <- exp_clients(zz); split <- exp_validate(clients)
  z <- matrix(zz,ncol=1); q <- matrix(model$score(target),ncol=1)
  qn <- floor(.8*nrow(q))
  fits <- app_ratio_fits(lapply(split$fit,function(ii) z[ii,,drop=FALSE]),
    q[seq_len(qn),,drop=FALSE],lapply(split$val,function(ii) z[ii,,drop=FALSE]),
    q[-seq_len(qn),,drop=FALSE],config,seed+7000L)
  fits$ClassifierSpline <- exp_spline_ratio(lapply(split$fit,function(ii) z[ii,,drop=FALSE]),
    q[seq_len(qn),,drop=FALSE],lapply(split$val,function(ii) z[ii,,drop=FALSE]),q[-seq_len(qn),,drop=FALSE])
  eval_z <- matrix(model$score(evaluation),ncol=1)
  ec <- exp_clients(as.numeric(eval_z)); pi <- lengths(split$fit)/sum(lengths(split$fit))
  if(!identical(names(ec),names(split$fit))) stop('Missing evaluation stratum')
  fits$DesignWeights <- list(predict=function(x) 1/propensity(x[,1]),validation=NA_real_)
  metrics <- list()
  for(loss in c('Brier','classification_error','log_loss')) {
    h <- function(ii) {
      prob <- model$probability(ii); y <- d$y[ii]
      switch(loss,Brier=(prob-y)^2,classification_error=as.numeric((prob>=.5)!=y),
             log_loss=-y*log(prob)-(1-y)*log1p(-prob))
    }
    truth <- mean(h(truth_ids)); values <- h(evaluation)
    for(method in names(fits)) {
      w <- pmax(0,fits[[method]]$predict(eval_z))
      # Source mean uses its empirical mixture; ratio estimates use fitting pi.
      estimate <- if(method=='Unweighted') mean(values) else
        exp_weighted_mean(values,w,ec,pi)
      metrics[[length(metrics)+1L]] <- data.frame(area='model_audit',dataset=d$name,
        scenario=scenario,seed=seed,method=method,metric=loss,estimate=estimate,
        truth=truth,error=estimate-truth,absolute_error=abs(estimate-truth))
    }
  }
  used <- c(model$training_ids,d$id[c(source,evaluation,target,truth_ids)])
  stopifnot(!anyDuplicated(used))
  list(metrics=do.call(rbind,metrics),diagnostics=exp_diagnostics(fits,eval_z),
       manifest=data.frame(id=used,role=c(rep('policy',length(model$training_ids)),
         rep('ratio',length(source)),rep('evaluation',length(evaluation)),
         rep('target_fit',length(target)),rep('target_truth',length(truth_ids)))),
       source_proportions=pi,
       qualification='Semi-synthetic X-only selective labeling on real observations; design weights are population-law references, finite split/stratum normalization is approximate.')
}

exp_logging <- function(z,client) {
  switch(as.character(client),
    '1'=.05+.25*plogis(z),
    '2'=.15+.7*plogis(-.7*z),
    '3'=.70+.25*plogis(z))
}
exp_target_policy <- function(z) .05+.90*plogis(2*z)
exp_bandit <- function(d,seed,config,scenario='new_policy') {
  p <- exp_partition(d,seed); model <- exp_policy_model(d,p$policy)
  log_data <- function(ii) {
    ii <- app_take(ii,config$sample_size)
    # Random assignment of contexts to historical logging policies. Client
    # heterogeneity is in P(action | context), hence in the joint source law.
    client <- sample(rep(1:3,length.out=length(ii)))
    z <- model$score(ii)
    b <- vapply(seq_along(ii),function(i) exp_logging(z[i],client[i]),numeric(1))
    a <- rbinom(length(ii),1,b)
    list(ids=ii,client=client,x=cbind(z,2*a-1),reward=as.numeric(a==d$y[ii]))
  }
  source <- log_data(p$ratio); evaluation <- log_data(p$evaluation)
  clients <- split(seq_along(source$ids),source$client); sp <- exp_validate(clients)
  pi <- lengths(sp$fit)/sum(lengths(sp$fit))
  target <- app_take(p$target_fit,config$sample_size)
  tz <- model$score(target); ta <- rbinom(length(target),1,exp_target_policy(tz))
  q <- cbind(tz,2*ta-1); qn <- floor(.8*nrow(q))
  fits <- app_ratio_fits(lapply(sp$fit,function(ii) source$x[ii,,drop=FALSE]),
    q[seq_len(qn),,drop=FALSE],lapply(sp$val,function(ii) source$x[ii,,drop=FALSE]),
    q[-seq_len(qn),,drop=FALSE],config,seed+7000L)
  fits$ClassifierSpline <- exp_spline_ratio(lapply(sp$fit,function(ii) source$x[ii,,drop=FALSE]),
    q[seq_len(qn),,drop=FALSE],lapply(sp$val,function(ii) source$x[ii,,drop=FALSE]),q[-seq_len(qn),,drop=FALSE])
  fits$KnownLogging <- list(predict=function(x) {
    target <- exp_target_policy(x[,1]); action <- x[,2]>0
    logging <- Reduce(`+`,lapply(1:3,function(c) pi[c]*exp_logging(x[,1],c)))
    ifelse(action,target/logging,(1-target)/(1-logging))
  },validation=NA_real_)
  truth_ids <- app_take(p$target_truth,2*config$sample_size)
  target_action <- exp_target_policy(model$score(truth_ids))
  truth <- mean(ifelse(d$y[truth_ids]==1,target_action,1-target_action))
  ec <- split(seq_along(evaluation$ids),evaluation$client)
  metrics <- lapply(names(fits),function(method) {
    w <- pmax(0,fits[[method]]$predict(evaluation$x))
    estimate <- exp_weighted_mean(evaluation$reward,w,ec,pi)
    data.frame(area='policy_evaluation',dataset=d$name,scenario=scenario,seed=seed,
      method=method,metric='policy_value',estimate=estimate,truth=truth,
      error=estimate-truth,absolute_error=abs(estimate-truth))
  })
  reward_fit <- exp_reward_model(source,sp)
  pt <- exp_target_policy(tz)
  dm <- mean(pt*reward_fit$predict(cbind(tz,1))+
               (1-pt)*reward_fit$predict(cbind(tz,-1)))
  add_value <- function(method,estimate) data.frame(area='policy_evaluation',
      dataset=d$name,scenario=scenario,seed=seed,method=method,metric='policy_value',
      estimate=estimate,truth=truth,error=estimate-truth,absolute_error=abs(estimate-truth))
  metrics[[length(metrics)+1L]] <- add_value('DirectRewardModel',dm)
  for(method in intersect(c('GD','RF','Kernel','ClassifierSpline'),names(fits))) {
    w <- pmax(0,fits[[method]]$predict(evaluation$x))
    residual <- evaluation$reward-reward_fit$predict(evaluation$x)
    metrics[[length(metrics)+1L]] <- add_value(paste0(method,'_DR'),
        dm+exp_weighted_mean(residual,w,ec,pi))
  }
  used <- c(model$training_ids,d$id[c(source$ids,evaluation$ids,target,truth_ids)])
  stopifnot(!anyDuplicated(used))
  list(metrics=do.call(rbind,metrics),diagnostics=exp_diagnostics(fits,evaluation$x),
       manifest=data.frame(id=used,role=c(rep('policy',length(model$training_ids)),
         rep('ratio',length(source$ids)),rep('evaluation',length(evaluation$ids)),
         rep('target_fit',length(target)),rep('target_truth',length(truth_ids)))),
       source_proportions=pi,
       qualification='Classification-to-bandit reduction: contexts/outcomes real; historical actions simulated. Reward correctness only observed for logged action. Known logging is a diagnostic unavailable in the intended missing-propensity use case.')
}

exp_pu <- function(d,seed,config,scenario='case_only') {
  p <- exp_partition(d,seed)
  source <- app_take(p$ratio,config$sample_size)
  evaluation <- app_take(p$evaluation,2*config$sample_size)
  positive_pool <- p$target_fit[d$y[p$target_fit]==1]
  target <- app_take(positive_pool,config$sample_size)
  if(length(target)<30) stop('Too few independently labeled positive cases')
  prep <- app_preprocess(list(d$x[source,,drop=FALSE]))
  x <- app_transform(d$x[source,,drop=FALSE],prep)
  signal <- x[,1]
  # Rank strata allocate all rows, including ties in sparse predictors, to
  # three nonempty sites. The pooled distribution is left unchanged.
  ordered <- order(signal)
  clients <- split(ordered,pmin(3L,ceiling(seq_along(ordered)/(length(ordered)/3))))
  sp <- exp_validate(clients)
  q <- app_transform(d$x[target,,drop=FALSE],prep); qn <- floor(.8*nrow(q))
  fits <- app_ratio_fits(lapply(sp$fit,function(ii) x[ii,,drop=FALSE]),
    q[seq_len(qn),,drop=FALSE],lapply(sp$val,function(ii) x[ii,,drop=FALSE]),
    q[-seq_len(qn),,drop=FALSE],config,seed+7000L)
  fits$ClassifierTuned <- exp_tuned_classifier(lapply(sp$fit,function(ii) x[ii,,drop=FALSE]),
    q[seq_len(qn),,drop=FALSE],lapply(sp$val,function(ii) x[ii,,drop=FALSE]),q[-seq_len(qn),,drop=FALSE])
  center <- colMeans(q[seq_len(qn),,drop=FALSE])
  fits$PositiveDistance <- list(predict=function(x) -sqrt(rowSums(sweep(x,2,center)^2)),validation=NA_real_)
  # Positive-density reference tunes bandwidth only with positive validation X.
  med <- median_pairwise_distance(q[seq_len(qn),,drop=FALSE]); best <- -Inf
  for(mult in c(.5,1,2)) {
    density <- .density_provider(q[seq_len(qn),,drop=FALSE],med*mult)
    criterion <- mean(log(pmax(1e-300,density(q[-seq_len(qn),,drop=FALSE]))))
    if(criterion>best) { best<-criterion; selected_density<-density }
  }
  fits$PositiveKDE <- list(predict=selected_density,validation=NA_real_)
  fits$Unweighted <- NULL
  xe <- app_transform(d$x[evaluation,,drop=FALSE],prep)
  metrics <- list()
  for(method in names(fits)) {
    score <- fits[[method]]$predict(xe)
    m <- app_rank_metrics(d$y[evaluation],score)
    for(k in names(m)) metrics[[length(metrics)+1L]] <- data.frame(area='positive_unlabeled',
       dataset=d$name,scenario=scenario,seed=seed,method=method,metric=k,
       estimate=m[k],truth=NA_real_,error=NA_real_,absolute_error=NA_real_)
  }
  used <- d$id[c(source,evaluation,target)]; stopifnot(!anyDuplicated(used))
  list(metrics=do.call(rbind,metrics),diagnostics=exp_diagnostics(fits,xe),
       manifest=data.frame(id=used,role=c(rep('ratio',length(source)),
         rep('evaluation',length(evaluation)),rep('positive_reference',length(target)))),
       qualification='Verified positive cases sampled randomly within the positive class (SCAR). Larger p_positive/p_unlabeled ranks positives. Class prior not required for ranking; no calibrated posterior claim.')
}
