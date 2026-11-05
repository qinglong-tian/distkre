# Numerical protocols. Raw-data pooling below occurs only in explicitly
# centralized offline baselines, not in the proposed ratio or predictor fit.
app_take <- function(ii,n) if(length(ii)>n) sample(ii,n) else ii

app_split <- function(d,seed,target=NULL,cap=80L) {
  roles <- app_roles(d,seed); set.seed(seed+1000L)
  manifest <- data.frame(id=d$id,block=d$block,site=d$site,role=roles,
                         used=FALSE,stringsAsFactors=FALSE)
  if(d$type=='anomaly') {
    # Selection uses labels solely to construct a controlled contamination
    # benchmark. Ratio estimation, validation and scoring never receive them.
    make <- function(role) which(roles==role)
    pools <- lapply(c('ratio','model','calibration','test'),make)
    names(pools) <- c('ratio','model','calibration','test')
    contaminate <- function(pool) {
      good <- sample(pool[d$y[pool]==0]); bad <- sample(pool[d$y[pool]==1])
      eps <- c(.05,.15,.25)
      n <- min(cap, floor(length(good)/sum(1-eps)),floor(length(bad)/sum(eps)))
      if(n<12) stop('Too few independent clean/fault observations for this split')
      out <- vector('list',3); gi <- bi <- 0L
      for(c in 1:3) {
        nb <- max(1L,round(n*eps[c])); ng <- n-nb
        out[[c]] <- c(good[gi+seq_len(ng)],bad[bi+seq_len(nb)])
        gi <- gi+ng; bi <- bi+nb
      }
      list(clients=out,remaining=setdiff(good,unlist(out)))
    }
    a <- contaminate(pools$ratio); b <- contaminate(pools$test)
    # Reserve clean ratio-role observations for Q before constructing clients.
    clean <- pools$ratio[d$y[pools$ratio]==0]
    q <- app_take(clean,max(12L,floor(length(clean)*.35)))
    a <- contaminate(setdiff(pools$ratio,q))
    source <- a$clients
    test <- unlist(b$clients)
    model <- list()
    calibration <- list(app_take(pools$calibration[d$y[pools$calibration]==0],3*cap))
    target_label <- 'trusted-reference'
  } else {
    if(!is.null(target)) {
      site <- d$site; target_label <- target
      istarget <- site==target
    } else {
      # Controlled X-only shift, bounded source/target selection propensity.
      # Predeclared physical transformations; no outcome or test-fitted scaling.
      v <- if(d$name=='spambase') log1p(d$x[,55]) else d$x[,1]
      g <- sin(v)+.5*cos(v/2)
      prob <- .15+.7*plogis(1.5*g)
      istarget <- runif(nrow(d$x))<prob
      client_signal <- sin(d$x[,2])
      site <- vapply(client_signal,function(v) {
        p <- exp(c(-v,0,v)); p <- .15+.55*p/sum(p)
        paste0('client',sample.int(3L,1L,prob=p))
      },character(1))
      # Every selected source observation is assigned; none discarded to make sites.
      target_label <- 'controlled-X-shift'
    }
    sites <- sort(unique(site[!istarget]))
    source <- lapply(sites,function(s) app_take(which(!istarget & site==s & roles=='ratio'),cap))
    model <- lapply(sites,function(s) app_take(which(!istarget & site==s & roles=='model'),cap))
    calibration <- lapply(sites,function(s) app_take(which(!istarget & site==s & roles=='calibration'),cap))
    # Exactly match model/calibration site proportions to the ratio estimand.
    n <- 4L*floor(pmin(lengths(source),lengths(model),lengths(calibration))/4L)
    if(any(n<5)) stop('Insufficient source observations in one or more roles')
    source <- Map(function(x,k) x[seq_len(k)],source,n)
    model <- Map(function(x,k) x[seq_len(k)],model,n)
    calibration <- Map(function(x,k) x[seq_len(k)],calibration,n)
    q <- app_take(which(istarget & roles=='ratio'),3*cap)
    test <- app_take(which(istarget & roles=='test'),3*cap)
  }
  if(length(q)<12 || length(test)<12) stop('Insufficient target fit/test observations')
  # Held-out ratio validation; fixed site allocation proportions preserved.
  frac <- .75
  fit <- lapply(source,function(ii) ii[seq_len(max(3,floor(length(ii)*frac)))])
  val <- Map(setdiff,source,fit)
  # Splitting each client can cause rounding changes; criterion explicitly uses
  # fitting proportions rather than validation sample proportions.
  qfit <- q[seq_len(floor(length(q)*frac))]; qval <- setdiff(q,qfit)
  selections <- list(source_fit=fit,source_validation=val,source_model=model,
                     source_calibration=calibration,target_fit=list(qfit),
                     target_validation=list(qval),target_test=list(test))
  manifest$role <- 'unused'
  for(nm in names(selections)) for(ii in selections[[nm]]) {
    stopifnot(!any(manifest$used[ii])); manifest$used[ii] <- TRUE
    manifest$role[ii] <- nm
  }
  # Dependence blocks cannot cross functional roles. Source and target fit are
  # one role family; test/calibration are isolated from all estimation roles.
  list(fit=fit,val=val,qfit=qfit,qval=qval,model=model,cal=calibration,
       test=test,manifest=manifest,target=target_label)
}

app_ratio_fits <- function(xs,z,xv,zv,config,seed) {
  bandwidth <- median_pairwise_distance(z)
  if(!is.finite(bandwidth)||bandwidth<=0) stop('Degenerate target bandwidth')
  pi <- vapply(xs,nrow,integer(1)); pi <- pi/sum(pi)
  validation <- function(fun) sum(pi*vapply(xv,function(x) mean(fun(x)^2)/2,numeric(1)))-mean(fun(zv))
  result <- list(); grid <- expand.grid(lambda=config$lambdas,mult=config$bandwidths)
  methods <- if(is.null(config$methods)) c('GD','RF','Kernel') else config$methods
  for(method in methods) {
    best <- Inf; cost <- 0
    for(i in seq_len(nrow(grid))) {
      bw <- bandwidth*grid$mult[i]; lambda <- grid$lambda[i]
      fit <- if(method=='Kernel') fit_kulsif(do.call(rbind,xs),z,lambda,bw) else
        dkre_rff(xs,z,lambda,config$features,bw,
                 method=if(method=='GD') 'gd' else 'closed',
                 iterations=config$iterations,seed=seed)
      fun <- local({f <- fit; function(x) as.numeric(predict(f,x))})
      loss <- validation(fun)
      if(method=='GD') cost <- cost+fit$metadata$communication_coordinates+
        length(xs)*(config$features+1) # candidate model/scalar validation
      if(loss<best) { best <- loss; chosen <- list(predict=fun,fit=fit,
           validation=loss,lambda=lambda,bandwidth=bw) }
    }
    chosen$tuning_coordinates <- cost
    result[[method]] <- chosen
  }
  # Centralized domain-classifier reference, class 1 = target.
  pooled <- do.call(rbind,xs); a <- rbind(pooled,z)
  yy <- c(rep(0,nrow(pooled)),rep(1,nrow(z)))
  fit <- app_predictor(list(a),list(yy),list(rep(1,nrow(a))),'binary',rho=.01)
  prior <- nrow(pooled)/nrow(z)
  fun <- function(x) { p <- pmin(1-1e-8,pmax(1e-8,fit$predict(x))); prior*p/(1-p) }
  result$Classifier <- list(predict=fun,validation=validation(fun),
                             tuning_coordinates=NA_real_)
  result$Unweighted <- list(predict=function(x) rep(1,nrow(x)),
                            validation=validation(function(x) rep(1,nrow(x))),
                            tuning_coordinates=0)
  result
}

app_predictor <- function(xs,ys,weights,type,rho=.01) {
  xs <- lapply(xs,function(x) cbind(1,x)); p <- ncol(xs[[1]])
  mass <- sum(vapply(weights,sum,numeric(1))); if(mass<=0) stop('Zero weight mass')
  if(type=='regression') {
    mu <- sum(vapply(ys,sum,numeric(1)))/sum(lengths(ys))
    sy <- sqrt(sum(vapply(ys,function(y) sum((y-mu)^2),numeric(1)))/sum(lengths(ys)))
    if(sy<1e-8) sy <- 1
    ys <- lapply(ys,function(y) (y-mu)/sy)
  } else { mu <- 0; sy <- 1 }
  calls <- 0L; gradcalls <- 0L
  penalty <- c(0,rep(rho,p-1))
  objective <- function(b) {
    calls <<- calls+1L
    sum(vapply(seq_along(xs),function(c) {
      eta <- as.numeric(xs[[c]]%*%b); y <- ys[[c]]
      loss <- if(type=='binary') pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta else (eta-y)^2/2
      sum(weights[[c]]*loss)
    },numeric(1)))/mass+sum(penalty*b^2)/2
  }
  gradient <- function(b) {
    gradcalls <<- gradcalls+1L
    Reduce(`+`,lapply(seq_along(xs),function(c) {
      eta <- as.numeric(xs[[c]]%*%b)
      residual <- if(type=='binary') plogis(eta)-ys[[c]] else eta-ys[[c]]
      as.numeric(crossprod(xs[[c]],weights[[c]]*residual))
    }))/mass+penalty*b
  }
  fit <- optim(numeric(p),objective,gradient,method='BFGS',control=list(maxit=500,reltol=1e-9))
  list(predict=function(x) {
    eta <- as.numeric(cbind(1,x)%*%fit$par)
    if(type=='binary') plogis(eta) else mu+sy*eta
  },convergence=fit$convergence,
  coordinates=length(xs)*(calls*(p+1)+gradcalls*2*p))
}

app_weighted_radius <- function(scores,weights,new_weight,alpha,tolerance=1e-7) {
  total <- sum(vapply(weights,sum,numeric(1)))
  needed <- (1-alpha)*(total+new_weight)
  if(total<=0 || needed>total) return(list(radius=Inf,queries=0L))
  lo <- 0; hi <- max(vapply(scores,max,numeric(1))); queries <- 0L
  # Aggregate CDF queries only; upper endpoint is a conservative numerical quantile.
  while(hi-lo>tolerance && queries<80L) {
    mid <- (lo+hi)/2
    below <- sum(vapply(seq_along(scores),function(c) sum(weights[[c]][scores[[c]]<=mid]),numeric(1)))
    if(below>=needed) hi <- mid else lo <- mid
    queries <- queries+1L
  }
  list(radius=hi,queries=queries)
}

app_rank_metrics <- function(y,score) {
  if(length(unique(y))!=2) stop('Both classes required for ranking metrics')
  # Group ties so AUPRC is invariant to row ordering.
  thresholds <- sort(unique(score),decreasing=TRUE)
  tp <- fp <- numeric(length(thresholds))
  for(i in seq_along(thresholds)) { pos <- score>=thresholds[i]; tp[i] <- sum(y[pos]==1); fp[i] <- sum(y[pos]==0) }
  recall <- tp/sum(y==1); precision <- tp/(tp+fp)
  auc <- (sum(rank(score)[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0))
  c(AUPRC=sum(diff(c(0,recall))*precision),AUROC=auc)
}
