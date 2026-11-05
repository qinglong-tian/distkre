# Deterministic simulation engine. Source the maintained package's R/ files first.
sim_seed <- function(config,case,replication,phase,key) {
  as.integer(config$seed + if(phase=='pilot') 10000000L else 20000000L) +
    as.integer(case$index*100000L+replication*100L+key)
}
sim_means <- function(d=5) {
  cbind(rbind(c(-1,-1),c(-1,1),c(1,-1),c(1,1)),matrix(0,4,d-2))
}
sim_allocate <- function(n,weights) {
  z<-floor(n*weights);rem<-n-sum(z)
  if(rem>0) {ix<-order(n*weights-z,decreasing=TRUE)[seq_len(rem)];z[ix]<-z[ix]+1L}
  as.integer(z)
}
sim_draw <- function(n,weights,means) {
  k<-sample.int(nrow(means),n,replace=TRUE,prob=weights)
  matrix(rnorm(n*ncol(means)),n,ncol(means))+means[k,,drop=FALSE]
}
sim_client_law <- function(case,c) {
  w<-(1-case$tau)*case$pi;w[c]<-w[c]+case$tau;w
}
sim_logmix <- function(x,weights,means) {
  L<-vapply(seq_len(nrow(means)),function(k) {
    rowSums(matrix(dnorm(sweep(x,2,means[k,],'-'),log=TRUE),nrow(x)))+log(weights[k])
  },numeric(nrow(x)))
  mx<-apply(L,1,max);mx+log(rowSums(exp(L-mx)))
}
sim_truth <- function(x,case,config) {
  mu<-sim_means(config$dimension)
  exp(sim_logmix(x,config$rho,mu)-sim_logmix(x,case$pi,mu))
}
sim_data <- function(case,replication,phase,config) {
  mu<-sim_means(config$dimension)
  total<-sim_allocate(case$N,case$pi);nf<-sim_allocate(round(.8*case$N),case$pi)
  nv<-total-nf; mf<-round(.8*case$M);mv<-case$M-mf
  set.seed(sim_seed(config,case,replication,phase,1))
  allx<-lapply(1:4,function(c)sim_draw(total[c],sim_client_law(case,c),mu))
  xs<-lapply(1:4,function(c)allx[[c]][seq_len(nf[c]),,drop=FALSE])
  xv<-lapply(1:4,function(c)allx[[c]][nf[c]+seq_len(nv[c]),,drop=FALSE])
  set.seed(sim_seed(config,case,replication,phase,2))
  allz<-sim_draw(case$M,config$rho,mu)
  set.seed(sim_seed(config,case,replication,phase,4))
  ns<-sim_allocate(config$test_points,case$pi)
  testlist<-lapply(1:4,function(c)sim_draw(ns[c],sim_client_law(case,c),mu))
  test<-do.call(rbind,testlist)
  wt<-rep(case$pi/ns,ns)
  truth<-sim_truth(test,case,config)
  stopifnot(all(is.finite(truth)),max(truth)<=max(config$rho/case$pi)+1e-10)
  set.seed(sim_seed(config,case,replication,phase,3))
  masterR<-max(config$features,config$mechanism_features)
  list(xs=xs,xv=xv,z=allz[seq_len(mf),,drop=FALSE],
    zv=allz[mf+seq_len(mv),,drop=FALSE],nf=nf,nv=nv,pi=nf/sum(nf),
    test=test,tw=wt,truth=truth,test_client=rep(1:4,ns),
    master=list(frequencies=matrix(rnorm(masterR*config$dimension),masterR),
                phases=runif(masterR,0,2*pi)),
    sigma=median_pairwise_distance(allz[seq_len(mf),,drop=FALSE]))
}
sim_spec <- function(data,R,sigma) rff_features(data$z,R,sigma,
  frequencies=data$master$frequencies[seq_len(R),,drop=FALSE]/sigma,
  phases=data$master$phases[seq_len(R)])
sim_metrics <- function(pred,truth,weights) {
  if(any(!is.finite(pred))) stop('Nonfinite test prediction')
  mse<-sum(weights*(pred-truth)^2)
  list(mse=mse,rmse=sqrt(mse),calibration=sum(weights*pred),negative=sum(weights*(pred<0)))
}
sim_batches <- function(n,size) split(seq_len(n),ceiling(seq_len(n)/size))
sim_rf_predict <- function(points,spec,coefficients,batch=500) {
  if(is.null(dim(coefficients))) coefficients<-matrix(coefficients,ncol=1)
  out<-matrix(0,nrow(points),ncol(coefficients))
  for(ix in sim_batches(nrow(points),batch))out[ix,]<-evaluate_rff(points[ix,,drop=FALSE],spec)%*%coefficients
  out
}
sim_gd <- function(target_mean,weights,exchange,lambda,T,checkpoints=T) {
  R<-length(target_mean);b<-numeric(R);eta<-1/(2+lambda)
  cp<-sort(unique(as.integer(checkpoints)));cp<-cp[cp>=0&cp<=T]
  path<-matrix(NA_real_,R,length(cp),dimnames=list(NULL,as.character(cp)))
  elapsed<-setNames(numeric(length(cp)),as.character(cp));start<-proc.time()[['elapsed']]
  if(0L%in%cp)path[,'0']<-b
  nextcp<-match(TRUE,cp>0,nomatch=length(cp)+1L)
  for(t in seq_len(T)) {
    products<-exchange(b)
    # Only vectors cross the client/server interface.
    if(!is.list(products)||length(products)!=length(weights)||
       !all(vapply(products,function(p)is.numeric(p)&&is.null(dim(p))&&length(p)==R,logical(1))))stop('Invalid client message')
    p<-Reduce('+',Map(function(x,w)x*w,products,weights))
    b<-b-eta*(p-target_mean+lambda*b)
    if(nextcp<=length(cp)&&t==cp[nextcp]) {
      path[,nextcp]<-b;elapsed[nextcp]<-proc.time()[['elapsed']]-start;nextcp<-nextcp+1L
    }
  }
  if(any(!is.finite(b)))stop('Nonfinite GD coefficients')
  list(coefficients=b,path=path,checkpoint_seconds=elapsed,iterations=T,exchanges=T,
       coordinates=2*length(weights)*R*T)
}
sim_rf_cache <- function(data,R,sigma) {
  spec<-sim_spec(data,R,sigma)
  clients<-lapply(data$xs,dkre_client,rff_spec=spec)
  val_moments<-lapply(data$xv,function(x)crossprod(evaluate_rff(x,spec))/nrow(x))
  target_mean<-colMeans(evaluate_rff(data$z,spec))
  val_target_mean<-colMeans(evaluate_rff(data$zv,spec))
  exchange<-function(b)lapply(clients,dkre_client_update,coefficients=b)
  validate<-function(b) {
    source<-vapply(val_moments,function(A)as.numeric(crossprod(b,A%*%b)),numeric(1))
    sum(data$pi*source)/2-sum(b*val_target_mean)
  }
  list(spec=spec,clients=clients,exchange=exchange,validate=validate,a=target_mean,
       H0=Reduce('+',Map(function(cl,w)cl$second_moment*w,clients,data$pi)))
}
sim_select <- function(table) {
  candidates<-which(is.finite(table$validation))
  if(!length(candidates))stop('All tuning candidates failed')
  best<-min(table$validation[candidates]);tol<-1e-10
  ix<-candidates[table$validation[candidates]<=best+tol]
  ix[order(-table$lambda[ix],-table$sigma[ix])][1]
}
sim_kernel_cache <- function(xs,z,xv,zv,weights,sigma) {
  X<-do.call(rbind,xs);V<-do.call(rbind,xv);points<-rbind(V,zv)
  n<-nrow(X);m<-nrow(z)
  list(X=X,z=z,n=n,m=m,sigma=sigma,
    K=.gaussian_kernel(X,X,sigma),q=rowMeans(.gaussian_kernel(X,z,sigma)),
    VX=.gaussian_kernel(points,X,sigma),qV=rowMeans(.gaussian_kernel(points,z,sigma)),
    nv=nrow(V),valweights=rep(weights/vapply(xv,nrow,integer(1)),vapply(xv,nrow,integer(1))))
}
sim_kernel_coeff <- function(cache,lambda) as.numeric(.solve_linear(
  cache$K/cache$n+diag(lambda,cache$n),-cache$q/(cache$n*lambda)))
sim_kernel_validate <- function(cache,alpha,lambda) {
  pred<-as.numeric(cache$VX%*%alpha+cache$qV/lambda)
  sum(cache$valweights*pred[seq_len(cache$nv)]^2)/2-mean(pred[-seq_len(cache$nv)])
}
sim_kernel_predict <- function(points,cache,alpha,lambda,batch=500) {
  pred<-numeric(nrow(points))
  for(ix in sim_batches(nrow(points),batch))pred[ix]<-as.numeric(
    .gaussian_kernel(points[ix,,drop=FALSE],cache$X,cache$sigma)%*%alpha+
      rowMeans(.gaussian_kernel(points[ix,,drop=FALSE],cache$z,cache$sigma))/lambda)
  pred
}
sim_kernel_fit <- function(xs,z,xv,zv,weights,base_sigma,config,tuned=TRUE) {
  sigmas<-if(tuned)base_sigma*config$bandwidth_multipliers else base_sigma
  lambdas<-if(tuned)config$lambdas else .01
  candidates<-list();best<-NULL
  for(sigma in sigmas) {
    cache<-sim_kernel_cache(xs,z,xv,zv,weights,sigma)
    for(lambda in lambdas) {
      alpha<-sim_kernel_coeff(cache,lambda)
      value<-sim_kernel_validate(cache,alpha,lambda)
      row<-data.frame(sigma=sigma,lambda=lambda,validation=value)
      candidates[[length(candidates)+1]]<-row
      if(is.null(best)||value<best$value-1e-10||
         (abs(value-best$value)<=1e-10 && (lambda>best$lambda||lambda==best$lambda&&sigma>best$sigma))) {
        best<-list(value=value,lambda=lambda,sigma=sigma,alpha=alpha)
      }
    }
    rm(cache);gc(FALSE)
  }
  # Rebuild only the selected representation. Cache construction is included in timing.
  cache<-sim_kernel_cache(xs,z,xv,zv,weights,best$sigma)
  alpha<-sim_kernel_coeff(cache,best$lambda)
  list(cache=cache,alpha=alpha,lambda=best$lambda,sigma=best$sigma,
       validation=best$value,candidates=do.call(rbind,candidates))
}
sim_kde_matrix <- function(points,x,bandwidths,batch) {
  out<-matrix(0,nrow(points),length(bandwidths))
  for(ix in sim_batches(nrow(points),batch)) {
    # Distances are shared across bandwidth candidates; no denominator floor is added.
    d2<-outer(rowSums(points[ix,,drop=FALSE]^2),rowSums(x^2),'+')-2*tcrossprod(points[ix,,drop=FALSE],x)
    d2<-pmax(d2,0)
    for(k in seq_along(bandwidths)) {
      h<-bandwidths[k]
      out[ix,k]<-rowMeans(exp(-d2/(2*h*h)))*(2*pi*h*h)^(-ncol(x)/2)
    }
  }
  out
}
sim_kde_fit_predict <- function(data,config,tuned=TRUE) {
  X<-do.call(rbind,data$xs);V<-do.call(rbind,data$xv);points<-rbind(V,data$zv)
  mult<-if(tuned)config$bandwidth_multipliers else 1
  hx<-.kde_bandwidth(X)*mult;hz<-.kde_bandwidth(data$z)*mult
  src<-sim_kde_matrix(points,X,hx,config$batch_size)
  tgt<-sim_kde_matrix(points,data$z,hz,config$batch_size)
  vw<-rep(data$pi/data$nv,data$nv);nv<-nrow(V)
  grid<-expand.grid(i=seq_along(hx),j=seq_along(hz))
  grid$validation<-vapply(seq_len(nrow(grid)),function(k) {
    r<-tgt[,grid$j[k]]/src[,grid$i[k]]
    sum(vw*r[seq_len(nv)]^2)/2-mean(r[-seq_len(nv)])
  },numeric(1))
  good<-which(is.finite(grid$validation));if(!length(good))stop('All KDE candidates failed')
  k<-good[which.min(grid$validation[good])]
  src_test<-sim_kde_matrix(data$test,X,hx[grid$i[k]],config$batch_size)[,1]
  tgt_test<-sim_kde_matrix(data$test,data$z,hz[grid$j[k]],config$batch_size)[,1]
  list(pred=tgt_test/src_test,source_bandwidth=hx[grid$i[k]],target_bandwidth=hz[grid$j[k]],
       validation=grid$validation[k],candidates=grid)
}
