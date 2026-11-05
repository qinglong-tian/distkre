#!/usr/bin/env Rscript
font_cache<-file.path(tempdir(),'simulation-font-cache')
dir.create(font_cache,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(XDG_CACHE_HOME=font_cache)
args<-commandArgs(trailingOnly=TRUE)
getarg<-function(name,default){x<-args[startsWith(args,paste0('--',name,'='))];if(length(x))sub(paste0('^--',name,'='),'',tail(x,1))else default}
run_id<-getarg('run-id','three-studies-20261102-v1')
run_dir<-getarg('run-directory',file.path('results','simulation-v2',run_id))
phase<-getarg('phase','production')
make_figures<-getarg('plots','true')=='true'
source('reproduction/engine/simulations/R/studies.R')
manifest<-readRDS(file.path(run_dir,'manifest.rds'))
kernel_cap<-manifest$config$exact_kernel_cap
kernel_coverage_note<-if(!is.null(manifest$kernel_extension))
 'Full-kernel reference at all source sample sizes' else paste('N <=',kernel_cap,'kernel cap')
files<-list.files(file.path(run_dir,phase),pattern='^rep-.*\\.rds$',recursive=TRUE,full.names=TRUE)
if(!length(files))stop('No completed jobs for ',phase)
set.seed(1729L) # Independent analysis-only bootstrap stream.
bootstrap_ci<-function(v,B=2000L) {
 if(length(v)<2L)return(c(NA_real_,NA_real_))
 as.numeric(quantile(replicate(B,mean(sample(v,length(v),replace=TRUE))),c(.025,.975),names=FALSE))
}
rows<-list();jobs<-list();fatal<-list()
for(f in files) {
 x<-readRDS(f)
 if(is.null(x$rows)) {fatal[[length(fatal)+1]]<-data.frame(path=f,error=x$error);next}
 rows[[length(rows)+1]]<-x$rows
 jobs[[length(jobs)+1]]<-data.frame(case_id=x$metadata$case$id,
   elapsed=x$metadata$elapsed,status=x$status,replication=x$rows$replication[1])
}
raw<-sim_bind(rows);jobtable<-do.call(rbind,jobs)
out<-file.path(run_dir,'summary',phase);dir.create(out,recursive=TRUE,showWarnings=FALSE)
# A Gaussian KDE ratio can have infinite L2(P) risk despite bounded truth.
# For common unit-covariance source components, a negative radial tail margin
# proves divergence. Preserve the fitted method and expose this diagnostic.
raw$timing_scope<-ifelse(raw$study=='communication','checkpoint_or_reference_compute',
  ifelse(raw$method%in%c('distributed_gd','rff_closed','rff_at_gd_parameters','gd_uniform_weight_ablation'),
         'fit_and_tuning_excludes_test_prediction','includes_test_prediction'))
# Initialization and optional final delivery are kept separate from 2CRT.
rfrow<-is.finite(raw$R)&raw$method%in%c('distributed_gd','rff_closed','rff_at_gd_parameters','gd_uniform_weight_ablation')
tunedrow<-raw$study%in%c('accuracy','imbalance')&raw$method%in%c('distributed_gd','rff_closed')
raw$sample_count_coordinates<-ifelse(rfrow,4,NA_real_)
raw$feature_initialization_coordinates<-ifelse(rfrow,4*5*ifelse(tunedrow,length(manifest$config$bandwidth_multipliers),1),NA_real_)
raw$optional_final_model_coordinates<-ifelse(rfrow,4*raw$R,NA_real_)
raw$validation_coordinates<-ifelse(tunedrow,length(manifest$config$lambdas)*length(manifest$config$bandwidth_multipliers)*(4*raw$R+4),0)
raw$tuning_optimization_coordinates<-raw$tuning_coordinates-raw$validation_coordinates
raw$dense_matrix_reference_coordinates<-ifelse(raw$method=='rff_closed',4*raw$R^2,NA_real_)
raw$kde_tail_margin<-NA_real_
if('target_bandwidth'%in%names(raw)) {
 i<-raw$method=='kde_pooled'&raw$status=='ok'
 raw$kde_tail_margin[i]<-.5-1/raw$sigma[i]^2+1/raw$target_bandwidth[i]^2
}
write.csv(raw,file.path(out,'replicate-metrics.csv'),row.names=FALSE)
write.csv(jobtable,file.path(out,'job-timings.csv'),row.names=FALSE)
if(length(fatal))write.csv(do.call(rbind,fatal),file.path(out,'fatal-jobs.csv'),row.names=FALSE)
keys<-c('case_id','study','N','M','tau','method','row_type','R')
raw$group_lambda<-ifelse(raw$study=='communication',raw$lambda,NA_real_)
raw$group_iterations<-ifelse(raw$study=='communication',raw$iterations,NA_real_)
raw$group_budget<-ifelse(raw$study=='communication',raw$budget,NA_real_)
keys<-c(keys,'group_lambda','group_iterations','group_budget')
key<-apply(raw[,keys],1,function(x)paste(x,collapse='|'))
groups<-split(seq_len(nrow(raw)),key)
summaries<-lapply(groups,function(ix) {
 z<-raw[ix,];ok<-z$status=='ok'&is.finite(z$rmse);v<-z$rmse[ok];n<-length(v)
 base<-z[1,keys,drop=FALSE]
 base$attempted<-nrow(z);base$successful<-n;base$failure_rate<-mean(!ok)
 base$rmse_mean<-if(n)mean(v)else NA_real_
 base$rmse_sd<-if(n>1)sd(v)else NA_real_
 base$rmse_mcse<-if(n>1)sd(v)/sqrt(n)else NA_real_
 ci<-bootstrap_ci(v);base$rmse_bootstrap_lower<-ci[1];base$rmse_bootstrap_upper<-ci[2]
 base$mse_mean<-if(n)mean(z$mse[ok])else NA_real_
 for(col in c('calibration','negative','seconds','optimization_gap','prediction_gap','fit_coordinates','tuning_coordinates','total_coordinates')) {
  vv<-z[[col]][ok];vv<-vv[is.finite(vv)]
  base[[paste0(col,'_mean')]]<-if(length(vv))mean(vv)else NA_real_
 }
 tailvalues<-z$kde_tail_margin[is.finite(z$kde_tail_margin)]
 base$kde_infinite_L2_fraction<-if(length(tailvalues))mean(tailvalues<0)else NA_real_
 base
})
summary<-do.call(rbind,summaries);rownames(summary)<-NULL
write.csv(summary,file.path(out,'method-summary.csv'),row.names=FALSE)
paired<-list()
accuracy<-raw[raw$study!='communication'&raw$status=='ok'&raw$row_type=='estimate',]
for(id in unique(accuracy$case_id)) {
 z<-accuracy[accuracy$case_id==id,];gd<-z[z$method=='distributed_gd',]
 for(method in setdiff(unique(z$method),'distributed_gd')) {
  ref<-z[z$method==method,];p<-merge(gd[,c('replication','rmse')],ref[,c('replication','rmse')],by='replication',suffixes=c('_gd','_ref'))
  if(!nrow(p))next
  delta<-p$rmse_gd-p$rmse_ref
  ci<-bootstrap_ci(delta)
  paired[[length(paired)+1]]<-data.frame(bootstrap_lower=ci[1],bootstrap_upper=ci[2],case_id=id,reference=method,n_pairs=nrow(p),
    mean_rmse_difference=mean(delta),sd_difference=if(length(delta)>1)sd(delta)else NA_real_,
    mcse_difference=if(length(delta)>1)sd(delta)/sqrt(length(delta))else NA_real_)
 }
}
if(length(paired))write.csv(do.call(rbind,paired),file.path(out,'paired-differences.csv'),row.names=FALSE)
selection<-raw[raw$study!='communication'&raw$status=='ok'&is.finite(raw$lambda),
              c('case_id','method','lambda','sigma','replication','boundary')]
write.csv(selection,file.path(out,'selected-parameters.csv'),row.names=FALSE)
write.csv(aggregate(elapsed~case_id,jobtable,function(x)c(mean=mean(x),median=median(x),total=sum(x))),
          file.path(out,'timing-summary.csv'),row.names=FALSE)
colors<-c(distributed_gd='#D55E00',rff_closed='#0072B2',kernel_pooled='#009E73',
          kde_pooled='#CC79A7',kernel_local_naive='#555555',constant_one='#999999',
          gd_uniform_weight_ablation='#A6761D')
labels<-c(distributed_gd='Distributed GD',rff_closed='RF closed form',kernel_pooled='Kernel KuLSIF',
          kde_pooled='Pooled KDE',kernel_local_naive='Naive local KuLSIF',constant_one='Constant 1',
          gd_uniform_weight_ablation='Uniform-weight ablation')
figure<-function(name,width,height,draw) {
 for(ext in c('svg','png')) {
  path<-file.path(out,paste0(name,'.',ext))
  if(ext=='svg')svg(path,width=width,height=height) else png(path,width=width,height=height,units='in',res=180,type='cairo')
  tryCatch(draw(),finally=dev.off())
 }
}
errorbars<-function(x,y,e,col) {
 good<-is.finite(x)&is.finite(y)&is.finite(e)&e>0
 if(any(good)) {
  xx<-x[good];lo<-pmax(y[good]-e[good],1e-8);hi<-y[good]+e[good]
  left<-grconvertX(grconvertX(xx,from='user',to='inches')-.025,from='inches',to='user')
  right<-grconvertX(grconvertX(xx,from='user',to='inches')+.025,from='inches',to='user')
  segments(xx,lo,xx,hi,col=col)
  segments(left,lo,right,lo,col=col);segments(left,hi,right,hi,col=col)
 }
}
if(make_figures && any(summary$study=='accuracy'))figure('study1-sample-size-accuracy',11,4.8,function() {
 layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(1,.16))
 par(mar=c(4.2,4.5,2.8,1),oma=c(0,0,2,0))
 for(sweep in c('source','target')) {
  s<-summary[summary$study=='accuracy'&summary$row_type=='estimate'&
    if(sweep=='source')summary$M==1000 else summary$N==2000,]
  methods<-intersect(c('distributed_gd','rff_closed','kernel_pooled','constant_one'),unique(s$method))
  s<-s[s$method%in%methods,];x<-if(sweep=='source')s$N else s$M
  yr<-range(c(s$rmse_mean,s$rmse_mean+s$rmse_mcse),finite=TRUE);yr[1]<-max(yr[1]*.8,1e-8);yr[2]<-yr[2]*1.2
  plot(range(x),yr,type='n',log='y',xlab=if(sweep=='source')'Available source N' else 'Available target M',
       ylab='Independent-test RMSE (log scale)',main=if(sweep=='source')'Source sample size' else 'Target sample size')
  for(m in methods) {
   z<-s[s$method==m,];xx<-if(sweep=='source')z$N else z$M;ii<-order(xx)
   lines(xx[ii],z$rmse_mean[ii],col=colors[m],lwd=2,type='b',pch=16,cex=.65)
   errorbars(xx,z$rmse_mean,z$rmse_mcse,colors[m])
  }
 }
 par(mar=rep(0,4));plot.new()
 legend('center',legend=labels[methods],col=colors[methods],lty=1,pch=16,cex=.8,bty='n',horiz=TRUE)
 mtext(paste(toupper(phase),' | +/- one MCSE;',kernel_coverage_note,'; raw KDE diagnostics are reported separately'),outer=TRUE,cex=.75)
})
if(make_figures && any(summary$study=='communication'))figure('study2-optimization-communication',14,4.6,function() {
 par(mfrow=c(1,3),mar=c(4.2,4.6,2.8,1),oma=c(0,0,2,0))
 rs<-c(50,200,400);pal<-c('#009E73','#0072B2','#D55E00')
 s<-summary[summary$study=='communication'&summary$group_lambda==.01,]
 cp<-s[s$row_type=='checkpoint',]
 plot(c(1,20001),c(1e-16,max(cp$optimization_gap_mean,na.rm=TRUE)*1.2),type='n',log='xy',
      xlab='Gradient updates T + 1',ylab='Training objective gap',main='Optimization convergence')
 for(i in seq_along(rs)) {z<-cp[cp$R==rs[i],];z<-z[order(z$group_iterations),]
  lines(z$group_iterations+1,pmax(z$optimization_gap_mean,1e-16),col=pal[i],lwd=2)}
 legend('topright',legend=paste('R =',rs),col=pal,lty=1,lwd=2,bty='n',cex=.75)
 cp<-cp[cp$group_iterations>0,];ref<-s[s$row_type=='reference',]
 plot(range(c(cp$fit_coordinates_mean,ref$fit_coordinates_mean),finite=TRUE),
      range(c(cp$rmse_mean,ref$rmse_mean),finite=TRUE)*c(.85,1.15),type='n',log='xy',
      xlab='Fitting coordinates transmitted',ylab='Independent-test RMSE',main='Accuracy versus communication')
 for(i in seq_along(rs)) {z<-cp[cp$R==rs[i],];z<-z[order(z$fit_coordinates_mean),]
  lines(z$fit_coordinates_mean,z$rmse_mean,col=pal[i],lwd=2)
  zz<-ref[ref$R==rs[i],];points(zz$fit_coordinates_mean,zz$rmse_mean,pch=17,col=pal[i],cex=1.3)}
 legend('topright',legend=c('GD trajectory','RF matrix aggregation'),pch=c(NA,17),lty=c(1,NA),bty='n',cex=.7)
 bu<-s[s$row_type=='budget',]
 plot(range(bu$group_budget,finite=TRUE),range(bu$rmse_mean,finite=TRUE)*c(.9,1.1),type='n',log='xy',
      xlab='Fitting coordinate budget',ylab='Independent-test RMSE',main='Spend on features or iterations?')
 for(i in seq_along(rs)){z<-bu[bu$R==rs[i],];z<-z[order(z$group_budget),];lines(z$group_budget,z$rmse_mean,col=pal[i],type='b',pch=16,lwd=2)}
 mtext(paste(toupper(phase),' | lambda = 0.01; all R/lambda conditions are in the tables; gap display floor = 1e-16'),outer=TRUE,cex=.8)
})
if(make_figures && any(summary$study=='heterogeneity'))figure('study3-heterogeneity-imbalance',12,5.5,function() {
 layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(1,.18))
 par(mar=c(4.5,4.5,2.8,1),oma=c(0,0,2,0))
 methods<-c('distributed_gd','rff_closed','kernel_pooled','kernel_local_naive')
 s<-summary[summary$study=='heterogeneity'&summary$method%in%methods&summary$row_type=='estimate',]
 plot(c(0,1),range(s$rmse_mean,finite=TRUE)*c(.8,1.2),type='n',log='y',
      xlab='Client heterogeneity tau',ylab='Independent-test RMSE',main='Same pooled P and target Q')
 for(m in methods){z<-s[s$method==m,];z<-z[order(z$tau),];lines(z$tau,z$rmse_mean,col=colors[m],lwd=2,type='b',pch=16);errorbars(z$tau,z$rmse_mean,z$rmse_mcse,colors[m])}
 s<-summary[summary$case_id%in%c('accuracy_N2000_M1000','imbalance')&
              summary$method%in%c(methods,'gd_uniform_weight_ablation')&summary$row_type=='estimate',]
 methods<-c(methods,'gd_uniform_weight_ablation')
 plot(c(.8,2.2),range(s$rmse_mean,finite=TRUE)*c(.8,1.2),type='n',log='y',xaxt='n',
      xlab='Source allocation',ylab='Independent-test RMSE',main='Equal versus unequal client sizes')
 axis(1,at=1:2,labels=c('Balanced','70/10/10/10'))
 for(m in methods){z<-s[s$method==m,];xx<-ifelse(z$case_id=='imbalance',2,1)
  lines(xx,z$rmse_mean,col=colors[m],type='b',pch=16,lwd=2);errorbars(xx,z$rmse_mean,z$rmse_mcse,colors[m])}
 par(mar=rep(0,4));plot.new()
 legend('center',legend=labels[methods],col=colors[methods],lty=1,pch=16,bty='n',cex=.75,horiz=TRUE)
 mtext(paste(toupper(phase),' | KDE and all method results are included in the summary tables'),outer=TRUE,cex=.8)
})
if(make_figures && any(summary$method=='kde_pooled'&summary$study=='accuracy'))figure('study1-raw-kde-diagnostic',10,4.8,function() {
 par(mfrow=c(1,2),mar=c(4.2,4.7,2.8,1),oma=c(0,0,2,0))
 for(sweep in c('source','target')) {
  z<-summary[summary$study=='accuracy'&summary$method=='kde_pooled'&
    if(sweep=='source')summary$M==1000 else summary$N==2000,]
  xx<-if(sweep=='source')z$N else z$M;ii<-order(xx)
  plot(xx[ii],z$rmse_mean[ii],type='b',log='y',pch=16,col=colors['kde_pooled'],lwd=2,
    xlab=if(sweep=='source')'Available source N' else 'Available target M',
    ylab='Raw KDE finite-test RMSE',main=if(sweep=='source')'KDE: source sweep' else 'KDE: target sweep')
  text(xx,z$rmse_mean,labels=paste0('divergent: ',round(100*z$kde_infinite_L2_fraction),'%'),
       pos=ifelse(xx==min(xx),4,ifelse(xx==max(xx),2,3)),cex=.65)
 }
 mtext(paste(toupper(phase),' | Divergent = selected fits with provably infinite population L2 risk; finite-test errors are descriptive'),outer=TRUE,cex=.7)
})
manifest<-readRDS(file.path(run_dir,'manifest.rds'))
expected<-if(phase=='pilot')length(manifest$cases)*manifest$config$pilot_replications else sum(vapply(manifest$cases,function(x)x$replications,integer(1)))
tailrows<-raw[is.finite(raw$kde_tail_margin),]
text<-c(paste('# Simulation results:',run_id),paste('Phase:',phase),
 paste('Completed dataset jobs:',length(files),'of',expected),
 paste('Numerical method failures:',sum(raw$status!='ok')),paste('Fatal job files:',length(fatal)),
 '', 'Tables: method-summary.csv, paired-differences.csv, selected-parameters.csv, timing-summary.csv.',
 'Figures are exported as SVG and PNG. Production summaries never include pilot data.',
 '', 'The primary metric is mean independent-test RMSE; Monte Carlo SE uses independent dataset replications.',
 'Tables also include percentile 95% bootstrap intervals from 2000 dataset-level resamples, with a fixed independent analysis seed.',
 'Reference computations and validation communication are recorded separately from GD updates.',
 'Initialization counts (source counts and five numeric feature-specification fields per client/spec) and optional final model delivery are separate columns; core total_coordinates includes fitting and tuning only.',
 'Local KuLSIF estimates client-specific ratios and is explicitly a naive baseline.',
 if(!is.null(manifest$kernel_extension))
  'The original N=2000 cap was extended to N=8000; 400 full-kernel fits use the same original production datasets.' else
  paste0('The kernel baseline is omitted above the prespecified N=',kernel_cap,' cap.'),
 if(!is.null(manifest$kernel_extension))
  'For extended jobs, elapsed times sum the original job and the later kernel fit. The kernel extension used process-local Accelerate BLAS; other methods retain their original timings. These are not directly comparable wall-clock benchmarks.' else
  'Job elapsed times include the complete dataset job; method timing scopes are explicitly labeled and differ in test-prediction inclusion.',
 '', 'KDE tail diagnostic:',
 paste('Selected KDE fits with provably infinite L2(P) risk:',sum(tailrows$kde_tail_margin<0),'of',nrow(tailrows)),
 'For these fits the finite-test RMSE is descriptive, not a consistent finite population-risk estimate.',
 'KDE empirical MCSE and bootstrap summaries are descriptive only where the fitted population risk is infinite; they are not confidence intervals for a finite population-risk parameter.',
 'Their raw results are preserved; no density floor or silent clipping was introduced.')
writeLines(text,file.path(out,'README.md'))
cat('Summarized',length(files),'jobs to',normalizePath(out),'\n')
