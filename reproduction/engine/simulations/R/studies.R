sim_bind <- function(rows) {
  if(!length(rows))return(data.frame())
  names_all<-unique(unlist(lapply(rows,names)))
  do.call(rbind,lapply(rows,function(x){for(n in setdiff(names_all,names(x)))x[[n]]<-NA;x[names_all]}))
}
sim_row <- function(method,metrics=NULL,status='ok',error='',R=NA_integer_,lambda=NA_real_,
  sigma=NA_real_,iterations=NA_integer_,validation=NA_real_,fit_coordinates=NA_real_,
  tuning_coordinates=0,seconds=NA_real_,row_type='estimate',budget=NA_real_,
  optimization_gap=NA_real_,coefficient_gap=NA_real_,prediction_gap=NA_real_,boundary=NA) {
  if(is.null(metrics))metrics<-list(mse=NA_real_,rmse=NA_real_,calibration=NA_real_,negative=NA_real_)
  as.data.frame(c(list(method=method,status=status,error=error,R=R,lambda=lambda,
    sigma=sigma,iterations=iterations,validation=validation,fit_coordinates=fit_coordinates,
    tuning_coordinates=tuning_coordinates,total_coordinates=fit_coordinates+tuning_coordinates,
    seconds=seconds,row_type=row_type,budget=budget,optimization_gap=optimization_gap,
    coefficient_gap=coefficient_gap,prediction_gap=prediction_gap,boundary=boundary),metrics),stringsAsFactors=FALSE)
}
sim_timed <- function(expr) {
  start<-proc.time()[['elapsed']];value<-force(expr)
  list(value=value,seconds=proc.time()[['elapsed']]-start)
}
sim_rf_study <- function(data,case,config,tuned=TRUE) {
  C<-length(data$xs);R<-config$features
  sigmas<-data$sigma*if(tuned)config$bandwidth_multipliers else 1
  lambdas<-if(tuned)config$lambdas else .01
  chosen_closed<-NULL;table<-list();tune_gd<-0;tune_closed<-0;init_time<-0
  for(sigma in sigmas) {
    tmp<-sim_timed(sim_rf_cache(data,R,sigma));cache<-tmp$value;init_time<-init_time+tmp$seconds
    for(lambda in lambdas) {
      cl<-sim_timed(.ridge_coefficients(cache$H0,cache$a,lambda))
      tune_closed<-tune_closed+cl$seconds
      value<-cache$validate(cl$value$coefficients)
      row<-data.frame(sigma=sigma,lambda=lambda,validation=value)
      if(is.null(chosen_closed)||sim_select(rbind(chosen_closed$row,row))==2L)
        chosen_closed<-list(row=row,b=cl$value$coefficients,spec=cache$spec)
      if(tuned) {
        gd<-sim_timed(sim_gd(cache$a,data$pi,cache$exchange,lambda,config$tuning_iterations))
        tune_gd<-tune_gd+gd$seconds
        table[[length(table)+1]]<-data.frame(sigma=sigma,lambda=lambda,
          validation=cache$validate(gd$value$coefficients),closed_validation=value)
      }
    }
    rm(cache);gc(FALSE)
  }
  pred_closed<-sim_rf_predict(data$test,chosen_closed$spec,chosen_closed$b,config$batch_size)[,1]
  nk<-length(sigmas)*length(lambdas)
  closed_tuning<-if(tuned)length(sigmas)*C*R*(R+1)/2+nk*(C*R+C) else 0
  closed_fit<-if(tuned)0 else C*R*(R+1)/2
  boundary<-function(lambda,sigma) lambda%in%range(config$lambdas)||sigma%in%range(sigmas)
  rows<-list(sim_row('rff_closed',sim_metrics(pred_closed,data$truth,data$tw),
    R=R,lambda=chosen_closed$row$lambda,sigma=chosen_closed$row$sigma,
    validation=chosen_closed$row$validation,fit_coordinates=closed_fit,
    tuning_coordinates=closed_tuning,seconds=init_time+tune_closed,
    boundary=if(tuned)boundary(chosen_closed$row$lambda,chosen_closed$row$sigma) else FALSE))
  if(tuned) {
    tab<-do.call(rbind,table);selected<-tab[sim_select(tab),]
  } else {tab<-data.frame();selected<-data.frame(sigma=data$sigma,lambda=.01)}
  tmp<-sim_timed(sim_rf_cache(data,R,selected$sigma));cache<-tmp$value
  T<-if(tuned)config$final_iterations else 1000L
  gd<-sim_timed(sim_gd(cache$a,data$pi,cache$exchange,selected$lambda,T))
  bstar<-.ridge_coefficients(cache$H0,cache$a,selected$lambda)$coefficients
  B<-cbind(gd$value$coefficients,bstar)
  pred<-sim_rf_predict(data$test,cache$spec,B,config$batch_size)
  delta<-gd$value$coefficients-bstar
  gap<-as.numeric(crossprod(delta,(cache$H0+diag(selected$lambda,R))%*%delta))/2
  tuning_coords<-if(tuned)nk*(2*C*R*config$tuning_iterations+C*R+C) else 0
  rows[[length(rows)+1]]<-sim_row('distributed_gd',sim_metrics(pred[,1],data$truth,data$tw),
    R=R,lambda=selected$lambda,sigma=selected$sigma,iterations=T,
    validation=cache$validate(gd$value$coefficients),fit_coordinates=gd$value$coordinates,
    tuning_coordinates=tuning_coords,seconds=init_time+tune_gd+tmp$seconds+gd$seconds,
    optimization_gap=gap,coefficient_gap=sqrt(sum(delta^2)),
    prediction_gap=sqrt(sum(data$tw*(pred[,1]-pred[,2])^2)),
    boundary=if(tuned)boundary(selected$lambda,selected$sigma) else FALSE)
  rows[[length(rows)+1]]<-sim_row('rff_at_gd_parameters',sim_metrics(pred[,2],data$truth,data$tw),
    R=R,lambda=selected$lambda,sigma=selected$sigma,validation=cache$validate(bstar),
    fit_coordinates=C*R*(R+1)/2,row_type='diagnostic_reference')
  if(case$study=='imbalance') {
    ab<-sim_timed(sim_gd(cache$a,rep(1/C,C),cache$exchange,selected$lambda,T))
    p<-sim_rf_predict(data$test,cache$spec,ab$value$coefficients,config$batch_size)[,1]
    rows[[length(rows)+1]]<-sim_row('gd_uniform_weight_ablation',sim_metrics(p,data$truth,data$tw),
      R=R,lambda=selected$lambda,sigma=selected$sigma,iterations=T,
      fit_coordinates=ab$value$coordinates,seconds=ab$seconds,
      validation=cache$validate(ab$value$coefficients))
  }
  list(rows=sim_bind(rows),tuning=tab,feature_initialization_coordinates=C*5*length(sigmas),sample_count_coordinates=C)
}
sim_accuracy <- function(data,case,config) {
  rows<-list();tuning<-list();tuned<-case$protocol=='tuned'
  safe<-function(method,fn) {
    start<-proc.time()[['elapsed']]
    tryCatch(fn(),error=function(e) {
      sim_row(method,status='failed',error=conditionMessage(e),seconds=proc.time()[['elapsed']]-start)
    })
  }
  rf<-tryCatch(sim_rf_study(data,case,config,tuned),error=function(e) {
    list(rows=sim_bind(lapply(c('rff_closed','distributed_gd'),function(m)
      sim_row(m,status='failed',error=conditionMessage(e)))),tuning=data.frame())
  })
  rows[[length(rows)+1]]<-rf$rows;tuning$rf<-rf$tuning
  if(case$N<=config$exact_kernel_cap)rows[[length(rows)+1]]<-safe('kernel_pooled',function() {
    start<-proc.time()[['elapsed']]
    fit<-sim_kernel_fit(data$xs,data$z,data$xv,data$zv,data$pi,data$sigma,config,tuned)
    p<-sim_kernel_predict(data$test,fit$cache,fit$alpha,fit$lambda,config$batch_size)
    sim_row('kernel_pooled',sim_metrics(p,data$truth,data$tw),lambda=fit$lambda,
      sigma=fit$sigma,validation=fit$validation,seconds=proc.time()[['elapsed']]-start,
      boundary=if(tuned)fit$lambda%in%range(config$lambdas)||fit$sigma%in%range(data$sigma*config$bandwidth_multipliers) else FALSE)
  })
  gc(FALSE)
  rows[[length(rows)+1]]<-safe('kde_pooled',function() {
    start<-proc.time()[['elapsed']];fit<-sim_kde_fit_predict(data,config,tuned)
    row<-sim_row('kde_pooled',sim_metrics(fit$pred,data$truth,data$tw),
      sigma=fit$source_bandwidth,validation=fit$validation,seconds=proc.time()[['elapsed']]-start)
    row$target_bandwidth<-fit$target_bandwidth;row
  })
  gc(FALSE)
  rows[[length(rows)+1]]<-safe('kernel_local_naive',function() {
    start<-proc.time()[['elapsed']];p<-numeric(nrow(data$test));settings<-list()
    for(c in seq_along(data$xs)) {
      fit<-sim_kernel_fit(list(data$xs[[c]]),data$z,list(data$xv[[c]]),data$zv,1,data$sigma,config,tuned)
      ix<-which(data$test_client==c)
      p[ix]<-sim_kernel_predict(data$test[ix,,drop=FALSE],fit$cache,fit$alpha,fit$lambda,config$batch_size)
      settings[[c]]<-c(lambda=fit$lambda,sigma=fit$sigma)
      rm(fit);gc(FALSE)
    }
    row<-sim_row('kernel_local_naive',sim_metrics(p,data$truth,data$tw),seconds=proc.time()[['elapsed']]-start)
    row$local_settings<-paste(vapply(settings,function(s)paste(format(s,digits=8),collapse=','),character(1)),collapse=';')
    row
  })
  rows[[length(rows)+1]]<-sim_row('constant_one',sim_metrics(rep(1,nrow(data$test)),data$truth,data$tw),seconds=0)
  list(rows=sim_bind(rows),tuning=tuning,
    metadata=list(feature_initialization_coordinates=rf$feature_initialization_coordinates,
                  sample_count_coordinates=length(data$xs)))
}
sim_communication <- function(data,case,config) {
  rows<-list();C<-length(data$xs)
  for(R in config$mechanism_features) {
    init<-sim_timed(sim_rf_cache(data,R,data$sigma));cache<-init$value
    for(lambda in config$mechanism_lambdas) {
      budget_T<-floor(config$communication_budgets/(2*C*R))
      cp<-sort(unique(c(config$checkpoints,budget_T)))
      cp<-cp[cp<=max(config$checkpoints)]
      ref<-sim_timed(.ridge_coefficients(cache$H0,cache$a,lambda))
      gd<-sim_timed(sim_gd(cache$a,data$pi,cache$exchange,lambda,max(cp),cp))
      pred<-sim_rf_predict(data$test,cache$spec,cbind(gd$value$path,ref$value$coefficients),config$batch_size)
      reference<-pred[,ncol(pred)]
      for(i in seq_along(cp)) {
        delta<-gd$value$path[,i]-ref$value$coefficients
        gap<-as.numeric(crossprod(delta,(cache$H0+diag(lambda,R))%*%delta))/2
        base<-sim_row('distributed_gd',sim_metrics(pred[,i],data$truth,data$tw),R=R,
          lambda=lambda,sigma=data$sigma,iterations=cp[i],fit_coordinates=2*C*R*cp[i],
          optimization_gap=gap,coefficient_gap=sqrt(sum(delta^2)),
          prediction_gap=sqrt(sum(data$tw*(pred[,i]-reference)^2)),
          seconds=unname(gd$value$checkpoint_seconds[as.character(cp[i])]),row_type='checkpoint')
        base$seconds_measurement<-'measured_at_checkpoint'
        if(cp[i]%in%config$checkpoints)rows[[length(rows)+1]]<-base
        for(b in which(budget_T==cp[i])) {
          budgetrow<-base;budgetrow$row_type<-'budget';budgetrow$budget<-config$communication_budgets[b]
          rows[[length(rows)+1]]<-budgetrow
        }
      }
      rows[[length(rows)+1]]<-sim_row('rff_closed',sim_metrics(reference,data$truth,data$tw),
        R=R,lambda=lambda,sigma=data$sigma,fit_coordinates=C*R*(R+1)/2,
        seconds=ref$seconds+init$seconds,row_type='reference')
    }
    rm(cache,pred);gc(FALSE)
  }
  list(rows=sim_bind(rows),tuning=list())
}
sim_run_replication <- function(case,replication,phase,config,source_hash='') {
  start<-proc.time()[['elapsed']]
  data<-sim_data(case,replication,phase,config)
  result<-if(case$protocol=='path')sim_communication(data,case,config) else sim_accuracy(data,case,config)
  result$rows$case_id<-case$id;result$rows$study<-case$study
  result$rows$N<-case$N;result$rows$M<-case$M;result$rows$tau<-case$tau
  result$rows$replication<-replication;result$rows$phase<-phase
  result$metadata<-c(result$metadata,list(config=config,case=case,n=data$nf,m=nrow(data$z),
    source_weights=data$pi,source_hash=source_hash,elapsed=proc.time()[['elapsed']]-start,
    seeds=vapply(1:4,function(k)sim_seed(config,case,replication,phase,k),integer(1))))
  result$status<-if(any(result$rows$status=='failed'))'method_failure' else 'ok'
  result
}
