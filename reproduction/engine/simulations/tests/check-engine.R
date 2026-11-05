for(f in list.files('R',full.names=TRUE))source(f)
source('simulations/configs/study.R')
source('simulations/R/core.R');source('simulations/R/studies.R')
config<-simulation_config(); config$features<-12L;config$mechanism_features<-c(8L,12L)
config$test_points<-200L;config$lambdas<-.01;config$bandwidth_multipliers<-1
config$tuning_iterations<-10L;config$final_iterations<-20L
config$checkpoints<-c(0L,1L,2L,10L,20L);config$mechanism_lambdas<-.01
config$communication_budgets<-c(128,640)
case<-simulation_cases(config)[[1]];case$N<-100;case$M<-100
D<-sim_data(case,1,'pilot',config)
stopifnot(abs(sum(D$tw)-1)<1e-14,all(D$truth>0),all(D$truth<=2),
  identical(D,sim_data(case,1,'pilot',config)))
null<-config;null$rho<-case$pi
stopifnot(max(abs(sim_truth(D$test,case,null)-1))<1e-14)
sp8<-sim_spec(D,8,D$sigma);sp12<-sim_spec(D,12,D$sigma)
stopifnot(max(abs(evaluate_rff(D$test,sp8)*sqrt(8)-evaluate_rff(D$test,sp12)[,1:8]*sqrt(12)))<1e-12)
cache<-sim_rf_cache(D,12,D$sigma)
pub<-dkre_server_fit(D$z,D$nf,.01,cache$spec,cache$exchange,iterations=20)
path<-sim_gd(cache$a,D$pi,cache$exchange,.01,20,c(0,1,2,10,20))
stopifnot(max(abs(pub$metadata$coefficients-path$coefficients))<1e-12,
  max(abs(path$path[,'1']-cache$a/(2+.01)))<1e-12)
b<-path$coefficients
sourceval<-sum(D$pi*vapply(D$xv,function(x)mean(as.numeric(evaluate_rff(x,cache$spec)%*%b)^2),numeric(1)))/2
rawval<-sourceval-mean(as.numeric(evaluate_rff(D$zv,cache$spec)%*%b))
stopifnot(abs(cache$validate(b)-rawval)<1e-12)
fit<-sim_kernel_fit(D$xs,D$z,D$xv,D$zv,D$pi,D$sigma,config,FALSE)
reference<-fit_kulsif(do.call(rbind,D$xs),D$z,.01,D$sigma)
p1<-sim_kernel_predict(D$test,fit$cache,fit$alpha,.01,17)
stopifnot(max(abs(p1-predict(reference,D$test)))<1e-10)
kde<-sim_kde_fit_predict(D,config,FALSE)
refkde<-fit_kde_ratio(D$xs,D$z,clip=FALSE)
stopifnot(max(abs(kde$pred-predict(refkde,D$test)))<1e-9)
# Partition invariance is checked with unequal client sizes.
X<-do.call(rbind,D$xs)
for(part in list(list(seq_len(nrow(X))),split(seq_len(nrow(X)),c(rep(1,10),rep(2,70))))) {
 clients<-lapply(part,function(ix)dkre_client(X[ix,,drop=FALSE],cache$spec))
 exchange<-function(b)lapply(clients,dkre_client_update,coefficients=b)
 sizes<-vapply(clients,function(cl)cl$n_source,integer(1))
 re<-sim_gd(cache$a,sizes/sum(sizes),exchange,.01,20)
 stopifnot(max(abs(re$coefficients-path$coefficients))<1e-12)
}
res<-sim_run_replication(case,1,'pilot',config)
stopifnot(all(res$rows$status=='ok'),all(is.finite(res$rows$rmse)))
communication<-case;communication$protocol<-'path';communication$study<-'communication'
cres<-sim_run_replication(communication,1,'pilot',config)
stopifnot(all(cres$rows$status=='ok'),all(is.finite(cres$rows$rmse)),
  all(cres$rows$optimization_gap[cres$rows$row_type=='checkpoint']>=0))
# The same seeded job gives identical statistical results despite timing noise.
res2<-sim_run_replication(case,1,'pilot',config)
stopifnot(identical(res$rows$rmse,res2$rows$rmse))
cat('Simulation engine checks passed: truth, RNG, shared features, GD, kernel, KDE, validation, partitions and checkpoints.\n')
