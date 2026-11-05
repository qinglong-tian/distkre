natural_data <- function(name,outcome=NULL) {
  if(name %in% c('heart','beijing')) return(app_dataset(name))
  fs<-list.files('applications/data/raw/turbine','gt_.*csv$',full.names=TRUE)
  z<-lapply(fs,function(f) {
    a<-read.csv(f); year<-sub('.*gt_([0-9]+)[.]csv','\\1',f)
    list(x=as.matrix(a[,c('AT','AP','AH','AFDP','GTEP','TIT','TAT','TEY','CDP')]),
         y=a[[outcome]],site=rep(year,nrow(a)),id=paste(year,seq_len(nrow(a)),sep=':'),
         block=as.character(ceiling(seq_len(nrow(a))/48)))
  })
  list(x=do.call(rbind,lapply(z,`[[`,'x')),y=unlist(lapply(z,`[[`,'y')),
       site=unlist(lapply(z,`[[`,'site')),id=unlist(lapply(z,`[[`,'id')),
       block=unlist(lapply(z,`[[`,'block')),type='regression',name=paste0('turbine_',outcome))
}

natural_roles <- function(d,seed) {
  if(d$name=='heart') {
    set.seed(seed); roles<-rep('',length(d$y))
    for(site in unique(d$site)) {
      ii<-which(d$site==site); ii<-sample(ii)
      roles[ii]<-as.character(cut(seq_along(ii)/length(ii),c(0,.25,.5,.65,1),
        labels=c('model','ratio','validation','evaluation')))
    }
    return(roles)
  }
  if(d$name=='beijing') {
    months<-sort(unique(d$block))
    r<-as.character(cut(seq_along(months)/length(months),c(0,.25,.6,.7,1),
      labels=c('model','ratio','validation','evaluation')))
    # Exclude a complete month at each role boundary.
    boundaries<-which(c(FALSE,r[-1]!=r[-length(r)])); r[boundaries]<-'gap'
    return(r[match(d$block,months)])
  }
  roles<-rep('',length(d$y))
  for(site in unique(d$site)) {
    ii<-which(d$site==site); blocks<-unique(d$block[ii])
    r<-as.character(cut(seq_along(blocks)/length(blocks),c(0,.25,.6,.7,1),
       labels=c('model','ratio','validation','evaluation')))
    boundaries<-which(c(FALSE,r[-1]!=r[-length(r)])); r[boundaries]<-'gap'
    roles[ii]<-r[match(d$block[ii],blocks)]
  }
  roles
}

natural_split <- function(d,target,seed,cap=80L) {
  roles<-natural_roles(d,seed); set.seed(seed+1000L)
  sources<-if(grepl('turbine',d$name)) c('2011','2012','2013') else setdiff(sort(unique(d$site)),target)
  take<-function(site,role,n) app_take(which(d$site==site & roles==role),n)
  split<-list(fit=lapply(sources,function(s) take(s,'ratio',cap)),
    val=lapply(sources,function(s) take(s,'validation',cap)),
    model=lapply(sources,function(s) take(s,'model',cap*2L)),
    evaluation=lapply(sources,function(s) take(s,'evaluation',cap*2L)),
    q=take(target,'ratio',3L*cap),qval=take(target,'validation',cap),
    test=take(target,'evaluation',6L*cap),sources=sources)
  groups<-c(split$fit,split$val,split$model,split$evaluation,list(split$q,split$qval,split$test))
  stopifnot(!anyDuplicated(unlist(groups)),min(lengths(groups))>=5)
  split$manifest<-data.frame(id=d$id,site=d$site,block=d$block,role='unused')
  for(key in c('fit','val','model','evaluation')) for(ii in split[[key]]) split$manifest$role[ii]<-key
  for(key in c('q','qval','test')) split$manifest$role[split[[key]]]<-key
  split
}
