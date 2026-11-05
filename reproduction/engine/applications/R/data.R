# Adapters return numeric predictors, outcomes, IDs, sites and dependence blocks.
app_dataset <- function(name) {
  raw <- file.path('applications', 'data', 'raw')
  if (name == 'heart') {
    sites <- c('cleveland', 'hungarian', 'switzerland', 'va')
    parts <- lapply(sites, function(s) {
      d <- read.csv(file.path(raw, 'heart', paste0('processed.', s, '.data')),
                    header=FALSE, na.strings='?')
      names(d) <- c('age','sex','cp','bp','chol','fbs','ecg','hr','exang',
                    'oldpeak','slope','ca','thal','disease')
      # Zero is a missing-value sentinel for these physiological measurements.
      for(k in c('bp','chol','hr')) d[[k]][d[[k]]==0] <- NA_real_
      # Fixed, documented category levels; no test-data-derived encoding.
      x <- cbind(as.matrix(d[c('age','sex','bp','chol','hr','exang','oldpeak')]),
                 sapply(2:4, function(k) as.numeric(d$cp == k)))
      list(x=x, y=as.numeric(d$disease > 0), site=rep(s,nrow(d)),
           id=paste0(s,':',seq_len(nrow(d))))
    })
    z <- list(x=do.call(rbind,lapply(parts,`[[`,'x')),
              y=unlist(lapply(parts,`[[`,'y')),
              site=unlist(lapply(parts,`[[`,'site')),
              id=unlist(lapply(parts,`[[`,'id')), type='binary')
    z$block <- z$id
  } else if (name == 'beijing') {
    cache <- 'applications/data/beijing-daily.rds'
    if (file.exists(cache)) { z <- readRDS(cache); z$name <- name; return(z) }
    files <- list.files(file.path(raw,'beijing'), 'PRSA_Data.*csv$',full.names=TRUE)
    stopifnot(length(files)==12)
    parts <- lapply(files,function(f) {
      d <- read.csv(f,check.names=FALSE)
      d$date <- as.Date(sprintf('%04d-%02d-%02d',d$year,d$month,d$day))
      vars <- c('PM2.5','TEMP','PRES','DEWP','RAIN','WSPM')
      # Require at least 18 measured hours per variable per daily mean.
      groups <- split(seq_len(nrow(d)),d$date)
      a <- t(vapply(groups,function(ii) vapply(d[ii,vars],function(v)
        if(sum(is.finite(v))>=18) mean(v,na.rm=TRUE) else NA_real_,numeric(1)),numeric(6)))
      date <- as.Date(names(groups)); season <- as.numeric(format(date,'%j'))
      list(x=cbind(a[,-1,drop=FALSE],sin(2*pi*season/365.25),cos(2*pi*season/365.25)),
           y=a[,1],site=rep(d$station[1],nrow(a)),id=paste(d$station[1],date,sep=':'),
           block=format(date,'%Y-%m'))
    })
    z <- lapply(c('x','y','site','id','block'),function(k)
      if(k=='x') do.call(rbind,lapply(parts,`[[`,k)) else unlist(lapply(parts,`[[`,k)))
    names(z) <- c('x','y','site','id','block'); z$type <- 'regression'
    # Calendar blocks shared across stations, with chronological holdouts.
    z$chronological <- TRUE
    keep <- is.finite(z$y)
    for(k in c('y','site','id','block')) z[[k]] <- z[[k]][keep]
    z$x <- z$x[keep,,drop=FALSE]
    saveRDS(z,cache)
  } else if(name == 'hydraulic') {
    cache <- 'applications/data/hydraulic-features.rds'
    if(file.exists(cache)) { z <- readRDS(cache); z$name <- name; return(z) }
    sensors <- c('FS1','FS2','TS1','TS2','TS3','TS4','VS1')
    features <- lapply(sensors,function(s) {
      a <- as.matrix(read.table(file.path(raw,'hydraulic',paste0(s,'.txt'))))
      cbind(mean=rowMeans(a),sd=apply(a,1,sd))
    })
    profile <- read.table(file.path(raw,'hydraulic','profile.txt'))
    z <- list(x=do.call(cbind,features),y=as.numeric(profile[,1]!=100),
              site=rep('test-rig',nrow(profile)),id=paste0('cycle:',seq_len(nrow(profile))),
              block=as.character(ceiling(seq_len(nrow(profile))/10)),type='anomaly')
    z$label_definition <- 'Cooler condition below 100%; other component conditions are nuisance variation.'
    saveRDS(z,cache)
  } else {
    if(name=='spambase') {
      d <- read.csv('applications/data/raw/spambase/spambase.data',header=FALSE)
      x <- as.matrix(d[,-58]); y <- d[,58]; type <- 'binary'
    } else if(name=='airfoil') {
      d <- read.table('DistKRE/conformal/tibshirani2019/airfoil.txt')
      x <- as.matrix(d[,1:5]); x[,c(1,5)] <- log(x[,c(1,5)])
      y <- d[,6]; type <- 'regression'
    } else if(name=='qsar') {
      d <- read.csv('DistKRE/qsar_aquatic_toxicity.csv',sep=';',header=FALSE)
      x <- as.matrix(d[,1:8]); y <- d[,9]; type <- 'regression'
    } else if(name=='banknote') {
      d <- read.csv('DistKRE/data_banknote_authentication.txt',header=FALSE)
      x <- as.matrix(d[,1:4]); y <- d[,5]; type <- 'anomaly'
      # Operational benchmark convention: class 0 is the reference class.
    } else if(name=='breast') {
      d <- read.csv('DistKRE/breast_cancer_bd.csv',na.strings='?')
      x <- as.matrix(d[,2:10]); storage.mode(x) <- 'double'
      y <- as.numeric(d[,11]==4); type <- 'anomaly'
    } else stop('Unknown dataset: ',name)
    z <- list(x=x,y=y,site=rep('constructed',length(y)),
              id=paste0(name,':',seq_along(y)),type=type)
    # Duplicate observations never cross split roles.
    z$block <- apply(cbind(x,y),1,paste,collapse='|')
  }
  z$name <- name
  stopifnot(nrow(z$x)==length(z$y),!anyDuplicated(z$id),all(is.finite(z$y)))
  z
}

app_roles <- function(d,seed) {
  set.seed(seed)
  blocks <- sort(unique(d$block))
  if(!isTRUE(d$chronological)) blocks <- sample(blocks)
  # Disjoint role blocks, shared across all sites. Target fitting uses ratio
  # blocks only; target evaluation uses test blocks only.
  role <- cut(seq_along(blocks)/length(blocks),c(0,.30,.60,.80,1),
              labels=c('ratio','model','calibration','test'))
  as.character(role[match(d$block,blocks)])
}

app_preprocess <- function(clients) {
  count <- Reduce(`+`,lapply(clients,function(x) colSums(is.finite(x))))
  total <- Reduce(`+`,lapply(clients,function(x) colSums(replace(x,!is.finite(x),0))))
  square <- Reduce(`+`,lapply(clients,function(x) colSums(replace(x,!is.finite(x),0)^2)))
  mu <- total/pmax(1,count)
  sd <- sqrt(pmax(0,square/pmax(1,count)-mu^2)); sd[sd<1e-8] <- 1
  list(mean=mu,sd=sd,observed=count)
}
app_transform <- function(x,p) {
  for(j in seq_len(ncol(x))) x[!is.finite(x[,j]),j] <- p$mean[j]
  sweep(sweep(x,2,p$mean),2,p$sd,'/')
}
