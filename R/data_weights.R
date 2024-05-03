data_weights <- function(data, X = NULL, strata = NULL, y, d, frailty = NULL, Random.Slope = NA,
                         ID = NULL, tr = NULL, weights = NULL, keep.names = TRUE){
  # does not consider Random slopes if !(X.tilde %in% X)

  if(!is.na(Random.Slope) & !Random.Slope%in%X){
    stop('function only works if random slope is in X.')
  }

  if(!is.null(frailty)){
    if(frailty %in% c(X,strata)) frailty <- NULL
  }
  if(!is.null(strata)){
    if(strata %in% X) strata <- NULL
  }

  if(keep.names){
    data.names <- c(ifelse(is.null(ID),'ID',ID),y,d,
                    ifelse(is.null(tr),'tr',tr))
    if(is.null(X)){
      data.names <- c(data.names,'zero')
    }else{
      data.names <- c(data.names,X)
    }
    data.names <- c(data.names,
                    ifelse(is.null(strata),'strata',strata),
                    ifelse(is.null(frailty),'frailty',frailty),
                    ifelse(is.null(weights),'weights',weights))
  }

  if(is.null(ID)){
    ID <- 1:dim(data)[1]
  }else{
    ID <- data[,ID]
  }
  data <- data[order(ID),]
  ID <- ID[order(ID)]

  # Retrieving/manipulating the data and plausability checks
  y <- data[,y]
  d <- data[,d]
  if(is.factor(d)){
    d <- as.numeric(d) - 1
  }
  if(is.null(tr)){
    tr <- y*0
  }else{
    tr <- data[,tr]
  }
  if(is.null(X)){
    X <- matrix(y*0, ncol = 1, dimnames = list(list(),list('zero')))
    K <- 1
  }else{
    K <- length(X)
    X <- data[,X, drop = FALSE]
    X <- as.matrix(X)
  }

  if(is.null(strata)){
    stratum <- factor(x = (y*0+1), levels = 1, labels = 'base')
  }else{
    stratum <- data[,strata]
  }

  if(is.null(frailty)){
    frailty2 <- factor(x = (y*0+1), levels = 1, labels = 'frail')
  }else{
    frailty2 <- data[,frailty]
  }

  # klassierte Daten (ber?cksichtigt cluster invariant weights)
  if(is.null(weights)){
    weights <- y*0+1
  }else{
    weights <- data[,weights]
  }

  ordr.command <- paste('order(ID','stratum',
                        paste('X','[,',1:dim(X)[2],']',sep='',collapse = ','),
                        'frailty2','y','tr','d','weights)',sep=',')
  idx.tmp <- eval(parse(text=ordr.command))
  max.var <- stats::aggregate(x = ID[idx.tmp], by = list(ID[idx.tmp]), FUN = length)[,-1]
  time.var <- unlist(sapply(max.var, FUN = function(x) 1:x,simplify = FALSE))

  long <- data.frame(ID[idx.tmp],y[idx.tmp],tr[idx.tmp],d[idx.tmp],weights[idx.tmp],stratum[idx.tmp],
                     frailty2[idx.tmp],X[idx.tmp,],time.var)
  colnames(long) <- c('ID','y','tr','d','weights','stratum','frailty',colnames(X),'time.var')
  wide <- stats::reshape(data = long, v.names = c('y','tr','d','stratum',colnames(X)), timevar = 'time.var',
                  idvar = c('ID','frailty','weights'),direction = 'wide')
  tmp <- unique(wide[,-c(1:2)], MARGIN = 1)
  n_xx <- apply(X = tmp, MARGIN = 1, FUN = function(x){
    apply(X = wide[,-c(1:2)], MARGIN = 1, FUN = function(z){
      all(z==x | (is.na(z)&is.na(x)))
    })
  })
  n_xx[is.na(n_xx)] <- FALSE
  n_xx <- apply(X = n_xx, MARGIN = 2, FUN = function(x) sum(wide$weights[x]))

  tmp <- cbind(1:dim(tmp)[1],tmp,n_xx)
  colnames(tmp)[c(1,dim(tmp)[2])] <- c('ID','weights')

  max.var <- max(max.var)
  Xnames <- matrix(NA,nrow=max.var,ncol=K)
  for(k in 1:K){
    Xnames[,k] <- paste(colnames(X)[k],1:max.var,sep='.')
  }
  col.tmp <- colnames(tmp)[-dim(tmp)[2]]
  Xnames <- t(cbind(col.tmp[seq(from=3,to=dim(wide)[2]-1,by=4+K)],
                    col.tmp[seq(from=4,to=dim(wide)[2]-1,by=4+K)],
                    col.tmp[seq(from=5,to=dim(wide)[2]-1,by=4+K)],
                    col.tmp[seq(from=6,to=dim(wide)[2]-1,by=4+K)],
                    Xnames))
  long <- reshape(tmp, varying=Xnames,
                  direction="long",  idvar=c("ID",'frailty','weights'),
                  v.names=c("y","tr", "d",'stratum',colnames(X)))
  long <- stats::na.omit(long)
  long <- long[order(long$ID),]
  max.var <- tapply(X = long$time, INDEX = long$ID, FUN = max)
  y <- as.vector(long$y)
  tr <- long$tr
  d <- as.vector(long$d)
  weights <- long$weights
  stratum <- long$stratum
  frailty2 <- long$frailty
  ID <- factor(x = long$ID, levels = 1:max(long$ID), labels = 1:max(long$ID))
  Xnames <- colnames(X)
  X <- as.matrix(long[,Xnames])
  colnames(X) <- Xnames

  data <- data.frame(ID = ID, y = y, d = d, tr = tr, X = X,
                     stratum = stratum, frailty = frailty2, weights = weights)
  if(keep.names){
    colnames(data) <- data.names
  }

  return(data)

}
