
maxll <- function(start_dist, start_lnlam, start_beta, c.weights, m.weights, overdisp, n.date,
                  baseline, type, data, level,y, d, X, ID, current.status, slope.model,H.singular,
                  sd_x, tr,stratum, breaks,frailty, K, N, slope, X.tilde,start_slope,start_pi,
                  iterlim, converge, multi.num, cs.options, frail.thresh, c.options, m.options, R.options,
                  scale.marginal, scale.cond, z, v, initial.z, prnt, num.score, num.hessian){

  n <- length(y)
  n_C <- length(unique(ID))
  n_i <- tapply(X = ID, INDEX = ID, FUN = length)
  max.n_i <- max(n_i)
  multi.num <- ifelse(!multi.num, max.n_i > 2, multi.num)

  strata <- levels(stratum)
  Q <- length(strata)
  if( Q != length(baseline[[1]]) ){
    baseline[[1]] <- rep(baseline[[1]], Q)
  }

  baseline <- c( baseline, replicate(NULL, n = 4 - length(baseline)) )
  names(baseline[[1]]) <- strata
  if(any(baseline[[1]] == 'Breslow')){
    baseline[2] <- 'Breslow'
    if(any(baseline[[1]] != 'Breslow')){
      stop('Breslow cannot be mixed with piecewise or GG. Maybe choose piecewise with Breslow breakpoints.')
    }
  }
  if(any(baseline[[1]] == 'location')){
    if(is.null(baseline[[4]])) {
      baseline[[4]] <- which(baseline[[1]] == 'location')
    }
    if(!all( which(baseline[[1]] == 'location') %in% sort(unlist(baseline[[4]])) |
             length(unlist(baseline[[4]])) == sum(baseline[[1]] == 'location') )){
      stop('Loacation stratification is incorrectly specified. Either one or more are missing
            or there are too many or some occur more than once in baseline[[4]]!')
    }
  }

  # set breaks for (piecewise) hazard and indices for hazards
  partial.rank <- GG.rank <- breslow.rank <- rep(NA, Q)
  partial.idx <- baseline[[1]] %in% c('piecewise','Breslow')
  partial.rank[partial.idx] <- rank(which(partial.idx))
  partial_i <- stratum %in% strata[which(partial.idx)]
  partial.stratum <- partial.rank[stratum[partial_i]] #individual level. Indicator of stratum, excluding all non-partial strata
  bline <- factor(x = rep(1,Q), levels = 1:16,
                  labels = c('piecewise', 'Breslow', 'GG', 'gamma', 'gamma^-1', 'Weibull', 'Weibull^-1',
                             'exponential', 'exponential^-1', 'ammag', 'ammag^-1', 'lognormal', 'location', 'N/2', 'GGpos', 'GGinv') )
  bline[] <- baseline[[1]]
  GG.idx <- bline != 'piecewise' & bline != 'Breslow'
  GG.special <- any(bline != 'piecewise' & bline != 'GG' & bline != 'Breslow')
  GG.rank[GG.idx] <- rank(which(GG.idx))
  GG_i <- stratum %in% strata[which(GG.idx)]
  GG.stratum <- GG.rank[stratum[GG_i]] #individual level. Indicator of stratum, excluding all non-GG strata
  breslow.idx <- baseline[[1]] %in% 'Breslow'
  breslow.rank[breslow.idx] <- rank(which(breslow.idx))
  breslow_i <- (stratum %in% strata[which(breslow.idx)])
  names(breslow.idx) <- names(partial.idx) <-  names(GG.idx) <-
    names(breslow.rank) <-  names(partial.rank) <- names(GG.rank) <-  strata
  n.p <- sum(partial_i)
  n.GG <- sum(GG_i)
  Q.partial <- Q - sum(GG.idx)
  Q.GG <- Q - Q.partial
  p.strata <- strata[partial.idx]
  GG.strata <- strata[GG.idx]
  if(is.null(breaks)){
    breaks <- replicate(list(0:3), n = Q) # placeholder for parametric hazard
    breaks[partial.idx] <- breaks_piece(y = y[partial_i], d = d[partial_i], strata = stratum[partial_i],
                                        Qstar = Q.partial, stratum = strata[partial.idx],
                                        style = baseline[[2]], by = baseline[[3]])
  }
  names(breaks) <- strata
  numh <- sapply(X = breaks, FUN = function(x) length(x)-1)

  deleted <- rep(0,Q)
  GG.op <- c()
  GG.sc <- 0
  sc.inv <- TRUE
  if(GG.special){
    # GG.sc: holt Elemente aus optimax Vektor und f?llt sie in GG Vektor
    # sc.inv: holt die Elemente aus dem GG Vektor und packt sie in den zu sch?tzenden Vektor
    # GG.op: f?hrt die Operationen aus (z.B. ^-1) die f?r GG Parameter n?tig sind
    # op.inv: inverse Operationen f?r GG.op. Nur f?r relevante, also zu sch?tzende.
    ## Alle anderen bleiben unver?ndert.
    ## piecewise bleibt hier grunds?tzlich unbeachtet!
    sc.inv <- 0
    add <- 0
    op.inv <- c()
    Q.loc <- length(baseline[4])
    loc.idx <- matrix(NA, nrow = Q.loc, ncol = 2)
    loc.1st <- rep(TRUE,Q.loc)
    for(q in 1:Q){
      if(bline[q] == 'piecewise') next
      mG <- max(GG.sc)
      mSC <- max(sc.inv) + add
      add <- deleted[q] <- 1
      if(bline[q] == 'GG'){
        GG.sc <- c(GG.sc,mG+1:3)
        sc.inv <- c(sc.inv,mSC+1:3)
        GG.op <- c(GG.op,'exp','identity','identity')
        op.inv <- c(op.inv,'log','identity','identity')
        add <- deleted[q] <- 0
      }else if(bline[q] == 'gamma'){
        GG.sc <- c(GG.sc,mG + c(1,2,1))
        sc.inv <- c(sc.inv,mSC + 1:2)
        GG.op <- c(GG.op,'exp','identity','exp')
        op.inv <- c(op.inv,'log','identity','identity')
      }else if(bline[q] == 'gamma^-1'){
        GG.sc <- c(GG.sc,mG+c(1,2,1))
        sc.inv <- c(sc.inv,mSC+1:2)
        GG.op <- c(GG.op,'exp','identity','negexp')
        op.inv <- c(op.inv,'log','identity','identity')
      }else if(bline[q] == 'Weibull'){
        GG.sc <- c(GG.sc,mG+c(1,2,1))
        sc.inv <- c(sc.inv,mSC+1:2)
        GG.op <- c(GG.op,'exp','identity','one')
        op.inv <- c(op.inv,'log','identity','identity')
      }else if(bline[q] == 'Weibull^-1'){
        GG.sc <- c(GG.sc,mG+c(1,2,1))
        sc.inv <- c(sc.inv,mSC+1:2)
        GG.op <- c(GG.op,'exp','identity','negone')
        op.inv <- c(op.inv,'log','identity','identity')
      }else if(bline[q] == 'exponential'){
        GG.sc <- c(GG.sc,mG+rep(1,3))
        sc.inv <- c(sc.inv,mSC+2)
        GG.op <- c(GG.op,'one','identity','one')
        op.inv <- c(op.inv,'identity','identity','identity')
        deleted[q] <- 2
      }else if(bline[q] == 'exponential^-1'){
        GG.sc <- c(GG.sc,mG+rep(1,3))
        sc.inv <- c(sc.inv,mSC+2)
        GG.op <- c(GG.op,'negone','identity','one')
        op.inv <- c(op.inv,'identity','identity','identity')
        deleted[q] <- 2
      }else if(bline[q] == 'ammag'){
        GG.sc <- c(GG.sc,mG+c(1,2,1))
        sc.inv <- c(sc.inv,mSC+1:2)
        GG.op <- c(GG.op,'exp','identity','invexp')
        op.inv <- c(op.inv,'log','identity','identity')
      }else if(bline[q] == 'ammag^-1'){
        GG.sc <- c(GG.sc,mG+c(1,2,1))
        sc.inv <- c(sc.inv,mSC+1:2)
        GG.op <- c(GG.op,'exp','identity','neginvexp')
        op.inv <- c(op.inv,'log','identity','identity')
      }else if(bline[q] == 'lognormal'){
        GG.sc <- c(GG.sc,mG+c(1,2,1))
        sc.inv <- c(sc.inv,mSC+1:2)
        GG.op <- c(GG.op,'exp','identity','zero')
        op.inv <- c(op.inv,'log','identity','identity')
      }else if(bline[q] == 'location'){
        for(q2 in 1:Q.loc){
          if(q%in%baseline[[4]][q2]) strat.group <- q2
        }
        if(loc.1st[strat.group]){
          loc.idx[strat.group,] <- mG+c(1,3)
          deleted[q] <- add <- 0
          sc.inv <- c(sc.inv,mSC+1:3)
        }else{
          deleted[q] <- 2
          sc.inv <- c(sc.inv, mSC+2)
        }
        GG.sc <- c(GG.sc,loc.idx[strat.group,1], mG+2-1*!loc.1st , loc.idx[strat.group,2])
        GG.op <- c(GG.op,'exp','identitiy','identity')
        op.inv <- c(op.inv,'log','identity','identity')
        loc.1st[strat.group] <- FALSE
      }else if(bline[q] == 'GGpos'){
        GG.sc <- c(GG.sc,mG+1:3)
        sc.inv <- c(sc.inv,mSC+1:3)
        GG.op <- c(GG.op,'exp','identity','exp')
        op.inv <- c(op.inv,'log','identity','log')
        add <- deleted[q] <- 0
      }else if(bline[q] == 'GGinv'){
        GG.sc <- c(GG.sc,mG+1:3)
        sc.inv <- c(sc.inv,mSC+1:3)
        GG.op <- c(GG.op,'exp','identity','negexp')
        op.inv <- c(op.inv,'log','identity','identity')
        add <- deleted[q] <- 0
      }else if(bline[q] == 'N/2'){
        GG.sc <- c(GG.sc,mG+rep(1,3))
        sc.inv <- c(sc.inv,mSC+2)
        GG.op <- c(GG.op,'sqrt2^-1','identity','sqrt2')
        op.inv <- c(op.inv,'identity','identity','identity')
        deleted[q] <- 2
      }
    }
    GG.sc <- GG.sc[-1]
    sc.inv <- sc.inv[-1]
    rm(list = c('mSC','mG','loc.1st','loc.idx','Q.loc','add'))
  }

  # count parameters, set start parameters

  numh <- numh - deleted
  numh2 <- GG.idx*3 - deleted
  numh.p <- numh[partial.idx]
  num.lam <- sum(numh2)
  num.lam.p <- sum(numh.p)
  cum.numh <- cumsum(c(0, numh))
  cum.numh.p <- cumsum(c(0, numh.p))
  names(cum.numh) <- c('0', strata)
  lam.idx <- rep(FALSE, cum.numh[Q+1])
  for(q in which(GG.idx)){
    lam.idx[(cum.numh[q]+1):cum.numh[q+1]] <- TRUE
    # lam.idx is on level of estimated parameter vector in current.status data case.
    # Also on level of cum.numh which is used for dt_1. Thats okay as non-piecewise entries are then deleted.
    # Is neither not on level of estimated or all-parameter vector if !current.status.
  }

  RE_stratum <- levels(frailty)
  J <- length(RE_stratum)
  if( (J > 1) & (length(type) == 1) ){
    type <- rep(type, J)
  }
  aNBp.idx <- type == 'aNB+p'
  aNB.idx <- type %in% c('aNB', 'NB+p')## add discrete ---
  aB.idx <- type %in% c('aB', 'B+p')## add discrete ---
  if(any(aB.idx)&!current.status) stop('aB, aB+p model only implemented for current status data so far!')
  Addams.idx <- type == 'Addams'
  ag.idx <- type == 'alpha_gamma' | Addams.idx | aNBp.idx | aNB.idx | aB.idx ## add discrete ---
  ag_idx <- frailty %in% RE_stratum[ag.idx]
  ag.type <- which(ag.idx)
  a.idx <- type %in% c('alpha','P+p')## add discrete ----
  a_idx <- frailty %in% RE_stratum[a.idx]
  a.type <- which(a.idx)
  g.idx <- type == 'gamma'
  g_idx <- frailty %in% RE_stratum[g.idx]
  g.type <- which(g.idx)
  IG.idx <- type == 'IG'
  IG_idx <- frailty %in% RE_stratum[IG.idx]
  IG.type <- which(IG.idx)
  cP.idx <- type == 'cP'
  H.idx <- type == 'Hougaard'
  PVF.idx <- type == 'PVF' | IG.idx | cP.idx | H.idx
  PVF_idx <-  frailty %in% RE_stratum[PVF.idx]
  PVF.type <- which(PVF.idx)
  no.frail <- all(type == 'nofrail')

  alpha.idx <- a.idx | (ag.idx & !aB.idx) | (PVF.idx & !IG.idx) ## add discrete ----
  gamma.idx <- g.idx | ag.idx | IG.idx | PVF.idx | aB.idx
  p.idx <- type %in% c('NB+p','P+p','B+p') ## add discrete ----
  alpha.rank <- gamma.rank <- alphagamma.rank <- rep(NA, J)
  alpha.rank[alpha.idx] <- rank(which(alpha.idx))
  gamma.rank[gamma.idx] <- rank(which(gamma.idx))
  alphagamma.rank[type == 'PVF' | ag.idx] <- rank(
    which(type == 'PVF' | ag.idx))
  num.a <- sum(alpha.idx)
  num.g <- sum(gamma.idx)
  num.p <- sum(p.idx) ## add discrete ----
  if(num.p>0 & (!current.status | cs.options$EM.cs)) stop('NB+p, P+p, B+p only implemented for current status data (non-EM optimization).')

  lnlambda <- NULL
  if(all(is.na(start_lnlam))){
    if(current.status){
      lnlambda <- rep(-10, cum.numh[Q+1])
      if(GG.special){
        for(q in which(GG.idx)){
          if(!(bline[q] %in% c('piecewise','GG','GGinv','GGpos','exponential','N/2','location'))){
            lnlambda[(cum.numh[q]+1):cum.numh[q+1]] <- c(0,log(max(y)))
          }else if(bline[q] %in%c('exponential','N/2')){
            lnlambda[(cum.numh[q]+1):cum.numh[q+1]] <- log(max(y))
          }else if(bline[q] == 'location'){
            tmp <- (cum.numh[q]+1):cum.numh[q+1]
            if(length(tmp)==3){
              lnlambda[(cum.numh[q]+1):cum.numh[q+1]] <- c(0, log(max(y)), 1)
            }else{
              lnlambda[(cum.numh[q]+1):cum.numh[q+1]] <- log(max(y))
            }
          }else if(bline[q] %in% c('GG','GGinv','GGpos')){
            lnlambda[(cum.numh[q]+1):cum.numh[q+1]] <- c(0, log(max(y)), 1)
          }
        }
      }else{
        lnlambda[lam.idx] <- rep(c(0, log(max(y)), 1), Q.GG)
      }
    }else{
      lnlambda <- rep( c(0, log(max(y)), 1), Q.GG)
      if(GG.special){
        lnlambda <- lnlambda[sc.inv]
      }
    }
  }else{
    lnlambda <- start_lnlam
  }

  if(all(is.na(start_dist))){
    alpha <- rep(-.2, num.a)
    lngamma <- rep(-2, num.g)
    lngamma[PVF.idx[gamma.idx]] <- 2
    lnp <- rep(0, num.p) ## add discrete ----
    p <- NULL
    theta <- c(alpha,lngamma,lnp)
    start_dist <- theta
  }else{
    p <- NULL # point mass only implemented for cs data, p should be implemented then after optimization. Here it is necassary for right censored data when parameter is necassary as parameter processing after optimizatzion has p in its arguments
    theta <- start_dist
  }

  slope.par <- 0
  pi <- NULL
  v.Domain <- NULL
  if(slope){
    if(is.na(start_slope)){
      if(slope.model == 'free-Willy'){
        lnphipsi <- rep(-.5, N)
        slope.par <- N+1
      }else{
        lnphipsi <- rep(-.5, 2)
        slope.par <- 3
      }
    }else{
      lnphipsi <- start_slope
      slope.par <- length(start_slope)+1
    }
    if(is.na(start_pi)){
      lnpi <- log(-log(0.5))
    }else{
      lnpi <- start_pi
    }
    X.tilde_C <- tapply(X = X.tilde, INDEX = ID, FUN = function(x) c(x, rep(0, max.n_i-length(x))), simplify = TRUE)
    X.tilde_C <- matrix(unlist(X.tilde_C), ncol = max.n_i, byrow = TRUE)
    rownames(X.tilde_C) <- unique(ID)
  }else{
    lnphipsi <- NULL
    lnpi <- NULL
    X.tilde_C <- NULL
  }

  if(all(is.na(start_beta))){
    beta <- rep(0, K)
  }else{
    beta <- start_beta
  }

  if(all(is.na(z))){
    z <- rep(1,n)
  }
  if(all(is.na(v))){
    v <- rep(0,n)
  }

  tau <- NULL
  tau.name <- NULL
  if(overdisp){
    tau <- matrix(NA,nrow=3,ncol=J)
    colnames(tau) <- RE_stratum
    rownames(tau) <- c('lntau','tau','omega')
    tau[1,] <- log(frail.thresh^-1/1000)
    tau.name <- paste('lntau_',RE_stratum,sep='')
  }

  # nlm options update
  if(current.status & !cs.options$EM.cs){
    tmp <- c(lnphipsi,lnpi,theta,beta,lnlambda)
  }else{
    tmp <- c(lnphipsi,lnpi,theta)
  }

  if(c.options$method == 'nlm'){
    c.options$typsize <- rep(c.options$typsize,
                             length(c(beta,lnlambda))+1 - length(c.options$typsize))
    c.options$stepmax <- ifelse(is.null(c.options$stepmax),
                                max(1000 * sqrt(sum((c(beta,lnlambda)/c.options$typsize)^2)), 1000)
                                ,c.options$stepmax)
  }

  if(m.options$method == 'nlm'){
    m.options$typsize <- rep(m.options$typsize,
                             length(c(tmp))+1 - length(m.options$typsize) + overdisp*J)
    m.options$stepmax <- ifelse(is.null(m.options$stepmax),
                                max(1000 * sqrt(sum((tmp/m.options$typsize)^2)), 1000),
                                m.options$stepmax)
  }

  # DEoptim options update
  if(c.options$method == 'DEoptim'){
    c.options$lower <- rep(c.options$lower, length(c(beta,lnlambda)))
    c.options$upper <- rep(c.options$upper, length(c(beta,lnlambda)))
  }
  if(m.options$method == 'DEoptim'){
    m.options$lower <- rep(m.options$lower, length(tmp) + overdisp*J)
    m.options$upper <- rep(m.options$upper, length(tmp) + overdisp*J)
  }
  # genoud options update
  if(c.options$method == 'genoud'){
    c.options$nvars <- length(c(beta,lnlambda))
  }
  if(m.options$method == 'genoud'){
    m.options$nvars <- length(tmp) + overdisp*J
  }
  rm(list = 'tmp')

  # Meassures/Matrices/indices needed in each optimaxization step
  dt <- lapply( X = breaks, FUN = function(x) diff(x))
  idx_1 <- mapply(FUN = function(y,str){
    sum(breaks[[str]] < y)
  },
  y = y, str = stratum)
  idx_1[y == 0] <- 1
  idx_1[GG_i] <- NA
  tdiff_1 <- mapply(FUN = function(y, idx, str){
    y-breaks[[str]][idx]
  },
  y = y, idx = idx_1, str = stratum)
  tr.idx_1 <- mapply(FUN = function(y,str){
    sum(breaks[[str]] < y)
  },
  y = tr, str = stratum)
  tr.idx_1[tr == 0] <- 1
  tr.idx_1[GG_i] <- NA
  tr.tdiff_1 <- mapply(FUN = function(y, idx, str){
    y-breaks[[str]][idx]
  },
  y = tr, idx = tr.idx_1, str = stratum)

  # Creating matrix of of interval lengths (baseline hazards). Needed for Matrix M.
  ## jth column has intervals of jth stratum in rows ... , 0 otherwise. # comment probably outdated. Created a vector instead
  dt_matrix <- unlist(dt)

  # Creating Selection Matrix (w.r.t. time and hazard stratum)
  zero_vec <- lam.idx + 0
  dt_1 <- tr.dt_1  <- m_1 <- tr.m_1 <- matrix(NA, nrow = length(lam.idx), ncol = n)

  for(i in 1:n){

    # cluster member (non-truncated)
    ## concerning whether individual is in stratum 1 and was still alive.
    ## Last non-zero entry is relative to survival time within the corresponding time interval
    mm_1 <- dt_i1 <- zero_vec
    left <- cum.numh[as.numeric(stratum[i])] + 1
    if(GG_i[i]){
    }else{
      idx_temp <- left:(left+idx_1[i]-1)
      dt_i1[idx_temp] <- dt_matrix[idx_temp]
      dt_i1[max(idx_temp)] <- tdiff_1[i]
      mm_1[max(idx_temp)] <- 1*(y[i]>0)
    }

    dt_1[,i] <- dt_i1
    m_1[,i] <- mm_1

    # cluster member (truncated)
    dt_i1 <- mm_1 <- zero_vec
    if(GG_i[i]){
    }else{
      idx_temp <- left:(left+tr.idx_1[i]-1)
      dt_i1[idx_temp] <- dt_matrix[idx_temp]

      dt_i1[max(idx_temp)] <- tr.tdiff_1[i]
      mm_1[max(idx_temp)] <- 1*(tr[i]>0) # neu ###
      tr.m_1[,i] <- mm_1 # vorher alleine
    }
    tr.dt_1[,i] <- dt_i1
  }

  dt_1 <- apply(X = dt_1[,partial_i,drop=FALSE], MARGIN = 2, FUN = function(dt) dt[!lam.idx])
  m_1 <- apply(X = m_1[,partial_i,drop=FALSE], MARGIN = 2,  FUN = function(m) m[!lam.idx])
  tr.dt_1 <- apply(X = tr.dt_1[,partial_i,drop=FALSE], MARGIN = 2,  FUN = function(dt) dt[!lam.idx])
  tr.m_1 <- apply(X = tr.m_1[,partial_i,drop=FALSE], MARGIN = 2,  FUN = function(m) m[!lam.idx])
  if(is.vector(dt_1)){
    dt_1 <- t(dt_1)
    m_1 <- t(m_1)
    tr.dt_1 <- t(tr.dt_1)
    tr.m_1 <- t(tr.m_1)
  }

  # Adaptions for breslow hazard. Would be more appropriate to set this in breaks_piece
  ## However, due to path dependence adaptions are made here.
  ## "dt, dt_1 etc"-part would have to be changed if changes are made in breaks_piece
  if(all(breslow.idx)){
    maxt <- stats::aggregate(y[d==1], list(stratum[d==1]), max)
    maxt[,1] <- as.numeric(maxt[,1])
    for(q in maxt[,1]){
      dt[[q]] <- c(dt[[q]][-length(dt[[q]])], maxt[q,2] - breaks[[q]][length(breaks[[q]])-1])
      breaks[[q]] <- c(breaks[[q]][-length(breaks[[q]])], maxt[q,2])
    }
    # maximum ?ndern wenn gr??te Zeit Zensierung selbe bei dt und letztem Zeitabstand
    for(i in which(breslow_i)){
      if(tdiff_1[i] < dt[[stratum[i]]][idx_1[i]]) dt_1[m_1[,i]==1,i] <- 0
      if(tr.tdiff_1[i] < dt[[stratum[i]]][tr.idx_1[i]]) tr.dt_1[tr.m_1[,i]==1,i] <- 0
    }
    dt_1 <- dt_1>0 + 0
    tr.dt_1 <- tr.dt_1>0 + 0

    rm(list = 'maxt')
  }

  breaks <- breaks[partial.idx]

  partial.cond <- all( partial.idx )
  no.piece <- all(breslow.idx) & !(all(X==0)&K==1)
  partial <- GG <- NULL

  rm(list = c('tr.m_1', 'dt', 'tdiff_1', 'tr.tdiff_1', 'idx_1', 'tr.idx_1','dt_matrix', 'left', 'dt_i1',
              'mm_1', 'breslow.idx', 'breslow.rank', 'breslow_i'))

  # EM-algorithm
  d.marg <- tapply(X = d, INDEX = ID, FUN = function(x) c(x, rep(0, max.n_i-length(x))), simplify = TRUE)
  d.marg <- matrix(unlist(d.marg), ncol = max.n_i, byrow = TRUE)
  sum.d_i <- rowSums(d.marg)
  tr.idx <- tapply(X = tr, INDEX = ID, FUN = function(x) any(x>0), simplify = TRUE)
  tr.idx <- matrix(tr.idx)

  if(current.status & !no.frail){
    # sign info matrix for sum of Probabilities in Likelihood.
    sign.info <- matrix( 1, nrow = n_C, ncol = 1+sum( choose(max(sum.d_i), 1:max(sum.d_i)) ) )
    for(i in 1:n_C){
      d_i <- sum.d_i[i]
      if(d_i>0){
        for(j in 1:d_i){
          if(j%%2!=0){
            sign.info[i,(2+sum(choose(d_i, 1:(j-1))*(j>1)) ) : (1+sum(choose(d_i, 1:j)))] <- -1
          }
        }
      }
    }

    # Versuch die SChleife in cs data zu ersetzen:
    # sign.info <- matrix( 1, nrow = n_C, ncol = 1+sum( choose(max(sum.d_i), 1:max(sum.d_i)) ) )
    # C <- matrix(0, nrow = n, ncol = n_C)
    # D <- matrix(NA, nrow = 0, ncol = n)
    # for(i in 1:n_C){
    #   col.idx <- sum(c(1,n_i)[1:i]):sum(n_i[1:i])
    #   C[col.idx,i] <- 1
    #   d_i <- sum.d_i[i]
    #   Ri <- matrix(0, nrow = d_i, ncol = n)
    #   if(d_i>0){
    #     row.idx <- sum(c(1,sum.d_i)[1:i]):sum(sum.d_i[1:i])
    #     Ri[row.idx,col.idx] <- diag(d.marg[i,])[d.marg[i,]==1,]
    #     D <- rbind(D,Ri)
    #     for(j in 1:d_i){
    #       haz.idx <- as.matrix(combn(x = d_i, m = j))
    #       sum(sum.d_i[1:i]) + haz.idx  # muesste Spaltenindex geben, Zeilen durch j
    #       if(j%%2!=0){
    #         sign.info[i,(2+sum(choose(d_i, 1:(j-1))*(j>1)) ) : (1+sum(choose(d_i, 1:j)))] <- -1
    #       }
    #     }
    #   }
    # }


  }else{sign.info <- NULL}

  # optimization:
  ## Select type of data and model and send it to specific loglihood:
  if(no.frail){ # univariate data
    opt <- loglihood.uni(X = X, iterlim = iterlim, tr = tr, y = y, d = d, stratum = stratum, stratum.GG = GG.stratum,
                         beta = beta, dt_1 = dt_1, tr.dt_1 = tr.dt_1, m_1 = m_1, GG_i = GG_i,
                         cum.numh.p = cum.numh.p, numh.p = numh.p, num.lam.p = num.lam.p, num.lam = num.lam,
                         breaks = breaks, p.strata = p.strata, Q.partial = Q.partial, Q.GG = Q.GG,
                         K = K, no.piece = no.piece, partial.cond = partial.cond, n.p = n.p, n.GG = n.GG,
                         lnlambda = lnlambda, GG.op = GG.op, c.options = c.options, R.options = R.options,
                         partial.rank=partial.rank, cum.numh=cum.numh, n_xx = c.weights,
                         sc.inv=sc.inv, current.status=current.status,lam.idx=lam.idx,
                         GG.sc = GG.sc, GG.special = GG.special, partial.idx = partial.idx, partial_i = partial_i,
                         scale.cond = scale.cond, GG.strata = GG.strata)

  }else if(current.status){ # multivariate current status data
    opt <- GG.cs(beta=beta,lnlambda=lnlambda,theta=theta,K=K,num.lam=num.lam,y=y[GG_i],tau=tau,
                 cum.numh=cum.numh,strata=strata,breaks=breaks,X=X,d=d,dt_1=dt_1,numh.p=numh.p,
                 z=z,n=n,n_C=n_C,n.GG=n.GG,ag_i=ag_idx,a_i=a_idx,g_i=g_idx,PVF_i=PVF_idx,
                 ag.idx=ag.idx,a.idx=a.idx,g.idx=g.idx,PVF.idx=PVF.idx,IG.idx=IG.idx,GG.idx=GG.idx,
                 cs.options=cs.options,frail.thresh=frail.thresh,initial.z=initial.z,start_dist=start_dist,
                 num.a=num.a,num.g=num.g,J=J,alpha.idx=alpha.idx,gamma.idx=gamma.idx, Addams.idx = Addams.idx,
                 aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,aB.idx=aB.idx,N=N,cP.idx=cP.idx,H.idx=H.idx,alpha.rank=alpha.rank,gamma.rank=gamma.rank,
                 idx=as.numeric(frailty),sum.d_i=sum.d_i,max.n_i=max.n_i,ID=ID,d.marg=d.marg,overdisp=overdisp,
                 prnt=prnt,iterlim=iterlim,partial.idx=partial.idx,c.n_xx=c.weights,n_xx=m.weights,n.date=n.date,
                 partial_i=partial_i,partial.rank=partial.rank,Q=Q.GG,stratum=GG.stratum,sign.info=sign.info,
                 lam.idx=lam.idx,GG_i=GG_i,num.lam.p=num.lam.p,GG.op=GG.op,GG.sc=GG.sc,GG.special=GG.special,
                 sc.inv=sc.inv,scale.cond=scale.cond,scale.marginal=scale.marginal,converge=converge,
                 c.options = c.options, m.options = m.options, R.options = R.options,
                 num.p = num.p, p.idx = p.idx## add discrete ----
    )
  }else{ # right censored, clustered data

    if(!multi.num){ # bivariate data
      if(partial.cond){ # only non-parametric or piecewise constant hazard
        if(no.piece){ # no-truncation model
          scale.cond <- 1
          opt <- partial.bi(X = X, stratum = stratum, beta = beta, theta = theta, y = y, tr = tr, d = d,
                            dt_1 = dt_1, tr.dt_1 = tr.dt_1, m_1 = m_1, cum.numh.p = cum.numh.p, K = K,
                            numh.p = numh.p, num.lam.p = num.lam.p, breaks = breaks, p.strata = p.strata,
                            Q.partial = Q.partial, n.p = n.p,partial.stratum = partial.stratum,
                            num.a = num.a, num.g = num.g,
                            d.marg = d.marg, d_1 = d.marg[,1], d_2 = d.marg[,2], tr.idx = tr.idx,
                            n_C = n_C, J = J, gamma.rank = gamma.rank, alpha.rank = alpha.rank, Addams.idx = Addams.idx,
                            aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                            ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx, IG.idx = IG.idx, PVF.idx = PVF.idx,
                            ag_idx = ag_idx, a_idx = a_idx, g_idx = g_idx, IG_idx = IG_idx, PVF_idx = PVF_idx,
                            idx = as.numeric(frailty), frailty = frailty, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                            ag.type = ag.type, a.type = a.type, g.type = g.type, z = z,
                            IG.type = IG.type, PVF.type = PVF.type,  ID = ID, iterlim = iterlim,
                            scale.marginal = scale.marginal,frail.thresh=frail.thresh,initial.z=initial.z,start_dist=start_dist,
                            converge = converge, prnt = prnt,
                            m.options = m.options, R.options = R.options)
        }else{ # profile loglihood including truncation
          opt <- profile.bi(X = X, stratum = stratum, beta = beta, theta = theta, tr = tr, y = y, d = d,
                            dt_1 = dt_1, tr.dt_1 = tr.dt_1, m_1 = m_1, cum.numh.p = cum.numh.p, K = K,
                            numh.p = numh.p, num.lam.p = num.lam.p, breaks = breaks, p.strata = p.strata,
                            Q.partial = Q.partial, n.p = n.p, partial.stratum = partial.stratum,
                            num.a = num.a, num.g = num.g, initial.z=initial.z, start_dist=start_dist,
                            d.marg = d.marg, d_1 = d.marg[,1], d_2 = d.marg[,2], tr.idx = tr.idx,
                            n_C = n_C, J = J, gamma.rank = gamma.rank, alpha.rank = alpha.rank,
                            ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx, IG.idx = IG.idx, PVF.idx = PVF.idx, Addams.idx = Addams.idx,
                            aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                            ag_idx = ag_idx, a_idx = a_idx, g_idx = g_idx, IG_idx = IG_idx, PVF_idx = PVF_idx,
                            idx = as.numeric(frailty), frailty = frailty, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                            ag.type = ag.type, a.type = a.type, g.type = g.type, ID = ID, z = z,
                            IG.type = IG.type, PVF.type = PVF.type,  scale.cond = scale.cond, prnt = prnt, iterlim = iterlim,
                            scale.marginal = scale.marginal, converge = converge,
                            frail.thresh=frail.thresh, c.options=c.options, m.options=m.options, R.options=R.options)
        }
      }else{ # GG-baseline
        opt <- GG.bi(theta = theta, X = X, partial_i = partial_i, GG_i =GG_i, beta = beta,y=y,
                     lnlambda = lnlambda, K = K, num.lam = num.lam, cum.numh.p = cum.numh.p,n = n,
                     numh.p = numh.p, num.lam.p = num.lam.p, breaks =breaks, p.strata = p.strata,
                     tr = tr, n.p = n.p, n.GG = n.GG, d = d, Q.p = Q.partial, Q.GG = Q.GG, stratum.GG = GG.stratum,
                     dt_1 = dt_1, tr.dt_1 = tr.dt_1, m_1 = m_1,
                     ID = ID, scale.cond = scale.cond, scale.marginal = scale.marginal,
                     d_1 = d.marg[,1], d_2 = d.marg[,2], d.marg = d.marg, GG.strata = GG.strata,
                     partial.idx = partial.idx, num.a = num.a, num.g = num.g, J = J, n_C = n_C, tr.idx = tr.idx,
                     gamma.rank = gamma.rank, alpha.rank = alpha.rank,
                     ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx, IG.idx = IG.idx, PVF.idx = PVF.idx, Addams.idx = Addams.idx,
                     ag_idx = ag_idx, a_idx = a_idx, g_idx = g_idx, IG_idx = IG_idx, PVF_idx = PVF_idx,
                     aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                     idx = as.numeric(frailty), frailty = frailty, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                     ag.type = ag.type, a.type = a.type, g.type = g.type, z = z,
                     GG.op = GG.op, GG.sc = GG.sc, GG.special = GG.special,initial.z=initial.z,start_dist=start_dist,
                     IG.type = IG.type, PVF.type = PVF.type, iterlim = iterlim,frail.thresh=frail.thresh,
                     converge = converge, prnt = prnt,
                     c.options = c.options, m.options = m.options, R.options = R.options)
      }

    }else{ # Multivariate case (numeric optimaxization only)
      #fac <- NULL
      #ag.L <- Laplace_derivatives_AG#data(Laplace_derivatives_AG)
      if(any(sum.d_i > 14 & type[frailty] %in% c('Addams', 'alpha_gamma'))){
        #fac <-  data(file = 'fac')
        if(any(sum.d_i > 400)){
          print('New factors for derivative of alpha-gamma Laplace are created as
                max(#events per cluster) is higher than 400. Might take forever and probably will not work anyway!!!')
          fac <- add.combos(order = max(sum.d_i+1), combos = fac)
        }
      }
      if(partial.cond){ # only non-parametric or piecewise constant hazard
        if(no.piece){ # no-truncation model
          scale.cond <- 1
          opt <- partial.multi(X = X, stratum = stratum, beta = beta, theta = theta, y = y, tr = tr, d = d,
                               dt_1 = dt_1, tr.dt_1 = tr.dt_1, m_1 = m_1, cum.numh.p = cum.numh.p, K = K,
                               numh.p = numh.p, num.lam.p = num.lam.p, breaks = breaks, p.strata = p.strata,
                               Q.partial = Q.partial, n.p = n.p, slope.model = slope.model,
                               partial.stratum = partial.stratum,N=N,X.tilde=X.tilde,X.tilde_C=X.tilde_C,slope=slope,lnphipsi=lnphipsi,lnpi=lnpi,
                               num.a = num.a, num.g = num.g, d.marg = d.marg, sum.d_i = sum.d_i, tr.idx = tr.idx,
                               n_C = n_C, J = J, gamma.rank = gamma.rank, alpha.rank = alpha.rank, Addams.idx = Addams.idx,
                               ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx, IG.idx = IG.idx, PVF.idx = PVF.idx,
                               ag_idx = ag_idx, a_idx = a_idx, g_idx = g_idx, IG_idx = IG_idx, PVF_idx = PVF_idx,
                               aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                               idx = as.numeric(frailty), frailty = frailty, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                               ag.type = ag.type, a.type = a.type, g.type = g.type,
                               IG.type = IG.type, PVF.type = PVF.type,  ID = ID, iterlim = iterlim,frail.thresh=frail.thresh,
                               scale.marginal = scale.marginal, max.n_i = max.n_i, z = z, v = v, fac = fac,initial.z=initial.z,start_dist=start_dist,
                               converge = converge, prnt = prnt,m.options = m.options, R.options = R.options)
        }else{ # profile loglihood including truncation
          opt <- profile.multi(X = X, stratum = stratum, beta = beta, theta = theta, tr = tr, y = y, d = d,
                               dt_1 = dt_1, tr.dt_1 = tr.dt_1, m_1 = m_1, cum.numh.p = cum.numh.p, K = K,
                               numh.p = numh.p, num.lam.p = num.lam.p, breaks = breaks, p.strata = p.strata,
                               Q.partial = Q.partial, n.p = n.p,v=v, X.tilde=X.tilde,X.tilde_C=X.tilde_C,slope=slope,
                               partial.stratum = partial.stratum,slope.model=slope.model,lnphipsi=lnphipsi,lnpi=lnpi,N=N,
                               num.a = num.a, num.g = num.g, d.marg = d.marg, sum.d_i = sum.d_i, tr.idx = tr.idx,
                               n_C = n_C, J = J, gamma.rank = gamma.rank, alpha.rank = alpha.rank, Addams.idx = Addams.idx,
                               ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx, IG.idx = IG.idx, PVF.idx = PVF.idx,
                               ag_idx = ag_idx, a_idx = a_idx, g_idx = g_idx, IG_idx = IG_idx, PVF_idx = PVF_idx,
                               idx = as.numeric(frailty), frailty = frailty, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                               ag.type = ag.type, a.type = a.type, g.type = g.type, ID = ID, max.n_i = max.n_i,
                               aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx, fac = fac,initial.z=initial.z,start_dist=start_dist,
                               IG.type = IG.type, PVF.type = PVF.type,  scale.cond = scale.cond, z = z, prnt = prnt, iterlim = iterlim,
                               scale.marginal = scale.marginal, converge = converge,frail.thresh=frail.thresh,
                               c.options = c.options, m.options = m.options, R.options = R.options)
        }
      }else{ # GG-baseline
        opt <- GG.multi(theta = theta, X = X, partial_i = partial_i, GG_i =GG_i,beta = beta,
                        lnlambda = lnlambda, K = K, num.lam = num.lam, cum.numh.p = cum.numh.p,n = n,
                        numh.p = numh.p, num.lam.p = num.lam.p, breaks =breaks, p.strata = p.strata,
                        tr=tr,Q.p=Q.partial,Q.GG=Q.GG,y=y,d=d,n.GG=n.GG,n.p=n.p,stratum.GG=GG.stratum,
                        dt_1=dt_1,tr.dt_1=tr.dt_1,m_1=m_1,X.tilde=X.tilde,X.tilde_C=X.tilde_C,lnpi=lnpi, lnphipsi=lnphipsi,
                        ID = ID, scale.cond = scale.cond, scale.marginal = scale.marginal,slope=slope,slope.model=slope.model,
                        sum.d_i = sum.d_i, d.marg = d.marg, GG.strata = GG.strata, max.n_i = max.n_i,v=v,N=N,
                        partial.idx = partial.idx, num.a = num.a, num.g = num.g, J = J, n_C = n_C, tr.idx = tr.idx,
                        gamma.rank = gamma.rank, alpha.rank = alpha.rank, Addams.idx = Addams.idx,
                        ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx, IG.idx = IG.idx, PVF.idx = PVF.idx,
                        aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                        ag_idx = ag_idx, a_idx = a_idx, g_idx = g_idx, IG_idx = IG_idx, PVF_idx = PVF_idx,  z = z,
                        idx = as.numeric(frailty), frailty = frailty, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                        ag.type = ag.type, a.type = a.type, g.type = g.type,  fac = fac,initial.z=initial.z,start_dist=start_dist,
                        IG.type = IG.type, PVF.type = PVF.type, GG.op = GG.op, GG.sc = GG.sc, GG.special = GG.special,
                        iterlim = iterlim, converge = converge, prnt = prnt,frail.thresh=frail.thresh,
                        c.options = c.options, m.options = m.options, R.options = R.options)
      }
    }
  }

  list2env(x = opt, envir = environment())
  rm(list = 'opt')
  #browser()

  if(any(Addams.idx)){
    g.idx[Addams.idx] <- abs(model_par[Addams.idx,'alpha'])< 1e-4
    a.idx[Addams.idx] <- abs(model_par[Addams.idx,'alpha']-model_par[Addams.idx,'gamma'])<1e-4 & !g.idx[Addams.idx]
    ag.idx[Addams.idx] <- !a.idx & !g.idx & Addams.idx
    ag_idx <- ag.idx[frailty]
    a_idx <- a.idx[frailty]
    g_idx <- g.idx[frailty]
    type[ag.idx] <- 'alpha_gamma'
    type[g.idx] <- 'gamma'
    type[a.idx] <- 'alpha'
  }

  # Gradient and Hessian
  if(no.frail){ # S,H: no frailty model
    if(current.status){ # current status data (+no.frail)
      if(!any(GG.idx)){ # piecewise hazard (+cs+no.frail)
        S <- dcs.dblam(par = c(beta,lnlambda), X=X,d=d,K=K,n_xx=c.weights,
                       num.lam.p=num.lam.p,dt_1=dt_1,scale.cond=abs(scale.cond))
        H <- numDeriv::hessian(func = cond.p.cs, x = c(beta,lnlambda), method = 'Richardson',z=z, z2=z,
                               method.args=R.options,n_xx=c.weights,
                               X=X,d=d,K=K,num.lam.p=num.lam.p,dt_1=dt_1,scale.cond=abs(scale.cond))
      }else{ # GG hazard (+cs+no.frail)
        D <- numDeriv::genD(func = cond.GG.cs, x = c(beta,lnlambda), method = 'Richardson',
                            method.args=R.options, d = d, y = y[GG_i], n_xx = c.weights,
                            stratum = GG.stratum, n = n, dt_1 = dt_1,
                            n.GG = n.GG, Q = Q.GG, X = X, z = rep(1,n), z2=z, K = K, num.lam = num.lam, num.lam.p =num.lam.p,
                            partial.idx = partial.idx, partial_i = partial_i, GG_i = GG_i, lam.idx = lam.idx,
                            GG.sc = GG.sc, GG.op = GG.op, GG.special = GG.special,
                            scale.cond = abs(scale.cond))$D
        np <- K+num.lam.p+num.lam
        S <- D[1:np]
        D <- D[-c(1:np)]
        H <- matrix(NA, nrow = np, ncol = np)
        H[upper.tri(H,diag = TRUE)] <- D
        tmp <- which(is.na(H), arr.ind = TRUE)
        tmp2 <- cbind(tmp[,2],tmp[,1])
        H[tmp] <- H[tmp2]
      }
    }else{ # right censored data
      if(partial.cond & no.piece){ # coxph
        S <- dpartial.db(beta = beta, z = rep(1,n), X = X, K = K, d = d, y = y,
                         dt = dt_1, tr.dt = tr.dt_1, m = m_1, n = n, scale.cond = abs(scale.cond))
        H <- solve(-cond$var)
      }else if(partial.cond & !no.piece){ # piecewise hazard or truncation (+rc)
        S <- dpartial.db(beta = beta, z = rep(1,n), X = X, K = K, d = d, y = y,
                         dt = dt_1, tr.dt = tr.dt_1, m = m_1, n = n, scale.cond = abs(scale.cond))

        H <- d2partial.db2(beta = beta, z = rep(1,n), X = X, K = K, d = d, y = y, tr = tr,
                           dt = dt_1, tr.dt = tr.dt_1, m = m_1, n = n, scale.cond = abs(scale.cond))
      }else{ # GG hazard (+ rc)
        D <- numDeriv::genD(func = cond.loglihood, x = c(beta,lnlambda), method = 'Richardson',
                            method.args=R.options,
                            X.GG = X[GG_i,,drop=FALSE], y.GG = y[GG_i], tr.GG = tr[GG_i], z.GG = rep(1,n.GG),
                            d.GG = d[GG_i],n.GG = n.GG, X.p=X[partial_i,,drop=FALSE], z.p=rep(1,n.p), d.p=d[partial_i],
                            n.p=n.p, Q.GG = Q.GG, stratum.GG = GG.stratum, dt = dt_1, m = m_1, tr.dt = tr.dt_1,
                            num.lam=num.lam,K=K,GG.op=GG.op,GG.sc=GG.sc,GG.special=GG.special,scale.cond=abs(scale.cond))$D
        np <- K+num.lam
        S <- D[1:np]
        D <- D[-c(1:np)]
        H <- matrix(NA, nrow = np, ncol = np)
        H[upper.tri(H,diag = TRUE)] <- D
        tmp <- which(is.na(H), arr.ind = TRUE)
        tmp2 <- cbind(tmp[,2],tmp[,1])
        H[tmp] <- H[tmp2]
      }
    }

  }else{
    rownames(model_par) <- RE_stratum
    if(current.status){
      D <- numDeriv::genD(func = dcs.multi, x = c(theta,lnp,beta,lnlambda,tau['lntau',]), method = 'Richardson',
                          method.args=R.options,n_xx=m.weights,sign.info=sign.info,
                          K=K,num.lam=num.lam,y=y[GG_i],X=X,dt_1=dt_1,overdisp=overdisp,n.date=n.date,
                          n=n,n_C=n_C,n.GG=n.GG,ag_i=ag_idx,a_i=a_idx,g_i=g_idx,PVF_i=PVF_idx,
                          ag.idx=ag.idx,a.idx=a.idx,g.idx=g.idx,PVF.idx=PVF.idx,IG.idx=IG.idx,
                          num.a=num.a,num.g=num.g,J=J,alpha.idx=alpha.idx,gamma.idx=gamma.idx,
                          aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,aB.idx=aB.idx,N=N,cP.idx=cP.idx,H.idx=H.idx,
                          idx=as.numeric(frailty),sum.d_i=sum.d_i,max.n_i=max.n_i,ID=ID,d.marg=d.marg,
                          partial.idx=partial.idx,GG.idx=GG.idx,
                          partial_i=partial_i,Q=Q.GG,stratum=GG.stratum,frail.thresh=frail.thresh,
                          lam.idx=lam.idx,GG_i=GG_i,num.lam.p=num.lam.p,GG.op=GG.op,GG.sc=GG.sc,GG.special=GG.special,
                          sc.inv=sc.inv,scale.marginal=abs(scale.marginal),num.p=num.p,p.idx=p.idx)$D ### add discrete ----
      np <- num.a+num.g+num.p+K+num.lam.p+num.lam+overdisp*J
    }else{ # right censored multivariate data
      D <- numDeriv::genD(func = dmarginal.multi, x = c(theta,beta,lnlambda,lnpi,lnphipsi), method = 'Richardson',
                          method.args=R.options,
                          X = X, partial_i = partial_i, GG_i =GG_i,GG.idx=GG.idx,
                          K = K, num.lam = num.lam, cum.numh.p = cum.numh.p,n = n, fac = fac,
                          numh.p = numh.p, num.lam.p = num.lam.p, breaks =breaks, p.strata = p.strata,
                          tr=tr,Q.p=Q.partial,Q.GG=Q.GG,y=y,d=d,n.GG=n.GG,n.p=n.p,stratum.GG=GG.stratum,
                          dt_1=dt_1,tr.dt_1=tr.dt_1,m_1=m_1,X.tilde=X.tilde,X.tilde_C=X.tilde_C,
                          ID = ID, scale.marginal = abs(scale.marginal),slope=slope,slope.model=slope.model,
                          sum.d_i = sum.d_i, d.marg = d.marg, GG.strata = GG.strata, max.n_i = max.n_i,v=v,N=N,
                          partial.idx = partial.idx, num.a = num.a, num.g = num.g, J = J, n_C = n_C, tr.idx = tr.idx,
                          gamma.rank = gamma.rank, alpha.rank = alpha.rank,
                          ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx, IG.idx = IG.idx, PVF.idx = PVF.idx,
                          aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                          ag_idx = ag_idx, a_idx = a_idx, g_idx = g_idx, IG_idx = IG_idx, PVF_idx = PVF_idx,  z = z,
                          idx = as.numeric(frailty), frailty = frailty, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                          ag.type = ag.type, a.type = a.type, g.type = g.type, frail.thresh=frail.thresh,
                          IG.type = IG.type, PVF.type = PVF.type, GG.op = GG.op, GG.sc = GG.sc, GG.special = GG.special)$D
      np <- num.a+num.g+K+num.lam+slope.par
    }
    S <- D[1:np]
    D <- D[-c(1:np)]
    H <- matrix(NA, nrow = np, ncol = np)
    H[upper.tri(H,diag = TRUE)] <- D
    tmp <- which(is.na(H), arr.ind = TRUE)
    tmp2 <- cbind(tmp[,2],tmp[,1])
    H[tmp] <- H[tmp2]
  }

  if(all(X==0)&K==1){
    S <- S[-(num.a+num.g+num.p+slope.par+K)]### add discrete ----
    H <- H[,-(num.a+num.g+num.p+slope.par+K),drop=FALSE]### add discrete ----
    H <- H[-(num.a+num.g+num.p+slope.par+K),,drop=FALSE]### add discrete ----
    K <- 0
  }

  Sigma2 <- tryCatch(expr = solve(-H),
                     error = function(e) {warning('Hessian is singular. Handling can be specified by H.singular!');
                       diag(Inf, nrow = dim(H)[1])} )
  H.check <- all(is.infinite(diag(Sigma2)))
  H.check2 <- diag(Sigma2<0)
  if((H.check|any(H.check2,na.rm=TRUE)) & H.singular != 'standard'){
    if(any(H.check2)){
      tmp <- which(H.check2)
      warning(paste(c('Sigma^2 had negative Variances. This was the case for parmeters',
                      paste(tmp,collapse=','),'according to order in H.'),collapse=' '))
    }else if(H.check){
      tmp <- which(abs(diag(H)) < 1e-7)
      warning(paste(c('The following diagonal elements in H are very close to zero:',
                      paste(tmp,collapse=','),'. Might be the reason for H being singular.',
                      'Often happens if parameter is close to boundary!'),collapse=' '))
    }

    if(H.singular == 'gcholesky'){
      Sigma2 <- solve(bdsmatrix::gchol(-H))
    }else if(H.singular == 'ginverse'){
      Sigma2 <- MASS::ginv(-H)
    }else{
      stop(paste('H.singular can take the values standard, gcholesky and ginverse.',
                 'Your incorrect choice was',H.singular))
    }
  }

  # del parameter of interest / del underlying parameter
  ## (underlying parameter=log(parameter of interest) for example):
  del <- model_par
  del[PVF.idx & !cP.idx & !H.idx,'alpha'] <- del[PVF.idx & !cP.idx & !H.idx,'alpha']  + 1
  #del[aNB.idx,'alpha'] <- del[aNB.idx,'alpha'] * log(del[aNB.idx,'alpha']/del[aNB.idx,'gamma']+0) # caution DR matrix not a diagonal matrix anymore!
  del[aNB.idx,'gamma'] <- del[aNB.idx,'gamma'] - del[aNB.idx,'alpha'] # caution DR matrix not a diagonal matrix anymore!, neue Parameterisierung
  del[H.idx,'alpha'] <- del[H.idx,'alpha']*log(del[H.idx,'alpha']+0) # negative sign cancels out
  del[ag.idx & !aNBp.idx & !aNB.idx,'alpha'] <- 1 + del[ag.idx & !aNBp.idx & !aNB.idx,'alpha']*0 # last summand is important for absence of frailty
  if(current.status){
    del.lam <- lambda
    # current.status:
    ## lambda includes all estimaterd parameters in transormed format, including piecewise,
    ## not in GG format but special case (i.e. no redunandant special case parameters)
    if(GG.special){
      tmp <- lam.idx
      tmp[tmp] <- GG.op[sc.inv] == 'identity'
      del.lam[tmp] <- 1
    }
    if(any(GG.idx) & !GG.special){
      del.lam[c(cum.numh[c(GG.idx,F)]+2,cum.numh[c(GG.idx,F)]+3)] <- 1
    }
  }else{
    del.lam <- lambda*rep(c(1,0,0),Q.GG) + rep(c(0,1,1),Q.GG)
    if(GG.special & any(GG.idx)){
      del.lam <- del.lam[sc.inv]
      # right censored
      ## lambda includes redundant GG parameters
    }
  }

  if(slope){
    del.B <- c(pi*log(pi), exp(lnphipsi))
  }else{
    del.B <- NULL
  }

  Delta.rule <- diag(
    c(del[alpha.idx,'alpha',drop=FALSE], del[gamma.idx, 'gamma',drop=FALSE], p[p.idx],
      del.B, sd_x^-1, del.lam, tau['lntau',]*0+1), nrow=dim(H)[1], ncol = dim(H)[2])### add discrete ---- !!! Fehler, findet p nicht bei Twin gamma Breslow model, no strat, no p!!!
  if(any(aNB.idx)){
    # tmp <- sum(alpha.idx <= which(aNB.idx))
    # tmp2 <- sum(gamma.idx <= which(aNB.idx))
    # Delta.rule[tmp, num.a+tmp2] <- del[aNB.idx,'alpha']

    tmp <- t(replicate(which(alpha.idx), n = sum(aNB.idx)))
    NB.idx <- replicate(which(aNB.idx), n = sum(alpha.idx))
    tmp <- rowSums(tmp <= NB.idx)
    tmp2 <- t(replicate(which(gamma.idx), n = sum(aNB.idx)))
    NB.idx <- replicate(which(aNB.idx), n = sum(gamma.idx))
    tmp2 <- rowSums(tmp2 <=NB.idx) # neue Parameterisierung
    for(j in 1:length(tmp)){
      Delta.rule[num.a+tmp2[j],tmp[j]] <- del[aNB.idx,'alpha'][j]# neue PArameterisierung
    }
  }
  if(slope & slope.model == 'free-Willy' & N>2){
    for(i in 1:(N+1)){
      if(i < N/2){
        Delta.rule[num.a+num.g+i+1,(num.a+num.g+i+2):(num.a+num.g+1+N/2)] <- -exp(lnphipsi[(i+1):(N/2)]) # der erste wird nie gebraucht
      }
      if(i>N/2+2){
        Delta.rule[num.a+num.g+i,(num.a+num.g+1+N/2+1):(num.a+num.g+i-1)] <- exp(lnphipsi[(N/2+1):(i-2)])# der letzte wird nie gebraucht
      }
    }
  }
  Sigma2.DR <- Delta.rule%*%Sigma2%*%t(Delta.rule)
  sd <- diag(as.matrix(Sigma2))
  sd.DR <- diag(as.matrix(Sigma2.DR))
  if(any(sd<0) | any(is.nan(sd))){
    tmp <- which(sd<0 | is.nan(sd))
    warning(paste(
      'Diagonal element(s) ', paste(tmp, collapse = ','),
      ' of Sigma^2 are < 0. Artificially set to zero.'),sep = ' ')
    sd[tmp] <- 0
    sd.DR[tmp] <- 0
  }
  sd <- sqrt(sd)
  sd.DR <- sqrt(sd.DR)

  if(K>0){
    beta_OUT <- matrix(NA, nrow = K, ncol = 5)
    colnames(beta_OUT) <- c('Estimate', 'SE', 'CI-', 'CI+', 'p-value')
    rownames(beta_OUT) <- colnames(X)
    if(any(sd_x!=1)){
      warning('baseline hazard parameters are based on scaled (and probably centered) data! See centerX, scaleX option.
              plot.predict element contains beta,X on the scale of the data during estimation!')
    }
    beta_OUT[,1] <- beta/sd_x
    beta_OUT[,2] <- sd.DR[(num.a+num.g+num.p+1):(num.a+num.g+num.p+K)] ### add discrete ----
    beta_OUT[,3:5] <-  CI.sig(par = beta_OUT[,1], Sigma = beta_OUT[,2], level = level)# 23.02.24: vorher war hier beta statt beta_OUT[,2]

    HR_OUT <- beta_OUT
    HR_OUT[,c(1,3,4)] <- exp(HR_OUT[,c(1,3,4)])
    HR_OUT[,2] <- exp(HR_OUT[,1])*HR_OUT[,2]
  }else{
    beta_OUT <- HR_OUT <- 'no-covariate model'
  }

  Hazard <- matrix(NA, nrow = num.lam + num.lam.p, ncol = 4)
  colnames(Hazard) <- c('Estimate', 'SE', 'CI-', 'CI+')
  tmp.name <- tmp.name2 <- NULL
  if(GG.special){
    if(current.status){
      tmp <- lnlambda[lam.idx]
      tmp.se <- sd[(num.a+num.g+num.p+K) + 1:(num.lam+num.lam.p)][lam.idx]### add discrete ----
    }else{
      tmp <- lnlambda
      tmp.se <- sd[(num.a+num.g+num.p+K) + 1:num.lam]
    }
    CI <- CI.sig(par = tmp, Sigma = tmp.se,
                 level = level, e = (GG.op[sc.inv] == 'exp'), p = FALSE)
    #lnlambda enth?lt alle bei current.status, aber keine Doppelten
    # current status piecwise hazards raus nehmen, werden eh unten gerehcnet!
    # Problem wird lnlambda, da ist current.status piecewise hazard drinne. [lam.idx] nutzen
    # auch bei sd
    # oben piecewise bei current status raus nehmen!
  }
  for(q in 1:Q){
    if(GG.idx[q]){
      idx <- GG.rank[q]
      tmp <- c(sigma[idx], eta[idx], nu[idx])
      tmp.tmpname <- paste(strata[q], c('sigma', 'eta', 'nu'), sep = ':')
      tmp.name2 <-  c(tmp.name2, tmp.tmpname)
      if(current.status){
        idx <- (cum.numh[q]+1):(cum.numh[q+1]) # necassry because current.status includes piecewise hazard pars
      }else{
        idx <- ((idx-1)*3+1):(idx*3)
      }
      if(GG.special){
        if(q>1 & !current.status) idx <- idx - sum(deleted[1:(q-1)])
        if(deleted[q]>0){
          if(!current.status){
            idx <- idx[-3*(deleted[q]==1) + 1*(deleted[q]==2)]
          }
          tmp <- tmp[-3*(deleted[q]==1) + 2*(deleted[q]==2)]
          tmp.tmpname <- tmp.tmpname[-3*(deleted[q]==1) + 2*(deleted[q]==2)]
        }
        tmp.se <- sd.DR[idx+num.a+num.g+num.p+K]### add discrete ----
        if(current.status & any(!is.na(partial.rank[1:q]))){
          idx <- idx - sum(numh.p[max(partial.rank[1:q],na.rm=TRUE)])
        }
        tmp.CI <- CI[idx,]
      }else{
        idx <- idx + num.a+num.g+num.p+K### add discrete ----
        tmp.se <- sd.DR[idx]
        tmp.CI <- CI.sig(par = c(log(tmp[1]), tmp[2:3]), Sigma = sd[idx],
                         level = level, e = c(TRUE, FALSE, FALSE), p = FALSE)
      }
      tmp.name <- c(tmp.name, tmp.tmpname)
    }else{
      idx <- partial.rank[q]
      tmp <- partial.lambda[[idx]]
      tmp.name <- c(tmp.name, paste(strata[q], names(partial.lambda[[idx]]), sep = ':'))
      if(current.status){
        idx <- (cum.numh[q]+1):(cum.numh[q+1]) + num.a+num.g+num.p+K### add discrete ----
        if(GG.special){
          if(deleted[q]>0){
            tmp <- tmp[-3*(deleted[q]==1) + 2*(deleted[q]==2)]
          }
        }
        tmp.se <- sd.DR[idx]
        tmp.CI <- CI.sig(par = log(tmp),
                         Sigma = sd[idx],
                         level = level, e = rep(TRUE, numh[q]), p = FALSE)
      }else{
        tmp.se <- NA
        tmp.CI <- NA
      }
    }
    Hazard[(1+cum.numh[q]):cum.numh[q+1],1] <- tmp
    Hazard[(1+cum.numh[q]):cum.numh[q+1],2] <- tmp.se
    Hazard[(1+cum.numh[q]):cum.numh[q+1],3:4] <- tmp.CI
  }
  if(current.status){
    tmp.name2 <- tmp.name
  }else if(GG.special){
    tmp.name2 <- tmp.name2[sc.inv]
  }
  rownames(Hazard) <- tmp.name

  if(all(!GG.idx) & !current.status){
    Hazard <- Hazard[,1,drop=FALSE]
  }

  theta_OUT <- matrix(NA, nrow = num.a+num.g+sum(IG.idx), ncol = 4)
  colnames(theta_OUT) <- c('Estimate', 'SE', 'CI-', 'CI+')
  tmp.name <- NULL
  for(j in 1:J){
    if(!is.na(alpha.rank[j]) | IG.idx[j]){
      tmp.name <- c(tmp.name, paste(RE_stratum[j], ': alpha', sep = ''))
    }
  }
  for(j in 1:J){
    if(!is.na(gamma.rank[j])){
      tmp.name <- c(tmp.name, paste(RE_stratum[j], ': gamma', sep = ''))
    }
  }
  rownames(theta_OUT) <- tmp.name
  alpha[aB.idx] <- NA

  theta_OUT[,1] <- c(alpha[!is.na(alpha)],gamma[!is.na(gamma)])
  for(j in 1:J){
    theta_OUT[alpha.rank[j],2] <- sd.DR[alpha.rank[j]]
    if(PVF.idx[j] & !IG.idx[j]){
      theta_OUT[alpha.rank[j],3:4] <- CI.sig(par = theta[alpha.rank[j]], Sigma = sd[alpha.rank[j]],
                                             level = level, e = TRUE, p = FALSE) - 1
    }
    if(IG.idx[j]){
      theta_OUT[sum(alpha.idx[1:j], na.rm = TRUE)+sum(IG.idx[1:j]),2:4] <- c(0,-0.5,-0.5)
    }
    if(ag.idx[j] & !aB.idx[j] & !aNB.idx[j]){
      theta_OUT[alpha.rank[j],3:4] <- CI.sig(par = theta[alpha.rank[j]], Sigma = sd[alpha.rank[j]], level = level,
                                             e = FALSE, p = FALSE)
    }
    if(a.idx[j] | aNB.idx[j]){ ## add discrete ----
      theta_OUT[alpha.rank[j],3:4] <- CI.sig(par = theta[alpha.rank[j]], Sigma = sd[alpha.rank[j]], level = level,
                                             e = TRUE, p = FALSE)
    }
  }
  if(num.g>0){
    theta_OUT[(num.a+sum(IG.idx)+1):(num.a+sum(IG.idx)+num.g),2] <- sd.DR[(num.a+1):(num.a+num.g)]
    theta_OUT[(num.a+sum(IG.idx)+1):(num.a+sum(IG.idx)+num.g),3:4] <- CI.sig(par = theta[(num.a+1):(num.a+num.g)],
                                                                             Sigma = sd[(num.a+1):(num.a+num.g)],
                                                                             level = level, e = TRUE, p = FALSE)
    if(any(aNB.idx)){
      warning('gamma CI for aNB, NB+p models is only for additive gamma component and ignores insecurity from alpha.
              Note, gamma = alpha + gamma*, gamma* = exp(gamma~). CI is CI for gamma*.')
    }
  }

  if(num.p==0){### add discrete ----
    name.p <- p_OUT <- NULL
  }else{
    p_OUT <- matrix(NA, nrow = num.p, ncol = 4)
    colnames(p_OUT) <- c('Estimate', 'SE', 'CI-', 'CI+')
    rownames(p_OUT) <- name.p <- paste('p:',RE_stratum[p.idx])
    p_OUT[,1] <- p[p.idx]
    p_OUT[,2] <- sd.DR[(num.a+num.g+1):(num.a+num.g+num.p)]
    p_OUT[,3:4] <- CI.sig(par = lnp,
                          Sigma = sd[(num.a+num.g+1):(num.a+num.g+num.p)],
                          level = level, e = TRUE, p = FALSE)
  }

  name.phi <- NULL
  if(slope){
    phi_OUT <- matrix(NA, nrow = N+1, ncol = 5)
    colnames(phi_OUT) <- c('Estimate', 'prob', 'SE', 'CI-', 'CI+')
    phi_OUT[,'Estimate'] <- v.Domain
    phi_OUT[N/2+1,] <- 0
    if(slope.model == 'free-Willy'){
      phi_OUT[-c(N/2+1),'SE'] <- sd.DR[(num.a+num.g+num.p+2):(num.a+num.g+num.p+1+N)]### add discrete ----
    }else{
      phi_OUT[c(N/2,N/2+2),'SE'] <- sd.DR[(num.a+num.g+num.p+2):(num.a+num.g+num.p+2+1)]### add discrete ----
      if(N>2){
        phi_OUT[c(1:(N/2-1)),'SE'] <- phi_OUT[N/2,'SE'] * ((N/2):2)
        phi_OUT[c((N/2+3):(N+1)),'SE'] <- phi_OUT[N/2+2,'SE'] * (2:(N/2))
      }
    }

    phi_OUT[-c(N/2+1),-(1:3)] <- CI.sig(par=v.Domain[-c(N/2+1)], Sigma=phi_OUT[-c(N/2+1),'SE'],
                                        level = level, e = FALSE, p = FALSE)
    phi_OUT[,'prob'] <- dbinom(x = 0:N, size = N, prob = pi)
    name.phi <- c(paste('phi_',(N/2):1,sep=''),'0', paste('psi_',1:(N/2),sep=''))
    rownames(phi_OUT) <- name.phi
    name.phi <- c('pi',name.phi[-c(N/2+2)])
    if(slope.model != 'free-Willy'){
      name.phi <- name.phi[c(1, 1+N/2, 1+N/2+1)]
    }
  }

  T.dist <- GG_tree(sigma = sigma, eta = eta, nu = nu, Q = Q, GG.idx = GG.idx, names.strata = strata)
  if(no.frail){
    Z.dist <- NULL
  }else{
    Z.dist <- which_dist(par = model_par,ag.idx=ag.idx, a.idx=a.idx, g.idx=g.idx,
                         PVF.idx=PVF.idx,IG.idx=IG.idx,frailty_names=RE_stratum,type=type,point_mass=p)
    if(num.p>0) Z.dist$shifted <- list(
      paste('Additionally domain shifted by',paste(p,collapse = ' ')), p = p )### add discrete ----
  }

  marg.logli <- marginal$value/scale.marginal
  marg.logli.nolam <- marginal$value.nolam/scale.marginal
  cond.logli <- cond$value/scale.cond
  AIC <- aic(loglihood = marg.logli, num_par = num.a+num.g+num.p+K+num.lam+num.lam.p+slope.par)
  AIC <- c(AIC,
           aic(loglihood = marg.logli.nolam, num_par = num.a+num.g+num.p+K+num.lam+num.lam.p+slope.par))
  if(!no.frail & !current.status) names(AIC) <- c('full', 'no sum(lam)')

  name.theta <- row.names(theta_OUT)
  if(any(IG.idx)){
    for(j in which(IG.idx)){
      tmp <- ifelse(is.finite(max(alpha.rank[1:j])), max(alpha.rank[1:j]), 0) + sum(IG.idx[1:j])
      name.theta[tmp] <- NA
    }
    name.theta <- name.theta[!is.na(name.theta)]
  }

  row.names(H) <- colnames(H) <- row.names(Sigma2.DR) <- colnames(Sigma2.DR) <- names(S) <-
    c(name.theta, name.p, name.phi, row.names(beta_OUT), tmp.name2, tau.name)

  if(any(GG.idx)){
    cat('Baseline Hazard: \n')
    print(Hazard)
    cat('\n')
  }

  if(K>0){
    cat('beta: \n')
    print(beta_OUT)
    cat('\n')

    cat('HR: \n')
    print(HR_OUT)
    cat('\n')
  }

  if(!no.frail){
    cat('theta: \n')
    print(theta_OUT)
    if(num.p>0){
      cat('point-mass: \n')
      print(p_OUT)
    }
    cat('\n')
  }


  if(slope){
    cat('phi, psi: \n')
    print(phi_OUT)
    cat('\n')
    cat('V* ~ Binom(',N,',',pi,')', sep = '')
    cat('\n')
  }

  if(!no.frail){
    cat('(full) univariate loglihood = ', uni.logli, '\n')
    if(!(current.status & !cs.options$EM.cs)) cat('conditional loglihood', cond.logli , '\n')
    cat('(full) marginal loglihood = ', marg.logli , '\n')
  }else{
    uni.logli <- NULL
    cat('univariate loglihood', cond.logli , '\n')
  }

  if(overdisp){
    cat(paste(RE_stratum,sep=','), 'overdispersion parameter(s)',
        paste('omega = ', round(tau['omega',],3), collapse = ';'),'\n')
  }

  if(iter==(iterlim+1)) print('Iteration limit exceeded!')

  OUT <- list(theta = theta_OUT, p = p_OUT, beta = beta_OUT, lam_0 = Hazard, model.par = model_par, v.Domain = v.Domain,
              Cov = Sigma2.DR, H = H, S = S, z = z, v=v, T.dist = T.dist, Z.dist = Z.dist,
              baseline = baseline, marginal = marginal, cond = cond, AIC = AIC, tau = tau,
              marginal.loglihood = marg.logli, conditional.loglihood = cond.logli, univariate.loglihood = uni.logli,
              iter  = iter, plot.predict = list(X = X, beta = beta, strata=strata,ag.idx=ag.idx,a.idx=a.idx,partial.lambda=partial.lambda,
                                                breaks=breaks,RE_stratum=RE_stratum,dt = dt_1, tr.dt = tr.dt_1, GG.stratum = GG.stratum, partial_i = partial_i))### add discrete. 22.02.24:dt and tr_dt added for oos cum hazard computation ----
  return(OUT)
}
