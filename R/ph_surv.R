#' Fit univariate or multivariate time-to-event models with Addams Family or PVF frailty or no frailty
#' with piecewise constant, Breslow or generalized gamma baseline hazard and proportional hazard factors.
#'
#' @param data data.frame. Only numeric values or logical or factors are allowed (latter will be coerced to numeric).
#' Note that multilevel factors as proportional hazard covariates are not yet supported
#' @param X vectors of covariate (column) names referring to data
#' @param strata string that refers to stratum variable (column) in data. Can be mutli-level.
#' Null if no stratified hazard
#' @param y name of time variable in data
#' @param d name of status indicator in data (should be 1 for event, 0 for censored)
#' @param frailty name stratum variable for frailty distrbution (models for frailty distribution).
#' NULL if frailty distribution is without model or no frailty mode is fitted
#' @param ID name of ID variable to indicate clusters
#' @param tr name of left-truncation time variable. Only used for right-censored data! can be 0 or NULL in case of no truncation.
#' @param weights name of weight variable. weighting factors twhitch which individual likelihood contribution is multiplied with
#' @param start_dist starting values of frailty parameters
#' @param start_lnlam starting values of hazard parameters (note that in some cases log hazard parameters are used, e.g. piecewise hazard)
#' @param start_beta starting values of PH coeffiecients
#' @param current.status logical. If TRUE current status data models will be fitted, right censored data if FALSE
#' @param Random.Slope logical (if random slopes should be imposed. experimental. Better ignore)
#' @param slope.model model for random slope. Was experimental. Should be ignored
#' @param binom.n if aB or B+p is chosen in type the this is the number of RC chosen.
#' Might also be number of RC for binomial random slopes (which is experimental, better don't model random slopes and better stratify frailty distribution by option frailty)
#' @param start_slope starting values ranodm slopes
#' @param start_pi starting values random slopes
#' @param baseline list up to length 4. Specifies the baseline hazard model.
#' Default is baseline = list('piecewise','Breslow',NULL,NULL)
#' The first element of the list can be a vector of length strata specifying the baseline hazards.
#' Options are piecewise, GG, and its special cases
#' c('gamma', 'gamma^-1', 'Weibull', 'Weibull^-1', 'exponential', 'exponential^-1', 'ammag', 'ammag^-1', 'lognormal', 'location', 'N/2', 'GGpos', 'GGinv' ).
#' Bseline hazard distibutions can differ for strata.
#' The second element of the list specifies handling of piecewise hazards. Can be Intervall or Breslow. In case of intervall the third element of the list is the step size of intervals.
#' Alternatively breaks.piece can be specified for customized time intervals for dsitnct hazard parameters.
#' Fourth elemnt of list specifies a groups for which all GG parameters except for the location parameter are shared.
#' Only relevant if location is chosen as baseline hazard (first element of list).
#' Groups are numbered strata, can be more than one group,
#' i.e. fourth element can be a list.
#' #' @param type c('alpha_gamma', 'alpha', 'gamma','PVF', 'cP', 'IG', 'H')
#' Frailty distibution (first is Addams, second scaled Poisson).
#' Can be of length(unique(data$frailty)), i.e. differnet distribution for varying strata can be chosen.
#' @param ci_level level on CIs
#' @param overdisp.cs logical. If true overdispersion model will be estimated (Compound multinomial Dirichilet).
#' Only available for current status data.
#' @param scaleX logical. Should covariates be scaled
#' @param centerX logical, should covariate be centered
#' @param breaks list of break points for distinct hazard parameters for piecewise hazard
#' @param iterlim maximum number of iterations
#' @param reltol
#' @param converge convergence threshold. Note that convergence criterion is abolute not relative increment in loglihood
#' @param method optimisazion mehtod. Essentiall everything that R hast to offer
#' can be chosen (optimize, every method of optim, nlm, nlminb, DEoptim, constrOptim, genoud...)
#' @param c.method same as method but for conditional model
#' @param m.method same as above butt for marignal model
#' @param Richardson logical. Indicates whther Richardson derivative as implementd in numDeriv should be used
#' @param R.options list of options for numDeriv partameters.
#' Exact naing is a must. Parameters that are not set here are default values
#' @param m.options list of options for the optimizer of maringal model.
#' Names must be exact the same as the options of chosen optimzation routines. All options that are not set will be set to default values.
#' @param c.options same as m.options for conditional model
#' @param multi.num whther numerical derivatives should be used.
#' For bivariate tight censored models derivatives and (partially) hessians are analytically derived
#' @param cs.options options for current statuts data optimization.
#' List list(EM.cs, univariate,uni.model=NA,uni.initial), First states whether EM algorithm should be used (set to FALSE, EM is horrible in this case), second whether univariate model should be also fitted, third can be the univariate model such that it does not have to be omputed more than once if you fit more than one model, fourth whether hazard parameters of maringal model should be initialized with univariate model (my impression is that this often works better for current status data).
#' @param frail.thresh threshold for borderline cases, i.e. when should we assume that the AF model is the gamma, Poisson etc model. Might be important for optimization purposes as this changes how survival function is computed in current optimization iteration
#' @param H.singular treatment if Hessian is singualr. Can be gcholesky, ginverse, standard.
#' @param scale.marginal scaling of maringal loglihood. Actually makes no sense.
#' Can better be handled via sspecifiying options for optimization routines.
#' @param scale.cond same as above
#' @param z initinal values of the frailties
#' @param v intial values of the random slopes
#' @param initial.z logical indicator. Specifies whther frailties should be intialized with intial frailty distribution paramers. If FALSE frailties are 1 in first iteration (only relevant for right censored data)
#' @param prnt logical. indicates whther current iteration should be printed.
#'
#' @return list of parameter estimates, choice of distribution and further estimates for plotting and prediction purposes
#' @export
#'
#' @examples 'Please see R markdown script that is also uploaded on git.
#' Documentaion and vignette is work in progress...'
ph.surv <- function(data, X = NULL, strata = NULL, y, d, frailty = NULL, ID = NULL, tr = NULL,
                    weights = NULL,start_dist = NA, start_lnlam = NA, start_beta = NA,
                    current.status = FALSE, Random.Slope = NA, slope.model = 'free-Willy', binom.n = 0, start_slope = NA, start_pi= NA,
                    baseline = list('piecewise','Breslow',NULL,NULL), type, ci_level = 0.95,
                    overdisp.cs = FALSE, scaleX = TRUE, centerX = TRUE, breaks = NULL,
                    iterlim = 1000, reltol = 1e-10, converge  = 1e-6, method = 'BFGS',
                    c.method = method, m.method = method, Richardson = FALSE, R.options = list(),
                    m.options = list(), c.options = list(),multi.num = FALSE,
                    cs.options = list(EM.cs=FALSE, univariate = TRUE, uni.model = NA, uni.initial=TRUE),
                    frail.thresh=1e-6, H.singular = 'standard',
                    scale.marginal = -1, scale.cond = -1, z = NA, v = NA, initial.z = TRUE,prnt = FALSE){

  # plausability checks and checking formats
  if('nofrail' %in% type){
    if(length(type) > 1){
      stop('No frailty model cannot be mixed with others.')
    }
  }

  if(!is.list(baseline)){
    baseline <- as.list(baseline)
  }

  if('Breslow' %in% baseline[[1]] & any(baseline[[1]] != 'Breslow')){
    stop('Either all hazards are non-parametric or none!
         Consider piecewise constant hazard with Breslow-style')
  }

  # optimax & Richardson options
  if(!('Breslow' %in% baseline[[1]]) | is.null(X)){
    c.options$Richardson <- Richardson
    if(c.method == 'nlm'){
      # nlm, elements not listed cannot be specified
      c.options$print.level <- ifelse(is.null(c.options$print.level),0,c.options$print.level)
      c.options$ndigit <- ifelse(is.null(c.options$ndigit),12,c.options$ndigit)
      c.options$gradtol <- ifelse(is.null(c.options$gradtol),reltol,c.options$gradtol)
      c.options$steptol <- ifelse(is.null(c.options$steptol),1e-6,c.options$steptol)
      c.options$iterlim <- ifelse(is.null(c.options$iterlim),iterlim,c.options$iterlim)
      c.options$typsize <- ifelse(is.null(c.options$typsize),1,c.options$typsize)
      # stepmax follows under initialization of parameters in maxll
    }else if(c.method %in% c("Nelder-Mead", "BFGS", "CG",
                             "L-BFGS-B", "SANN","Brent",
                             "constr.optim", "genoud")){
      # optim
      if(is.null(c.options$control)) c.options$control <- list()
      c.options$control$maxit <- ifelse(is.null(c.options$control$maxit),iterlim,c.options$control$maxit)
      c.options$control$reltol <- ifelse(is.null(c.options$control$reltol),reltol,c.options$control$reltol)
      c.options$lower <- ifelse(is.null(c.options$lower),-Inf,c.options$lower)
      c.options$upper <- ifelse(is.null(c.options$upper),Inf,c.options$upper)
      if(c.method == 'constr.optim'){
        warning('constr.optim ist noch nicht getestet')
        # ci, ui must be set by user in c.options!
        c.options$optim.method <- ifelse(is.null(c.options$optim.method),'BFGS',
                                         c.options$optim.method)
        c.options$mu <- ifelse(is.null(c.options$mu),1e-04,c.options$mu)
        c.options$outer.iterations <- ifelse(is.null(c.options$outer.iterations),
                                             100,c.options$outer.iterations)
        c.options$outer.eps <- ifelse(is.null(c.options$outer.eps),
                                      1e-05,c.options$outer.eps)
      }
      if(method == 'genoud'){
        # nvars has to be set further below
        gen <- list(pop.size=1000, max.generations=100,
                    wait.generations=10, hard.generation.limit=TRUE, starting.values=NULL,
                    MemoryMatrix=TRUE, Domains=NULL, default.domains=10,
                    solution.tolerance=0.001, boundary.enforcement=0, lexical=FALSE,
                    gradient.check=TRUE, BFGS=TRUE, data.type.int=FALSE, hessian=FALSE,
                    unif.seed=round(stats::runif(1, 1, 2147483647L)),
                    int.seed=round(stats::runif(1, 1, 2147483647L)),print.level=2, share.type=0,
                    instance.number=0,
                    project.path=NULL, P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50,
                    P8=50, P9=0, P9mix=NULL, BFGSburnin=0, BFGSfn=NULL, BFGShelp=NULL,
                    control=list(),
                    transform=FALSE, debug=FALSE, cluster=FALSE, balance=FALSE)
        gen$optim.method <- ifelse(gen$boundary.enforcement < 2, "BFGS",
                                   "L-BFGS-B")
        for(i in names(gen)){
          if(is.null(c.options[[i]])) c.options[[i]] <- gen[[i]]
        }
      }
    }else if(c.method == 'optimize'){
      # optimize
      c.options$lower <- ifelse(is.null(c.options$lower),-1e2,c.options$lower)
      c.options$upper <- ifelse(is.null(c.options$upper),1e2,c.options$upper)
      c.options$tol <- ifelse(is.null(c.options$tol),reltol,c.options$tol)
      c.options$tol <- reltol
    }else if(c.method == 'nlminb'){
      if(is.null(c.options$control)) c.options$control <- list(eval.max=200,iter.max=iterlim,trace=0,rel.tol=reltol)
      c.options$control$iter.max <- ifelse(is.null(c.options$control$iter.max),iterlim,c.options$control$iter.max)
      c.options$control$rel.tol <- ifelse(is.null(c.options$control$rel.tol),reltol,c.options$control$rel.tol)
      c.options$lower <- ifelse(is.null(c.options$lower),-1e4,c.options$lower)
      c.options$upper <- ifelse(is.null(c.options$upper),1e4,c.options$upper)
    }else if(c.method %in% c('DEoptim')){
      warning('Bei DEoptim ist irgendwas in Argument falsch vermute ich')
      # DEoptim
      if(is.null(c.options$control)) c.options$control <- list()
      tmp <- stats::DEoptim.control(itermax = iterlim, reltol = reltol)
      for(i in names(tmp)){
        if(is.null(c.options$control[[i]])) c.options$control[[i]] <- tmp[[i]]
      }
      c.options$lower <- ifelse(is.null(c.options$lower),-100,c.options$lower)
      c.options$upper <- ifelse(is.null(c.options$upper),100,c.options$upper)
      # has to be adjusted to parameter vector length
      # fnMap is either set or NULL is correct argument
      rm(list = 'tmp')
    }

  }
  c.options$method <- c.method

  if(!('nofrail' %in% type)){
    m.options$Richardson <- Richardson
    if(m.method == 'nlm'){
      # nlm, elements not listed cannot be specified
      m.options$print.level <- ifelse(is.null(m.options$print.level),0,m.options$print.level)
      m.options$ndigit <- ifelse(is.null(m.options$ndigit),12,m.options$ndigit)
      m.options$gradtol <- ifelse(is.null(m.options$gradtol),reltol,m.options$gradtol)
      m.options$steptol <- ifelse(is.null(m.options$steptol),1e-6,m.options$steptol)
      m.options$iterlim <- ifelse(is.null(m.options$iterlim),iterlim,m.options$iterlim)
      m.options$typsize <- ifelse(is.null(m.options$typsize),1,m.options$typsize)
      # stepmax follows under initialization of parameters in maxll
    }else if(m.method %in% c("Nelder-Mead", "BFGS", "CG",
                             "L-BFGS-B", "SANN","Brent",
                             "constrOptim", "genoud")){
      # optim
      if(is.null(m.options$control)) m.options$control <- list()
      m.options$control$maxit <- ifelse(is.null(m.options$control$maxit),iterlim,m.options$control$maxit)
      m.options$control$reltol <- ifelse(is.null(m.options$control$reltol),reltol,m.options$control$reltol)
      m.options$lower <- ifelse(is.null(m.options$lower),-Inf,m.options$lower)
      m.options$upper <- ifelse(is.null(m.options$upper),Inf,m.options$upper)
      if(m.method == 'constr.optim'){
        # ci, ui must be set by user in m.options!
        m.options$optim.method <- ifelse(is.null(m.options$optim.method),'BFGS',
                                         m.options$optim.method)
        m.options$mu <- ifelse(is.null(m.options$mu),1e-04,m.options$mu)
        m.options$outer.iterations <- ifelse(is.null(m.options$outer.iterations),
                                             100,m.options$outer.iterations)
        m.options$outer.eps <- ifelse(is.null(m.options$outer.eps),
                                      1e-05,m.options$outer.eps)
      }
      if(method == 'genoud'){
        # nvars has to be set further below
        for(i in names(gen)){
          if(is.null(m.options[[i]])) m.options[[i]] <- gen[[i]]
        }
        rm(list = 'gen')
      }
    }else if(m.method == 'optimize'){
      # optimize
      m.options$lower <- ifelse(is.null(m.options$lower),-1e2,m.options$lower)
      m.options$upper <- ifelse(is.null(m.options$upper),1e1,m.options$upper)
      m.options$tol <- ifelse(is.null(m.options$tol),reltol,m.options$tol)
      m.options$tol <- reltol
    }else if(m.method == 'nlminb'){
      if(is.null(m.options$control)) m.options$control <- list(eval.max=200,iter.max=iterlim,trace=0,rel.tol=reltol)
      m.options$control$iter.max <- ifelse(is.null(m.options$control$iter.max),iterlim,m.options$control$iter.max)
      m.options$control$rel.tol <- ifelse(is.null(m.options$control$rel.tol),reltol,m.options$control$rel.tol)
      m.options$lower <- ifelse(is.null(m.options$lower),-1e4,m.options$lower)
      m.options$upper <- ifelse(is.null(m.options$upper),1e4,m.options$upper)
    }else if(c.method %in% c('DEoptim')){
      # DEoptim
      if(is.null(m.options$control)) m.options$control <- list()
      tmp <- DEoptim::DEoptim.control(itermax = iterlim, reltol = reltol)
      for(i in names(tmp)){
        if(is.null(m.options$control[[i]])) m.options$control[[i]] <- tmp[[i]]
      }
      m.options$lower <- ifelse(is.null(m.options$lower),-100,m.options$lower)
      m.options$upper <- ifelse(is.null(m.options$upper),100,m.options$upper)
      # has to be adjusted to parameter vector length
      # fnMap is either set or NULL is correct argument
      rm(list = 'tmp')
    }
  }
  m.options$method <- m.method

  R.options$eps <- ifelse(is.null(R.options$eps),1e-4,R.options$eps)
  R.options$d <- ifelse(is.null(R.options$d),0.0001,R.options$d)
  R.options$zero.tol <- ifelse(is.null(R.options$zero.tol),sqrt(.Machine$double.eps/7e-7),
                               R.options$zero.tol)
  R.options$r <- ifelse(is.null(R.options$r),4,R.options$r)
  R.options$v <- ifelse(is.null(R.options$v),2,R.options$v)

  n.date <- NULL
  if(overdisp.cs){
    if(!current.status){
      stop('Overdispersion can only implemented for current status data!')
    }
    if('nofrail' %in% type){
      stop('Overdispersion only implemented for frailty model!')
    }
    if(cs.options$EM.cs){
      stop('Overdispersion not implemented for EM algorithm')
    }
    if(any(dim(unique(data[,c(X,strata,y,d,frailty,ID,tr,weights)])) !=
           dim(data[,c(X,strata,y,d,frailty,ID,tr,weights)]))){
      warning('Data must be aggregated for overdispersion model!
           You might want to use data_weights beforehand.
              Function input has to include weights!')
      data <- data_weights(data=data, X = X, strata = strata, y=y, d=d, frailty = frailty,
                           Random.Slope = NA,ID = ID, tr = tr, weights = weights, keep.names = TRUE)
    }
    data2 <- data
    data2[,d] <- 1
    n.date <- data_weights(data=data2, X = X, strata = strata, y=y, d=d, frailty = frailty,
                           Random.Slope = NA,ID = ID, tr = tr, weights = weights, keep.names = TRUE)[,c(ID,frailty,weights)]
    n.date[,c(ID)] <- as.numeric(n.date[,c(ID)])
    n.date[,c(frailty)] <- as.numeric(n.date[,c(frailty)])
    n.date <- as.matrix(n.date)
    n.date <- sapply(X = unique(n.date[,ID]),FUN = function(i){n.date[n.date[,ID]==i,c(frailty,weights),drop=FALSE][1,]}, simplify = TRUE)
    if(is.vector(n.date)){
      n.date <- cbind(1,n.date)
    }else{
      n.date <- t(n.date)
    }
    colnames(n.date) <- c('frailty','n')
    rm(list = c('data2'))
  }

  if(is.null(ID)){
    ID <- 1:dim(data)[1]
    if(!('nofrail' %in% type)){
      stop('You have to set an ID variable!')
    }
  }else{
    ID <- data[,ID]
  }
  if(!is.factor(ID)){
    ID <- factor(ID)
  }
  if(length(unique(ID)) != length(levels(ID))){
    warning('Level in ID has been dropped as data not present!')
    ID <- droplevels(x = ID,
                     exclude = levels(ID)[!(levels(ID)%in%unique(ID))])
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
    sd_x <- NULL
    mean_x <- NULL
  }else{
    K <- length(X)
    X <- data[,X, drop = FALSE]
    sd_x <- rep(1, K)
    mean_x <- rep(0,K)
    for(k in 1:K){
      if(is.factor(X[,k])){
        X[,k] <- as.numeric(X[,k]) - 1
      }else{
        stdX <- scale(x = X[,k], center = centerX, scale = scaleX)
        if(scaleX) sd_x[k] <- attr(stdX, 'scaled:scale')
        if(centerX) mean_x[k] <- attr(stdX, 'scaled:center')
        X[,k] <- stdX
      }
    }
    X <- as.matrix(X)
  }
  attr(X,'sd') <- sd_x# 24.02.2024
  attr(X,'mean') <- mean_x# 24.02.2024

  slope <- FALSE
  X.tilde <- 0
  if(!is.na(Random.Slope)){
    slope <- TRUE
    multi.num <- TRUE
    if(Random.Slope %in% colnames(X)){
      X.tilde <- X[,Random.Slope]
    }else{
      X.tilde <- data[,Random.Slope]
      X.tilde <- abs(X.tilde)
    }

    if(binom.n != round(binom.n)){
      stop('binom.n must be an even integer! There are as many negative as positive values on the domain.')
    }
    if(current.status){
      stop('Random slope not implemented for curent status data.')
    }

    if(!(slope.model %in% c('free-Willy', 'proportional'))) stop('slope.model must be \"free-Willy\" or \"proportional\"!')
  }

  if(is.null(strata)){
    stratum <- factor(x = (y*0+1), levels = 1, labels = 'base')
  }else{
    stratum <- data[,strata]
  }
  if(!is.factor(stratum)){
    stratum <- factor(x = stratum, levels = unique(stratum),
                      labels = paste('Str.', unique(stratum)))
    warning('Highly advised to set factors (for strata) oneself as order of
             baseline[[1]] has to be in order of levels.')
  }
  if(length(unique(stratum)) != length(levels(stratum))){
    warning('Level has been dropped as not present in dataset. Better do this yourself!')
    stratum <- droplevels(x = stratum,
                          exclude = levels(stratum)[!(levels(stratum)%in%unique(stratum))])
  }

  if(is.null(frailty)){
    frailty <- factor(x = (y*0+1), levels = 1, labels = 'frail')
  }else{
    frailty <- data[,frailty]
  }
  frailty <- tapply(X = frailty, INDEX = ID, FUN = function(x) x[1], simplify = FALSE)
  frailty <- sapply(X = frailty, FUN = function(x) x)
  if(!is.factor(frailty)){
    frailty <- factor(x = frailty, levels = unique(frailty),
                      labels = paste('Pop.', unique(frailty)))
    warning('Highly advised to set factors (for frailty) oneself as order of
             type has to be in order of levels.')
  }
  if(length(unique(frailty)) != length(levels(frailty))){
    warning('Level has been dropped as not present in dataset. Better do this yourself!')
    frailty <- droplevels(x = frailty,
                          exclude = levels(frailty)[!(levels(frailty)%in%unique(frailty))])
  }

  if(is.null(weights)){
    c.weights <- y*0+1
  }else{
    c.weights <- data[,weights]
  }
  m.weights <- tapply(X = c.weights, INDEX = ID, FUN = function(x) x[1], simplify = FALSE)
  m.weights <- sapply(X = m.weights, FUN = function(x) x)

  OUT <- maxll(start_dist = start_dist, start_lnlam = start_lnlam, start_beta = start_beta,
               baseline = baseline, K = K, current.status = current.status, v = v, slope.model = slope.model,
               type = type, level = ci_level, z = z, slope = slope, X.tilde = X.tilde, N=binom.n, overdisp=overdisp.cs,n.date=n.date,
               y = y, d = d, X = X, ID = ID, tr = tr, frailty = frailty, start_slope = start_slope,
               start_pi = start_pi, sd_x = sd_x, stratum = stratum, breaks = breaks, initial.z = initial.z,
               iterlim = iterlim, converge  = converge, multi.num = multi.num,cs.options=cs.options,
               H.singular = H.singular, c.options = c.options, m.options = m.options, R.options = R.options,
               scale.marginal = scale.marginal, scale.cond = scale.cond,frail.thresh=frail.thresh,
               c.weights = c.weights, m.weights = m.weights, num.score = num.score, num.hessian  = num.hessian, prnt = prnt)

  return(OUT)

}
