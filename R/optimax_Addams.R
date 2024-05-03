
optimax <- function(par,fn,gr=NULL,Hessian=NULL,method,options,Richardson,R.options,...){

  if(method == 'nlm'){
    # nlm
    if(!is.null(gr)){
      m <- fn
      fn <- function(par,...){
        y <- m(par,...)
        attr(y,'gradient') <- gr(par, ...)
        if(!is.null(Hessian)){
          attr(y,'hessian') <- Hessian(par, ...)
        }
        return(y)
      }
    }else if(Richardson){
      m <- fn
      fn <- function(par,...){
        y <- m(par,...)
        del <- numDeriv::genD(func = m, x = par, method = 'Richardson', method.args = R.options,...)$D
        np <- length(par)
        gr <- del[1:np]
        del <- del[-c(1:np)]
        Hessian <- matrix(NA, nrow = np, ncol = np)
        Hessian[upper.tri(Hessian,diag = TRUE)] <- del
        tmp <- which(is.na(Hessian), arr.ind = TRUE)
        tmp2 <- cbind(tmp[,2],tmp[,1])
        Hessian[tmp] <- Hessian[tmp2]
        attr(y,'gradient') <- gr
        attr(y,'hessian') <- Hessian
        return(y)
      }
    }
    opt <- stats::nlm(f = fn, p = par, ...,hessian = FALSE,
               fscale = 1,print.level = options$print.level,ndigit = options$ndigit,
               gradtol = options$gradtol,steptol = options$steptol, stepmax = options$stepmax,
               typsize = options$typsize, iterlim = options$iterlim,check.analyticals = FALSE)

    OUT <- list(par=opt$estimate,value=opt$minimum,counts=opt$iterations,
                convergence=opt$code,gradient=opt$gradient,method=method,
                note='Caution: nlm convergence code, not optim!')
  }else if(method %in% c("Nelder-Mead", "BFGS", "CG",
                         "L-BFGS-B", "SANN","Brent")){
    # optim
    OUT <- stats::optim(par = par,fn = fn,gr = gr, ...,method = method,lower = options$lower,
                 upper = options$upper,control = options$control,hessian = FALSE)
    OUT$method <- method

  }else if(method == 'optimize'){
    # optimize
    opt <- stats::optimize(f = fn, ...,lower = options$lower,upper = options$upper,tol = options$tol)
    OUT <- list(value=opt$objective,par=opt$minimum,counts=NA,
                convergence=NA,method=method)
  }else if(method == 'nlminb'){
    opt <- stats::nlminb(start = par, objective =fn, gradient = gr, hessian = Hessian, ...,
                  scale = 1, control = options$control, lower = options$lower, upper = options$upper)
    OUT <- list(par=opt$par,value=opt$objective,counts=opt$iterations,
                convergence=opt$convergence,evaluations=opt$evaluations,method=method,
                note=paste(opt$message,'Caution: nlminb convergence code, not optim!'))
  }else if(method == 'DEoptim'){
    opt <- DEoptim::DEoptim(fn = fn, lower = options$lower,upper = options$upper,control = options$control,...,
                   fnMap = options$fnMap)
    OUT <- list(value=opt$optim$bestval,par=opt$optim$bestmem,counts=opt$optim$iter,
                convergence=NA,method=method,DEoptim=opt)
  }else if(method == 'constrOptim'){
    opt <- stats::constrOptim(theta = par, f = fn, grad = gr, ui = options$ui, ci = ui$ci, mu = options$mu,
                       control = options$control,
                       method = options$optim.method,
                       outer.iterations = options$outer.iterations, outer.eps = options$outer.eps,
                       ...,hessian = FALSE)
    opt$method <- method
    opt$method.optim <- options$optim.method
    OUT <- opt
  }else if(method == 'genoud'){
    opt <- rgenoud::genoud(fn=fn, nvars=options$nvars, max=FALSE, pop.size=options$pop.size, max.generations=options$max.generations,
                  wait.generations=options$wait.generations, hard.generation.limit=options$hard.generation.limit, starting.values=par,
                  MemoryMatrix=options$MemoryMatrix, Domains=options$Domains, default.domains=options$default.domains,
                  solution.tolerance=options$solution.tolerance, gr=gr, boundary.enforcement=options$boundary.enforcement, lexical=options$lexical,
                  gradient.check=options$gradient.check, BFGS=options$BFGS, data.type.int=options$data.type.int, hessian=FALSE,
                  unif.seed=options$unif.seed,
                  int.seed=options$int.seed,print.level=options$print.level, share.type=options$share.type,
                  instance.number=options$instance.number, output.path="stdout", output.append=FALSE,
                  project.path=options$project.path, P1=options$P1, P2=options$P2, P3=options$P3, P4=options$P4, P5=options$P5, P6=options$P6, P7=options$P7,
                  P8=options$P8, P9=options$P9, P9mix=options$P9mix, BFGSburnin=options$BFGSburnin, BFGSfn=options$BFGSfn, BFGShelp=options$BFGShelp,
                  control=options$control, optim.method=options$optim.method, transform=options$transform, debug=options$debug, cluster=options$cluster, balance=options$balance,
                  ...)
    opt$optim.method <- options$optim.method
    opt$method <- opt$method
    opt$convergence <- NA
    OUT <- opt
  }


  return(OUT)
}
