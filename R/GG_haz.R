GG.haz <- function(s,e,nu,y,str,n){
  # hier d?rfen nur GG Individuen rein.

  Lambda_ind <- matrix(NA, nrow = n, ncol = 2)
  colnames(Lambda_ind) <- c('Lambda', 'lambda')

  tmp <- y>0
  Lambda_ind[!tmp,] <- c(rep(0,sum(!tmp)),rep(.Machine$double.xmin,sum(!tmp)))
  y <- y[tmp]
  str <- str[tmp]
  idx.pos <- nu > 0
  idx.0 <- nu == 0
  idx.neg <- nu < 0
  str.pos <- str[idx.pos[str]]
  str.0 <- str[idx.0[str]]
  str.neg <- str[idx.neg[str]]
  tmp <- 1:sum(tmp)
  tmp.pos <- tmp[idx.pos[str]]
  tmp.0 <- tmp[idx.0[str]]
  tmp.neg <- tmp[idx.neg[str]]

  #shape <- nu^-2
  #rate <- shape * exp(-e*nu/s)
  #x <- y^(nu[str]/s[str])
  #  lnconst <- log(abs(nu[str])) - log(s[str]) + log(y)*(nu[str]/s[str]-1)
  shape <- rate <- nu^-2
  x <- (exp(-e[str])*y)^(nu[str]/s[str])
  lnconst <- log(abs(nu[str])) - log(s[str]) + (log(y)-e[str])*(nu[str]/s[str]-1) - e[str]

  if(any(tmp.pos)){
    y.p <- x[tmp.pos]
    shape.p <- shape[str.pos]
    rate.p <- rate[str.pos]
    ln.p <- lnconst[tmp.pos]
    Lambda_ind[tmp.pos,1] <- - stats::pgamma(q = y.p, shape = shape.p,
                                      rate = rate.p,log.p = TRUE, lower.tail = FALSE)
    Lambda_ind[tmp.pos,2] <- exp( stats::dgamma(x = y.p, shape = shape.p,
                                         rate = rate.p, log = TRUE) +
                                    ln.p +
                                    Lambda_ind[tmp.pos,1] )
  }

  if(any(tmp.0)){
    y.0 <- y[tmp.0]
    s.0 <- s[str.0]
    e.0 <- e[str.0]
    Lambda_ind[tmp.0,1] <- -stats::plnorm(q = y.0, meanlog = e.0, sdlog = s.0,
                                   log.p = TRUE, lower.tail = FALSE)
    Lambda_ind[tmp.0,2] <- exp( stats::dlnorm(x = y.0, meanlog = e.0, sdlog = s.0, log = TRUE) +
                                  Lambda_ind[tmp.0,1] )
  }

  if(any(tmp.neg)){
    y.n <- x[tmp.neg]
    shape.n <- shape[str.neg]
    rate.n <- rate[str.neg]
    ln.n <- lnconst[tmp.neg]
    Lambda_ind[tmp.neg,1] <- - stats::pgamma(q = y.n, shape = shape.n,
                                      rate = rate.n,log.p = TRUE)
    Lambda_ind[tmp.neg,2] <- exp( stats::dgamma(x = y.n, shape = shape.n,
                                         rate = rate.n, log = TRUE) +
                                    ln.n +
                                    Lambda_ind[tmp.neg,1] )
  }

  return(Lambda_ind)

}
