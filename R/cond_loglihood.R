cond.loglihood <- function(par,num.lam,K,scale.cond=-1,GG.op,GG.sc,GG.special,X.GG,y.GG,
                           tr.GG,z.GG,d.GG,n.GG,X.p,z.p,d.p,n.p,Q.GG,stratum.GG,dt,m,tr.dt){

  beta <- par[1:K]
  lnlambda <- par[(K+1):(K+num.lam)]

  if(n.p>0){
    cox <- partial.loglihood(beta=beta, z = z.p, d = d.p, X = X.p, dt = dt,
                             tr.dt = tr.dt, m = m, n = n.p, K = K, scale.cond = 1)
  }else{
    cox <- 0
  }

  full <- GG.loglihood(beta=beta, lnlambda = lnlambda, z = z.GG, d = d.GG, X = X.GG, y = y.GG,
                       tr = tr.GG, n = n.GG, Qstar = Q.GG, stratum = stratum.GG,
                       GG.op = GG.op, GG.sc = GG.sc, GG.special = GG.special)

  loglihood <- full + cox

  return(scale.cond*loglihood)
}
