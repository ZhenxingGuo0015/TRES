### functions to solve NB problem
M6Apeak <- function(mat,
                    sf = NULL,
                    cutoff = NULL,
                    update = "Joint",
                    trans = NULL,
                    optM = "L-BFGS-B",
                    myfscale = -1e+6){
  ### get statistic for each peak
  require(stats)
  mu.hat = EstiMu(counts = mat, sf = sf)
  para = EstiPhi(counts = mat, sf = sf, update, trans = trans,
                 optM = optM, myfscale = myfscale)
  mu.var = EstiVarofMu(counts = mat,
                       sf = sf,
                       mu = mu.hat,
                       phi = para$phi,
                       theta = para$theta)
  if(is.null(cutoff)){
    ix1 = which(rowMeans(mat) > 200)
    ix2 = which(!is.na(mu.hat))
    ix = intersect(ix1, ix2)
    tmp = mix.2norm(mu.hat[ix], var.eq=FALSE)
    cutoff = (tmp$mu1+tmp$mu2) / 2
  }
  tstat = getStats(mu = mu.hat, mu0 = cutoff, var.mu = mu.var)
  pval = 1 - pnorm(tstat)
  return(data.frame(mu = mu.hat,
                    mu.var = mu.var,
                    stats = tstat,
                    shrkPhi = para$phi,
                    pvals = pval,
                    shrkTheta = para$theta
                    ))
}


EstiMu <- function(counts, sf){
  ### function to estimate mu
  ### counts: expression count, input1, ip1, input2, ip2, input3, ip3, ...
  require(stats)
  x = counts[, seq(1, ncol(counts), 2)] # input
  y = counts[, seq(2, ncol(counts), 2)] # ip
  n = ncol(y)
  if(is.null(sf)){
    sf = colSums(counts, na.rm = TRUE)/median(colSums(counts, na.rm = TRUE))
  }
  sf.x = sf[seq(1, ncol(counts), 2)]
  sf.y = sf[seq(2, ncol(counts), 2)]
  lambda_x.hat = sweep(x, 2, sf.x, FUN = "/")
  lambda_y.hat = sweep(y, 2, sf.y, FUN = "/")

  #### version 1
  # p.hat = lambda_y.hat/(lambda_x.hat + lambda_y.hat)
  # mu.hat = rowMeans(p.hat)
  #
  #### version 2
  mu.hat = rowMeans(lambda_y.hat)/rowMeans(lambda_x.hat + lambda_y.hat)

  return(mu.hat)
}

EstiPhi <- function(counts, sf, update = c("OnlyPhi", "Iterative", "Joint"), trans, optM, myfscale){
  require(stats)
  x = counts[, seq(1, ncol(counts), 2)]
  y = counts[, seq(2, ncol(counts), 2)]
  n = ncol(y)
  if(is.null(sf)){
    sf = colSums(counts)/median(colSums(counts))
  }
  sf.x = sf[seq(1, ncol(counts), 2)]
  sf.y = sf[seq(2, ncol(counts), 2)]
  lambda_x.hat = sweep(x, 2, sf.x, FUN = "/")
  lambda_y.hat = sweep(y, 2, sf.y, FUN = "/")
  p.hat = lambda_y.hat/(lambda_x.hat + lambda_y.hat)
  mu.hat = rowMeans(p.hat)

  require(matrixStats)
  phi.mom = ((n-1)/n)*rowVars(p.hat)/(mu.hat*(1 - mu.hat))
  theta.mom = rowMeans(lambda_x.hat + lambda_y.hat)/(phi.mom^{-1} -1)
  phi.bar = mean(phi.mom, na.rm = TRUE)
  s2phi = var(phi.mom, na.rm = TRUE)
  mlphi = log(phi.bar + 0.00001) - (1/2)*log(s2phi/(phi.bar)^2 + 1)
  sdlphi = sqrt(log(s2phi/(phi.bar)^2 + 1))

  ### bayes estimate of phi and theta
  update <- match.arg(update)
  switch(update,
         OnlyPhi = EstiPhi.SingleUpdate(y = y, x = x, sy = sf.y, sx = sf.x,
                                        mu = mu.hat, theta = theta.mom,
                                        mlphi = mlphi, sdlphi = sdlphi),
         Iterative = EstiPhi.IterativeUpdate(y = y, x = x, sy = sf.y, sx = sf.x,
                                             mu = mu.hat, theta0 = mean(theta.mom[theta.mom < Inf],
                                                                        na.rm = TRUE),
                                             mlphi = mlphi, sdlphi = sdlphi),
         Joint = EstiPhi.JointUpdate(y, x, sy = sf.y, sx = sf.x,
                                     mu = mu.hat, mlphi = mlphi, sdlphi = sdlphi,
                                     trans = trans, optM = optM, myfscale = myfscale))
}


EstiPhi.JointUpdate <- function(y, x, sy, sx, mu, mlphi, sdlphi, trans, optM, myfscale,
                                eps = 1e-05){
  ### function to get estimate of phi and theta with optim function
  require(stats)
  ### log likelihood
  log.lik <- function(para, yy, xx, mmu, sx, sy, mlphi, sdlphi){
    s = para[1]
    t = para[2]
    loglik = sum(dnbinom(yy, size = mmu*(s^{-1} - 1),
                         prob = 1/(1+sy*t), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(s^{-1} - 1),
                           prob = 1/(1+sx*t), log = TRUE)
    ) + dlnorm(s, mean = mlphi, sd = sdlphi,
               log = TRUE)
    loglik
  }

  log.lik_1 <- function(para, yy, xx, mmu, sx, sy, mlphi, sdlphi, trans){
    s = para[1]
    t = para[2]
    if(trans == "sin"){
      loglik = sum(dnbinom(yy, size = mmu*(((sin(s)+1)/2)^{-1} - 1),
                           prob = 1/(1+sy*exp(t)), log = TRUE) +
                     dnbinom(xx, size = (1 - mmu)*(((sin(s)+1)/2)^{-1} - 1),
                             prob = 1/(1+sx*exp(t)), log = TRUE)
      ) + dlnorm((sin(s)+1)/2, mean = mlphi, sd = sdlphi,
                 log = TRUE)
    }else if(trans == "exp"){
      loglik = sum(dnbinom(yy, size = mmu*((exp(-exp(s)))^{-1} - 1),
                           prob = 1/(1+sy*exp(t)), log = TRUE) +
                     dnbinom(xx, size = (1 - mmu)*((exp(-exp(s)))^{-1} - 1),
                             prob = 1/(1+sx*exp(t)), log = TRUE)
      ) + dlnorm(exp(-exp(s)), mean = mlphi, sd = sdlphi,
                 log = TRUE)
    }
    loglik
  }

  ### estimate with optim function
  res = matrix(0, nrow = nrow(y), ncol = 4)
  ix = which(is.na(mu))
  for(i in 1:length(mu)){
    # cat(i, sep = "\n")
    if(!is.na(mu[i])){
      ## initial (-0.0001, 0.0001)
      if(optM != "L-BFGS-B"){
        tmp = optim(c(-0.01, 0.001), log.lik_1,   #(-0.01, 0.001), (-0.0001, 0.0001)
                    yy = y[i, ], xx = x[i, ],
                    mmu = mu[i], sx = sx, sy = sy,
                    mlphi = mlphi, sdlphi = sdlphi,
                    trans = trans,
                    method = optM,
                    control = list(fnscale = -1),
                    hessian = FALSE)
      }else if(optM == "L-BFGS-B"){
        tmp = optim(c(exp(-6.9), exp(2.3)), log.lik,# (0.001, 10) ~= (exp(-6.9), exp(2.3)) >(exp(-4), exp(2.3))
                    yy = y[i, ], xx = x[i, ],
                    mmu = mu[i], sx = sx, sy = sy,
                    mlphi = mlphi, sdlphi = sdlphi,
                    method = optM,
                    lower = rep(0.00001, 2),  # rep(0.00001, 2)
                    upper = c(0.999, exp(700)), #c(0.999, exp(700))
                    control = list(fnscale = myfscale),
                    hessian = FALSE)
      }

      res[i, ] = c(tmp$par, tmp$value, tmp$convergence)
    }
  }
  res[ix, ] = NA

  ### grasp results
  if(optM != "L-BFGS-B"){
    if(trans == "sin"){
      return(data.frame(phi = (sin(res[, 1]) + 1)/2,
                        theta = exp(res[, 2]),
                        obj = res[, 3],
                        convergence = res[, 4]))
    }else if(trans == "exp"){
      return(data.frame(phi = exp(-exp(res[, 1])),
                        theta = exp(res[, 2]),
                        obj = res[, 3],
                        convergence = res[, 4]))
    }
  }else if(optM == "L-BFGS-B"){
    return(data.frame(phi = res[, 1],
                      theta = res[, 2],
                      obj = res[, 3],
                      convergence = res[, 4]))
  }


}

EstiPhi.SingleUpdate <- function(y, x, sy, sx, mu, theta, mlphi, sdlphi){
  require(stats)
  ### only update phi while theta is fixed as its moment estimate
  log.lik <- function(phi, yy, xx, sy, sx, mmu, ttheta, mlphi, sdlphi){
    loglik = sum(dnbinom(yy, size = mmu*(phi^{-1} - 1),
                         prob = 1/(1+sy*ttheta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(phi^{-1} - 1),
                           prob = 1/(1+sx*ttheta), log = TRUE)
    ) + dlnorm(phi, mean = mlphi, sd = sdlphi,
               log = TRUE)
    loglik
  }

  theta[is.na(theta)] = mean(theta, na.rm = TRUE)
  res = matrix(0, nrow = nrow(y), ncol = 3)
  ix = which(is.na(mu))
  for (i in 1:length(mu)) {
    cat(i, sep = "\n")
    if(!is.na(mu[i])){
      tmp = optimize(log.lik, interval = c(0, 1),
                     yy = y[i, ], xx = x[i, ],
                     sy = sy, sx = sx,
                     mmu = mu[i], ttheta = theta[i],
                     mlphi = mlphi, sdlphi = sdlphi,
                     maximum = TRUE, tol = 1e-05)
      res[i, ] = c(tmp$maximum, theta[i], tmp$objective)
    }
    res[ix, ] = NA
  }
  colnames(res) = c("phi", "theta", "obj")
  res = as.data.frame(res)
  res$convergence = 0
  return(res)
}

EstiPhi.IterativeUpdate <- function(y, x, sy, sx, mu, theta0, mlphi, sdlphi,
                                    maxit = 100, eps = 1e-03){
  ### iteratively update phi and theta
  require(stats)
  log.lik.phi <- function(phi, yy, xx, sy, sx, mmu, ttheta, mlphi, sdlphi){
    loglik = sum(dnbinom(yy, size = mmu*(phi^{-1} - 1),
                         prob = 1/(1+sy*ttheta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(phi^{-1} - 1),
                           prob = 1/(1+sx*ttheta), log = TRUE)
    ) + dlnorm(phi, mean = mlphi, sd = sdlphi,
               log = TRUE)
    loglik
  }

  log.lik.theta <- function(theta, yy, xx, sy, sx, mmu, pphi, mlphi, sdlphi){
    loglik = sum(dnbinom(yy, size = mmu*(pphi^{-1} - 1),
                         prob = 1/(1+sy*theta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(pphi^{-1} - 1),
                           prob = 1/(1+sx*theta), log = TRUE)
    ) + dlnorm(pphi, mean = mlphi, sd = sdlphi,
               log = TRUE)
    loglik
  }

  res = matrix(0, nrow = nrow(y), ncol = 4)
  ix = which(is.na(mu))
  for (i in 1:length(mu)) {
    cat(i, sep = "\n")
    if(!is.na(mu[i])){
      delta = 1
      theta.old = theta0
      iter = 0
      while(delta > eps && iter <= maxit){
        iter = iter + 1
        ### update phi while fixing theta
        tmp = optimize(log.lik.phi, interval = c(0, 1),
                       yy = y[i, ], xx = x[i, ],
                       sy = sy, sx = sx,
                       mmu = mu[i], ttheta = theta.old,
                       mlphi = mlphi, sdlphi = sdlphi,
                       maximum = TRUE, tol = 1e-05)
        phi.old = tmp$maximum
        obj.old = tmp$objective
        tmp = optimize(log.lik.theta, interval = c(0, exp(700)),
                       yy = y[i, ], xx = x[i, ],
                       sy = sy, sx = sx,
                       mmu = mu[i], pphi = phi.old,
                       mlphi = mlphi, sdlphi = sdlphi,
                       maximum = TRUE, tol = 1e-05)
        theta.old = tmp$maximum
        obj.new = tmp$objective
        delta = obj.new - obj.old
      }
      res[i, ] = c(phi.old, theta.old, obj.new, as.numeric(iter >= maxit))
    }
    res[ix, ] = NA
  }
  colnames(res) = c("phi", "theta", "obj", "convergence")
  res = as.data.frame(res)
  return(res)
}

# EstiVarofMu <- function(counts, sf, mu, phi, theta, method = c("lnnb", "lnbb")){
#
#   method <- match.arg(method)
#   switch(method,
#          lnnb = EstiVarofMu.NB(counts, sf, mu, phi, theta),
#          lnbb = EstiVarofMu.BB(counts, mu, phi, theta))
# }
# EstiVarofMu.NB <- function(counts, sf, mu, phi, theta){
#   ### variance estimate of mu.hat
#   n = ncol(counts)/2
#   if(is.null(sf)){
#     sf = colSums(counts)/median(colSums(counts))
#   }
#   sf.x = sf[seq(1, ncol(counts), 2)]
#   sf.y = sf[seq(2, ncol(counts), 2)]
#
#   # if(!useTylor){
#   tmp = sum(sf.y^{-1}) + n*theta
#   mu.var = (mu/(n^2*theta))*(phi/(1 -  phi))*tmp
#   # }else if(useTylor){
#   #   ## approximation 2: by taylor's expansion
#   #   tmp = (mu*sum(sf.x^{-1}) + (1-mu)*sum(sf.y^{-1}) + n*theta)
#   #   mu.var = ((mu*(1-mu))/(n^2*theta))*tmp*(phi/(1 -  phi))
#   # }
#
#
#   return(mu.var)
# }
# EstiVarofMu.BB <- function(counts, mu, phi, theta){
#   ### variance estimate of mu.hat, approximating by beta-binomial
#   n = ncol(counts)/2
#   mu.var = n^{-1}*mu*(1- mu)*(1 + theta^{-1})*phi
#   return(mu.var)
# }

EstiVarofMu <- function(counts, sf, mu, phi, theta){
  require(stats)
  ### variance estimate of mu.hat
  n = ncol(counts)/2
  if(is.null(sf)){
    sf = colSums(counts)/median(colSums(counts))
  }
  sf.x = sf[seq(1, ncol(counts), 2)]
  sf.y = sf[seq(2, ncol(counts), 2)]

  tmp = sum(sf.y^{-1}) + n*theta
  mu.var = (mu/(n^2*theta))*(phi/(1 -  phi))*tmp
  return(mu.var)
}



getStats <- function(mu, mu0, var.mu){
  ### formulate statistics
  require(stats)
  stats = (mu - mu0)/sqrt(var.mu)
  return(stats)
}

EstiCutoff <- function(counts, sf, CountThreshold = 100){
  ### estimate cutoff for mu.hat
  require(stats)
  counts.norm = sweep(counts, 2, sf, FUN="/") + 0.5
  ss = colMeans(counts.norm)
  obsRatios = rowSums(counts.norm[,seq(2, length(ss), by=2),drop=FALSE]) / rowSums(counts.norm)

  ix = rowMeans(counts.norm) > CountThreshold ## pick some sites with large counts, because the ratio estimation is more accurate for those sites.
  tmp = mix.2norm(obsRatios[ix], var.eq=FALSE)
  cutoff = (tmp$mu1+tmp$mu2) / 2
  return(cutoff)
}

